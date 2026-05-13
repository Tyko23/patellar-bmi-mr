# =====================================================================
# Mendelian randomization analysis: BMI -> recurrent patellar dislocation
#
# Companion code for:
#   Kody MT, Ornelas L, Banffy M, Limpisvasti O.
#   Genetically Predicted Body Mass Index Causally Increases the Risk of
#   Recurrent Patellar Dislocation: A Mendelian Randomization Study.
#   Submitted to Knee Surg Sports Traumatol Arthrosc (2026).
#
# Run from a clean R session. Estimated runtime: 10-15 minutes.
#
# Prerequisites (see README.md):
#   - R >= 4.3
#   - OpenGWAS API token set as environment variable OPENGWAS_JWT
#   - FinnGen R12 M13_PATELLADISLOCATIO summary statistics file
#   - Packages: TwoSampleMR, MendelianRandomization, MRPRESSO,
#               ieugwasr, data.table
# =====================================================================

# ---- User configuration --------------------------------------------------
# Point this to the FinnGen R12 summary statistics file for
# M13_PATELLADISLOCATIO. The file is publicly available from
# https://www.finngen.fi/en/access_results
FINNGEN_FILE <- "~/Desktop/summary_stats_release_finngen_R12_M13_PATELLADISLOCATIO.gz"

# Outcome metadata (from FinnGen R12 endpoint definition)
N_CASE <- 2346
N_CTRL <- 517864
N_TOTAL <- N_CASE + N_CTRL
PREVALENCE <- N_CASE / N_TOTAL

# Exposure sample size (Locke 2015)
N_EXPOSURE <- 339224

# Set the OpenGWAS token (uncomment and paste yours, or set in .Renviron):
# Sys.setenv(OPENGWAS_JWT = "your_token_here")

# ---- Libraries -----------------------------------------------------------
library(TwoSampleMR)
library(MendelianRandomization)
library(MRPRESSO)
library(data.table)

# Verify API status
stopifnot(!identical(Sys.getenv("OPENGWAS_JWT"), ""))
ieugwasr::api_status()

# ---- 1. Extract BMI instruments (Locke 2015 GIANT, OpenGWAS ieu-a-2) -----
set.seed(42)
exp_bmi <- extract_instruments(
  outcomes = "ieu-a-2",
  clump    = TRUE,
  r2       = 0.001,
  kb       = 10000,
  p1       = 5e-8
)
exp_bmi$exposure    <- "Body mass index"
exp_bmi$id.exposure <- "bmi_locke"
cat("BMI instruments after clumping:", nrow(exp_bmi), "\n")  # 78

# ---- 2. Read FinnGen outcome (memory-efficient) -------------------------
raw <- fread(
  FINNGEN_FILE,
  select = c("rsids", "alt", "ref", "beta", "sebeta", "af_alt", "pval")
)
raw <- as.data.frame(raw[rsids %in% exp_bmi$SNP])
cat("SNPs in FinnGen outcome with BMI-instrument SNPs:", nrow(raw), "\n")

out <- format_data(
  raw, type = "outcome",
  snp_col           = "rsids",
  beta_col          = "beta",
  se_col            = "sebeta",
  effect_allele_col = "alt",
  other_allele_col  = "ref",
  eaf_col           = "af_alt",
  pval_col          = "pval"
)
out$outcome    <- "Patellar dislocation (recurrent)"
out$id.outcome <- "patelladisloc_finngen_r12"

# ---- 3. Harmonise (action = 2: drop ambiguous palindromic SNPs) ---------
dat <- harmonise_data(exp_bmi, out, action = 2)
cat("SNPs retained for analysis (mr_keep = TRUE):", sum(dat$mr_keep), "\n")  # 76

# ---- 4. Primary MR: five complementary methods --------------------------
res <- mr(dat)
print(res[, c("method", "nsnp", "b", "se", "pval")])

# ---- 5. Sensitivity battery ---------------------------------------------
het <- mr_heterogeneity(dat)            # Cochran's Q (IVW) and Rucker's Q' (Egger)
ple <- mr_pleiotropy_test(dat)          # MR-Egger intercept test
loo <- mr_leaveoneout(dat)              # leave-one-out
sin <- mr_singlesnp(dat)                # single-SNP Wald ratios (also drives funnel plot)

# I^2_GX statistic (Bowden 2016) for MR-Egger reliability
isq <- MendelianRandomization::Isq(
  dat$beta.exposure[dat$mr_keep],
  dat$se.exposure[dat$mr_keep]
)
cat("I^2_GX =", round(isq, 3), "\n")  # 0.98

# MR-PRESSO with 10,000 simulations
set.seed(42)
presso <- mr_presso(
  BetaOutcome     = "beta.outcome",
  BetaExposure    = "beta.exposure",
  SdOutcome       = "se.outcome",
  SdExposure      = "se.exposure",
  OUTLIERtest     = TRUE,
  DISTORTIONtest  = TRUE,
  data            = dat[dat$mr_keep, ],
  NbDistribution  = 10000,
  SignifThreshold = 0.05
)

# ---- 6. Instrument strength --------------------------------------------
F_mean <- mean((exp_bmi$beta.exposure / exp_bmi$se.exposure)^2)
R2     <- sum(2 * exp_bmi$eaf.exposure * (1 - exp_bmi$eaf.exposure) * exp_bmi$beta.exposure^2)
cat("Mean F-statistic:", round(F_mean, 2), "\n")  # 66.29
cat("Cumulative R^2: ", round(R2 * 100, 2), "%\n")  # 2.35%

# ---- 7. Save primary results -------------------------------------------
saveRDS(
  list(
    harmonised    = dat,
    univariable   = res,
    heterogeneity = het,
    pleiotropy    = ple,
    presso        = presso,
    isq           = isq,
    loo           = loo,
    singlesnp     = sin,
    F_mean        = F_mean,
    R2            = R2,
    exposure_data = exp_bmi
  ),
  "patellar_paper_results.rds"
)
cat("Saved: patellar_paper_results.rds\n")

# ---- 8. Steiger directionality test ------------------------------------
# Steiger requires the appropriate columns for a binary outcome.
dat$samplesize.exposure <- N_EXPOSURE
dat$samplesize.outcome  <- N_TOTAL
dat$ncase.outcome       <- N_CASE
dat$ncontrol.outcome    <- N_CTRL
dat$prevalence.outcome  <- PREVALENCE
dat$units.outcome       <- "log odds"
dat$units.exposure      <- "SD"

stei <- directionality_test(dat)
print(stei)
saveRDS(stei, "steiger_patellar.rds")
cat("Saved: steiger_patellar.rds\n")

# ---- 9. Capture R session info -----------------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("Saved: sessionInfo.txt\n")

cat("\nAll analyses complete.\n")
