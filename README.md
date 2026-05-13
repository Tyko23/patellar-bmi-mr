# patellar-bmi-mr

Reproducibility repository for the Mendelian randomization analysis of body mass index (BMI) and recurrent patellar dislocation reported in:

> Kody MT, Ornelas L, Banffy M, Limpisvasti O. **Genetically Predicted Body Mass Index Causally Increases the Risk of Recurrent Patellar Dislocation: A Mendelian Randomization Study.** Manuscript submitted to *Knee Surgery, Sports Traumatology, Arthroscopy* (2026).

## Project summary

Two-sample summary-level Mendelian randomization testing whether genetically predicted BMI causally increases the risk of recurrent patellar dislocation. The exposure is BMI (Locke et al. 2015 GIANT consortium GWAS, n = 339,224). The outcome is the FinnGen R12 endpoint `M13_PATELLADISLOCATIO` (2,346 cases and 517,864 controls). A sensitivity replication uses the larger Yengo et al. 2018 BMI GWAS.

## Repository contents

| File | Description |
|---|---|
| `analysis.R` | Complete reproducible R script that produces `patellar_paper_results.rds` and `steiger_patellar.rds` from publicly available data sources. |
| `patellar_paper_results.rds` | Saved analysis output: harmonised data, all five univariable MR estimates, sensitivity diagnostics (Cochran's Q, MR-Egger intercept, MR-PRESSO with 10,000 simulations, leave-one-out, single-SNP), I²GX statistic, instrument F-statistic and R². |
| `steiger_patellar.rds` | Steiger directionality test output. |
| `sessionInfo.txt` | R session and package version information used for the published analysis. |
| `LICENSE` | MIT License. |
| `README.md` | This file. |

## Prerequisites

- **R** version ≥ 4.3 (analysis was conducted on R 4.5.3).
- **OpenGWAS API token** for accessing the GIANT BMI summary statistics via the `ieugwasr` package. Sign up at <https://api.opengwas.io/> and set the token as an environment variable: `Sys.setenv(OPENGWAS_JWT = "your_token_here")`.
- **FinnGen R12 summary statistics** for the `M13_PATELLADISLOCATIO` endpoint. Available from the FinnGen public release at <https://www.finngen.fi/en/access_results>. The file used in this analysis is named `summary_stats_release_finngen_R12_M13_PATELLADISLOCATIO.gz`.

### R packages

```r
install.packages(c("data.table", "remotes"))
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages(c("MendelianRandomization", "MRPRESSO"))
remotes::install_github("MRCIEU/ieugwasr")
```

## How to reproduce

1. Set your OpenGWAS API token as described above.
2. Download the FinnGen R12 summary statistics file for `M13_PATELLADISLOCATIO` from the FinnGen public release and place it in a known directory.
3. Open `analysis.R` and update the `FINNGEN_FILE` variable to point to the downloaded file.
4. Run `analysis.R`. The script will:
   - Extract 76 BMI instruments from Locke 2015 via OpenGWAS
   - Filter and harmonise against the FinnGen outcome
   - Run all five univariable MR methods
   - Run the full sensitivity battery
   - Run the Steiger directionality test
   - Save the results to `patellar_paper_results.rds` and `steiger_patellar.rds`

Typical runtime: 10–15 minutes (the MR-PRESSO step with 10,000 simulations is the slowest).

## Citation

If you use code or data from this repository, please cite the manuscript above and the original source GWAS:

- Locke AE, et al. Genetic studies of body mass index yield new insights for obesity biology. *Nature*. 2015;518:197–206.
- Kurki MI, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature*. 2023;613:508–518.

## License

This repository is released under the MIT License (see `LICENSE`). The underlying genome-wide association study summary statistics are subject to the licensing terms of the GIANT consortium and FinnGen.

## Contact

For questions about the analysis, contact Michael T. Kody (michael.kody@gmail.com).
