# Covariate-Dependent Loss-Weighted Kappa for Comparing Paired Binary Diagnostic Tests

**R source code and reproducibility materials | Public software release v1.3.0**

**Author:** Ismail Radi Ismail Fkhadah  
**Affiliation:** University of Granada, Granada, Spain  
**ORCID:** [0009-0004-4092-7836](https://orcid.org/0009-0004-4092-7836)  
**Archived, fixed release:** [Zenodo DOI 10.5281/zenodo.23628904](https://doi.org/10.5281/zenodo.23628904)

## Overview

This repository provides the original R source code and compact prepared analytical datasets associated with the research manuscript *Covariate-Dependent Loss-Weighted Kappa for Comparing Paired Binary Diagnostic Tests: Orthogonal Estimation and Policy-Robust Dominance*.

The methods concern covariate-dependent misclassification-loss-weighted kappa for paired binary diagnostic tests. The deposited programs cover cross-fitted estimation, influence-function-based comparisons, finite-grid and boundary investigations, sensitivity analyses, secondary diagnostic applications, and figure generation, subject to the conditions and limitations in the research materials.

**This GitHub repository contains readable source code and compact input files. The fixed software release, archived simulation outputs, tables, figures, and historical execution records are preserved on [Zenodo](https://doi.org/10.5281/zenodo.23628904).**

## R source files

All seven R programs currently reside in the **repository root** (not in a separate `R/` directory):

| Program | Purpose |
|---|---|
| [`JBS_V13_3_ONE_CLICK_ALL_STAGES.R`](JBS_V13_3_ONE_CLICK_ALL_STAGES.R) | Main analysis and simulation pipeline |
| [`JBS_V14_Boundary_Null_Simulation.R`](JBS_V14_Boundary_Null_Simulation.R) | Boundary-null simulation experiments |
| [`V14_TARGETED_CONFIRM_1000.R`](V14_TARGETED_CONFIRM_1000.R) | Targeted confirmation experiments |
| [`JBS_V15_NUISANCE_RATE_DIAGNOSTIC.R`](JBS_V15_NUISANCE_RATE_DIAGNOSTIC.R) | Nuisance-estimation rate diagnostics |
| [`JBS_V18_PREPARE_HEART_DATA.R`](JBS_V18_PREPARE_HEART_DATA.R) | Heart disease data preparation |
| [`JBS_V18_SECOND_APPLICATION_AND_SENSITIVITY.R`](JBS_V18_SECOND_APPLICATION_AND_SENSITIVITY.R) | Additional applications and sensitivity analyses |
| [`JBS_V18_FIGURES_ALL.R`](JBS_V18_FIGURES_ALL.R) | Figure generation |

The version identifiers in individual R filenames refer to historical development stages. **The first public software release is v1.3.0.** The R programs are preserved unchanged from the archived release.

## Analytical inputs

The four prepared CSV files are also in the repository root:

- [`cervical_858_age_biopsy_hinselmann_schiller.csv`](cervical_858_age_biopsy_hinselmann_schiller.csv): cervical screening, 858 records.
- [`heart_cleveland_301_age_cad_exercise_ecg_thallium.csv`](heart_cleveland_301_age_cad_exercise_ecg_thallium.csv): Cleveland heart disease, 301 complete analytical records.
- [`heart_cleveland_301_extended.csv`](heart_cleveland_301_extended.csv): additional heart disease variables and alternative indicator construction.
- [`crc_diabetes_168_fit_fobt_biopsy.csv`](crc_diabetes_168_fit_fobt_biopsy.csv): 168-row reconstruction from published aggregate counts, not individual patient records.

See [DATA_PROVENANCE.md](DATA_PROVENANCE.md) for detailed variable definitions, upstream sources, data transformations, counts, and attribution. **The MIT License for original R code does not supersede third-party data licensing.**

## Getting started in RStudio

Download the repository using **Code → Download ZIP**, extract the files, open RStudio, and set the working directory to the extracted repository folder. The main analysis program uses the historical environment variable `JBS_EXCEL_PATH`, which can also point to a CSV file. An example for cervical screening is:

```r
install.packages(c("readxl", "nnet", "brglm2"))  # if needed
Sys.setenv(
  JBS_ONE_CLICK = "0",
  JBS_STAGE = "analysis",
  JBS_EXCEL_PATH = normalizePath(
    "cervical_858_age_biopsy_hinselmann_schiller.csv",
    mustWork = TRUE
  ),
  JBS_OUTPUT_DIR = file.path(getwd(), "new_results_cervical"),
  JBS_BOOTSTRAP_B = "999",
  JBS_MULTIPLIER_B = "1999"
)
source("JBS_V13_3_ONE_CLICK_ALL_STAGES.R", encoding = "UTF-8")
```

Consult the individual scripts and the documentation in the [archived software release](https://doi.org/10.5281/zenodo.23628904) for their detailed parameters, dependencies, computational requirements, and limitations. Do not overwrite the archived results when rerunning analyses.

## Archived results and reproducibility scope

Zenodo includes five separate archives for (1) R code, (2) analytical data, (3) archived R results, tables, figures, simulations and logs, (4) documentation, and (5) supplementary numerical outputs from a 324-scenario prevalence-by-loss-weight experiment, with 400 recorded replications per scenario.

**The supplementary 324-scenario experiment is distributed as archived numerical results only. An executable R implementation for regenerating that specific experiment is not included in this release.** A saved report documents the planned 120-cell main simulation grid. The presence of archived outputs does not assert an independently completed rerun of every historical analysis.

## Citation

Please cite the fixed Zenodo software version:

> Fkhadah, I. R. I. (2026). *Covariate-Dependent Loss-Weighted Kappa for Comparing Paired Binary Diagnostic Tests: R Code and Reproducibility Materials* (Version v1.3.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23628904

See [`CITATION.cff`](CITATION.cff) for the machine-readable citation metadata. The manuscript is a separate research work from this software release.

## License

Original R code: [MIT License](LICENSE). UCI-derived datasets: CC BY 4.0 with attribution. The categorical-covariate example is reconstructed from published aggregate counts; its source and reuse conditions are described in [DATA_PROVENANCE.md](DATA_PROVENANCE.md). Clinical loss weights are illustrative and the software is intended for methodological research.
