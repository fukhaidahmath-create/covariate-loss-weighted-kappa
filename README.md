# Analysis inputs

All analysis files have the columns `X` (covariate), `D` (reference), `T1` (test 1) and `T2` (test 2), coded 0/1 for the binary variables. Extra columns (such as `source_row`) are ignored by the R program. Rows are in the order used by the analyses, so the cross-fitting folds are reproduced with seed 20261017.

## 1. Cervical screening (Section 6): `cervical_858_age_biopsy_hinselmann_schiller.csv`

| Column | Meaning |
|---|---|
| `X` | Age in years (integer, 13-84) |
| `D` | Biopsy result, recorded reference (1 = positive) |
| `T1` | Hinselmann test (1 = positive) |
| `T2` | Schiller test (1 = positive) |

Counts: 55 biopsy positive, 35 Hinselmann positive, 74 Schiller positive, 41 discordant test pairs; no missing values.

**Source and license.** Derived from: Fernandes, K., Cardoso, J. S. and Fernandes, J. (2017). *Cervical Cancer (Risk Factors)* [Data set]. UCI Machine Learning Repository. https://doi.org/10.24432/C5Z310. The original dataset is licensed under Creative Commons Attribution 4.0 International (CC BY 4.0); this derivative keeps four of its variables unchanged and is redistributed under the same license with attribution. The data were collected at Hospital Universitario de Caracas, Venezuela; the dataset documentation does not describe the biopsy verification protocol.

The original analysis read the same four columns from the workbook `ccCervical_Cancer_858_Paired_Kappa_RStudio.xlsx` (sheet `RStudio_Input`, MD5 eccd75c3bec8928d428dfcf142cbe1e8). The R pipeline accepts either the workbook or this CSV file.

## 2. Coronary artery disease (Section 7): `heart_cleveland_301_age_cad_exercise_ecg_thallium.csv`

| Column | Meaning |
|---|---|
| `X` | Age in years (integer, 29-77) |
| `D` | Angiographic coronary artery disease: more than 50% diameter narrowing in at least one major vessel (UCI `num` >= 1) |
| `T1` | Exercise electrocardiogram: exercise-induced ST depression >= 1 mm (UCI `oldpeak` >= 1) |
| `T2` | Thallium scintigraphy: fixed or reversible defect (UCI `thal` = 6 or 7) |
| `source_row` | Row number in the 303-record UCI file |

Counts: 301 records (2 of 303 excluded for missing `thal`), 138 with CAD, 136 exercise-ECG positive, 135 thallium positive, 103 discordant test pairs.

`heart_cleveland_301_extended.csv` contains, for the same records, the stricter exercise-ECG definition used in Table S8 (`T1_flat_down`: ST depression >= 1 mm with a flat or downsloping peak-exercise ST segment), sex, `oldpeak`, `slope`, `thal` and `num`.

**Source and license.** Derived from: Janosi, A., Steinbrunn, W., Pfisterer, M. and Detrano, R. (1989). *Heart Disease* [Data set]. UCI Machine Learning Repository. https://doi.org/10.24432/C52P4X (Cleveland Clinic Foundation series; Detrano et al., 1989, *American Journal of Cardiology* 64:304-310). Licensed under CC BY 4.0; this derivative is redistributed under the same license with attribution.

**Derivation.** `R/JBS_V18_PREPARE_HEART_DATA.R` builds both files from the UCI file `processed.cleveland.data` or from `heart_disease_cl.rda` of the R data package `ucidata` (same records); either source gives files identical to those in this folder. The records were also checked against the independently distributed copy `heart.rda` of the R package `kmed`: the 297 common complete records were identical.

## 3. Colorectal cancer with a categorical covariate (Section 7.7): `crc_diabetes_168_fit_fobt_biopsy.csv`

| Column | Meaning |
|---|---|
| `X` | Diabetes (1 = diabetic, 0 = non-diabetic) |
| `D` | Biopsy (1 = positive) |
| `T1` | Faecal immunochemical test, FIT (1 = positive) |
| `T2` | Faecal occult blood test, FOBT (1 = positive) |

The 168 records are fully determined by the counts in Table 9 of Roldán-Nofuentes, J. A. and Regad, S. B. (2025), Simultaneous comparison of the predictive values of two binary diagnostic tests in the presence of categorical covariates, *Journal of Biopharmaceutical Statistics*, https://doi.org/10.1080/10543406.2025.2547589. The categorical example is implemented by the R code `R/JBS_V18_SECOND_APPLICATION_AND_SENSITIVITY.R` in part `crc`, which reconstructs the 168 expanded rows directly from published aggregate cell counts.
