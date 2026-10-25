################################################################################
# JBS_V18_PREPARE_HEART_DATA.R
#
# Builds the two coronary analysis files used in Section 7 of the manuscript
#   data/heart_cleveland_301_age_cad_exercise_ecg_thallium.csv  (X, D, T1, T2, source_row)
#   data/heart_cleveland_301_extended.csv                       (adds T1_flat_down and source variables)
# from the UCI Heart Disease data, Cleveland Clinic series (303 records;
# Janosi, Steinbrunn, Pfisterer and Detrano, 1989; https://doi.org/10.24432/C52P4X;
# CC BY 4.0).
#
# Author : Ismail Radi Ismail Fkhadah, Department of Biostatistics, School of
#          Medicine, University of Granada, Spain (ismailfk2021@correo.ugr.es)
# License: MIT
#
# HOW TO RUN (RStudio). Set the working directory to the repository root and
# point JBS_HEART_SOURCE to ONE of the following files:
#   (a) processed.cleveland.data, downloaded from the UCI repository (comma
#       separated, no header, missing values coded "?"), or
#   (b) heart_disease_cl.rda from the R data package 'ucidata'
#       (https://github.com/coatless-rpkg/ucidata), which contains the same records.
#   Sys.setenv(JBS_HEART_SOURCE = "path/to/processed.cleveland.data")
#   Sys.setenv(JBS_HEART_OUT = "new_data")          # output folder (default "new_data")
#   source("R/JBS_V18_PREPARE_HEART_DATA.R", encoding = "UTF-8")
# The files written are identical to those in data/.
#
# Coding (Section 7.1): X = age (years); D = 1 if num >= 1 (more than 50% diameter
# narrowing in at least one major vessel at angiography); T1 = 1 if oldpeak >= 1
# (exercise-induced ST depression of at least 1 mm); T2 = 1 if thal = 6 or 7 (fixed or
# reversible thallium defect). The two records with missing thal are excluded.
################################################################################

src <- Sys.getenv("JBS_HEART_SOURCE", "processed.cleveland.data")
out <- Sys.getenv("JBS_HEART_OUT", "new_data")
if (!file.exists(src)) stop("Source file not found: ", src, ". Set JBS_HEART_SOURCE (see the header).")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
vars <- c("age", "sex", "cp", "trestbps", "chol", "fbs", "restecg", "thalach",
          "exang", "oldpeak", "slope", "ca", "thal", "num")

if (grepl("\\.rda$", src, ignore.case = TRUE)) {
  e <- new.env(); load(src, envir = e)
  u <- get("heart_disease_cl", envir = e)
  # factor labels used by 'ucidata' -> original UCI codes
  u$sex   <- as.numeric(u$sex == "Male")
  u$slope <- c(upsloping = 1, flat = 2, downsloping = 3)[as.character(u$slope)]
  u$thal  <- c("normal" = 3, "fixed defect" = 6, "reversable defect" = 7)[as.character(u$thal)]
  u <- as.data.frame(u)[, vars]
} else {
  u <- utils::read.csv(src, header = FALSE, col.names = vars, na.strings = "?")
}
stopifnot(nrow(u) == 303L)
cat("num:", paste(names(table(u$num)), table(u$num), sep = "=", collapse = " "),
    "| men:", sum(u$sex), "| mean age:", sprintf("%.2f", mean(u$age)), "\n")

u$source_row <- seq_len(nrow(u))
a <- u[!is.na(u$thal), ]
a$X  <- as.integer(a$age)
a$D  <- as.integer(a$num >= 1)
a$T1 <- as.integer(a$oldpeak >= 1)
a$T2 <- as.integer(a$thal %in% c(6, 7))
a$T1_flat_down <- as.integer(a$oldpeak >= 1 & a$slope %in% c(2, 3))
a$sex_male <- as.integer(a$sex)
a$slope <- c("upsloping", "flat", "downsloping")[a$slope]
a$thal  <- c("3" = "normal", "6" = "fixed defect", "7" = "reversable defect")[as.character(a$thal)]
cat("analysis records:", nrow(a), "| CAD:", sum(a$D), "| exercise ECG positive:", sum(a$T1),
    "| thallium positive:", sum(a$T2), "\n")

num_fmt <- function(v) { v <- as.numeric(v); ifelse(v == round(v), sprintf("%.1f", v), as.character(v)) }
f1 <- file.path(out, "heart_cleveland_301_age_cad_exercise_ecg_thallium.csv")
utils::write.csv(a[, c("X", "D", "T1", "T2", "source_row")], f1, row.names = FALSE, quote = FALSE)
ext <- a[, c("X", "D", "T1", "T2", "T1_flat_down", "sex_male", "oldpeak", "slope", "thal", "num", "source_row")]
ext$oldpeak <- num_fmt(ext$oldpeak); ext$num <- num_fmt(ext$num)
f2 <- file.path(out, "heart_cleveland_301_extended.csv")
utils::write.csv(ext, f2, row.names = FALSE, quote = FALSE)
cat("written:", f1, "and", f2, "\n")
