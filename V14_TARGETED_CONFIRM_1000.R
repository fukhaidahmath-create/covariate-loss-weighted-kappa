# V14 - targeted boundary-null confirmation, preserving all previous results.
# Runs 1,000 simulations PER DIRECTION at n=1500, rho=0.25;
# includes all three prespecified smoothing parameters and both test directions.
# The first 200 repetitions reproduce earlier seeds; this gives 1,000 TOTAL
# repetitions per direction, not 1,000 additional independent repetitions.
# This script is a launcher for the existing V14_Boundary_Null_Simulation.R.

code_file <- normalizePath(file.choose(), winslash="/", mustWork=TRUE)
if (basename(code_file) != "JBS_V14_Boundary_Null_Simulation.R") {
  stop("Select JBS_V14_Boundary_Null_Simulation.R, not the Excel file or V13 code.")
}

# V16 note: the original run used C:/Users/user/Downloads/JBS_V14_CONFIRM_1000_N1500_RHO025
out_dir <- Sys.getenv("JBS_V14_CONFIRM_DIR", unset = file.path(getwd(), "JBS_V14_CONFIRM_1000_N1500_RHO025"))
if (dir.exists(out_dir)) {
  stop("The confirmation output directory already exists. Preserve it; do not overwrite. ",
       "Inspect previous run before choosing a NEW directory.")
}

Sys.setenv(
  JBS_V14_ESTIMATOR="crossfit",
  JBS_V14_REPS="1000",
  JBS_V14_NS="1500",
  JBS_V14_RHOS="0.25",
  JBS_V14_TAUS="0.0001,0.0005,0.002",
  JBS_V14_DIR=out_dir
)
cat("Starting targeted independent boundary validation.\n",
    "Output directory:",out_dir,"\n")
source(code_file, encoding="UTF-8")

raw_files <- list.files(out_dir, pattern="^V14_RAW_.*_crossfit\\.csv$", full.names=TRUE)
if (length(raw_files) != 2L) {
  stop("Expected two completed boundary-direction RAW files; found ",length(raw_files))
}
a <- do.call(rbind, lapply(raw_files, utils::read.csv, stringsAsFactors=FALSE))
stopifnot(nrow(a)==6000L, all(a$success), length(unique(a$tau))==3L)
a$any_cert <- as.integer(a$T1_certificate==1L | a$T2_certificate==1L)
family <- aggregate(any_cert ~ direction + n + rho + rep, data=a, FUN=max)
summary <- aggregate(any_cert ~ direction + n + rho, data=family,
                     FUN=function(x)c(errors=sum(x),total=length(x),rate=mean(x)))
result <- data.frame(direction=summary$direction,n=summary$n,rho=summary$rho,
                     errors=summary$any_cert[,"errors"],
                     replications=summary$any_cert[,"total"],
                     familywise_error_rate=summary$any_cert[,"rate"])
utils::write.csv(result,file.path(out_dir,"V14_TARGETED_FWER_1000.csv"),row.names=FALSE)
cat("\nFamilywise error across all three taus AND both directions:\n")
print(result, row.names=FALSE)
cat("\nFINISHED; upload V14_TARGETED_FWER_1000.csv and V14_BOUNDARY_SIZE_SUMMARY_crossfit.csv.\n")
