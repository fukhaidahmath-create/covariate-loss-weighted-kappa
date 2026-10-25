# JBS V13.2 FULL RESEARCH AND SCIENTIFIC AUDIT PIPELINE
# Primary statistical contribution: continuous covariate, covariate-dependent
# misclassification loss-weighted KAPPA for paired binary diagnostic tests.
# Context-only categorical PPV/NPV comparisons from a DIFFERENT manuscript.
# See README and methodological map before using paper's name or claims.
# Run: Sys.setenv(JBS_STAGE="selftest"); source(file.choose(), encoding="UTF-8")
# Stages: selftest, protocol, analysis, generate_data, sensitivity,
#         mc_pilot, mc_full, mc_sequence, pv_mc_pilot, pv_mc_full,
#         local_calibration, decision, all_pilot, all_full, all_stages.
# Monte Carlo only stages DO NOT require the Excel file.
# WARNING: all_stages can run tens of thousands of simulated datasets.
# Reference paper: Roldan-Nofuentes & Regad (2025), DOI 10.1080/10543406.2025.2547589;
# main predictive-value estimand/methods are DIFFERENT, not reproduced or relabeled.
# DO NOT interpret self-tests or successful MC fits as proof of journal acceptance.

# ==============================================================================
# V12 COMPREHENSIVE RESEARCH PIPELINE | BASE V11.1 REPRODUCIBLE INFERENCE / REPRODUCIBLE TABLES AND FIGURES — COVARIATE-DEPENDENT LOSS-WEIGHTED KAPPA
# Paired binary diagnostic tests | cervical-cancer UCI real data | RStudio
# ==============================================================================
# الترخيص البحثي: كود قابل للمراجعة وإعادة الإنتاج، وليس برنامجاً طبياً معتمداً.
# هذا الملف نصي TXT يحتوي على أوامر R صالحة؛ يمكن فتحه في RStudio ثم Source.
# INPUT: Cervical_Cancer_858_Paired_Kappa_RStudio.xlsx / sheet RStudio_Input
# COLUMNS: X=Age, D=Biopsy, T1=Hinselmann, T2=Schiller.
# المصدر: https://archive.ics.uci.edu/dataset/383/cervical+cancer+risk+factors
# BACKGROUND: V5 local splines, V7 integrated fixed-policy orthogonal estimation
#            and fixed-tau conservative lower-envelope dominance inference.
# IMPORTANT: لا يوجد ضمان بقبول المخطوطة في أي مجلة أو بأن صلاحية R تعني صلاحية الاستدلال.
# IMPORTANT: لا تفترض هذه المنهجية تحقق المرجع التشخيصي بلا خطأ أو اختيار عينات ممثل.
# IMPORTANT: خسائر FN/FP لم تُقَس في البيانات؛ القيم الافتراضية أمثلة تحليل حساسية فقط.
# IMPORTANT: الجزء المحلي على finite grid ولا يدعي uniform continuum coverage.
# IMPORTANT: bootstrap model now uses prespecified ridge penalty and fixed spline basis;
#            passing a 95% computational success screen does NOT certify 95% coverage.
# IMPORTANT: plots/tables support manuscript preparation, NOT acceptance guarantees.
# IMPORTANT: اختبار هيمنة tau ثابت هو pointwise asymptotically conservative؛
#            ليس uniform-in-P ولا إجراءً دقيقاً على الحدود غير المنتظمة.
# IMPORTANT: افترض أن جميع السجلات مستقلة بين الأشخاص ومرجع المرض معلوم لكل شخص.
# V10.1: Newton + line-search + two warm-start fallbacks for rare-event resamples.
# V10.1: same ridge likelihood, penalty, basis, tolerance and inferential policy as V10.
# V10.1: NEVER convert a nonstationary fit to success or alter failed replicates.
# V10: consistent penalized local regression, prespecified lambda, fixed knots.
# V10: transparent bootstrap acceptance/failure summaries, journal-ready tables & figures.
# V10: finite-grid bands remain exploratory until MC calibration.
# V9.1 PATH HOTFIX (17 Oct 2026):
# - Select Excel BEFORE establishing the results directory. This avoids a Windows
#   file chooser / working-directory switch invalidating relative output paths.
# - Store a normalized, absolute workbook path and resolve output next to it.
# - Do not silently ignore dir.create() failures: check existence and write access.
# - Print the resolved Excel and output locations before producing any CSV files.
# - No change to formulas, estimands, influence functions, or inferential settings.
# V9 CHANGELOG (17 Oct 2026):
# - جميع النماذج اللوجستية: bias-reduced adjusted scores (brglm2::brglmFit).
#   هذا يعالج مشكلة complete/quasi separation دون إخفاء التحذيرات.
# - كل محاولة bootstrap لها حالة وسبب فشل أو عدد تحذيرات، حتى عند فشل band.
# - الاعتماد على نسبة النجاح المحددة مسبقاً؛ لا خفض اصطناعي لها.
# - فحص multivariate Delta Method مقابل influence-function SE آلياً.
# - DGM Monte Carlo يستخدم إعداد nuisance مطابقاً للتحاليل الرئيسية.
# - مقارنة مشتركة أخرى باستخدام نموذج multinomial تبقى استطلاعية وحساسة للخلايا النادرة.
# - ليس هناك ضمان أن biased-reduced local estimator has nominal bootstrap coverage:
#   تغطية finite-grid bands تحتاج Monte Carlo فعلي. لم يُجر تشغيل R في بيئة إعداد الملف.
# ==============================================================================

# 00 | إعدادات قابلة للتعديل — افحصها قبل أي تحليل أو نشر.
CFG <- list(
  excel = "Cervical_Cancer_858_Paired_Kappa_RStudio.xlsx",
  sheet = "RStudio_Input",
  out_dir = "JBS_V12_RESULTS",
  seed = 20261017L,
  alpha = 0.05,
  folds = 5L,
  nuisance_model = "spline",    # "spline" or "linear"; predetermination is advised
  nuisance_df = 3L,             # leave unchanged for primary prespecified analysis
  ref_denominator_tol = 1e-8,
  local_df_prev = 3L,
  local_df_pair = 2L,
  local_primary = "marginal",    # bootstrap stable with sparse joint cells
  joint_multinomial_sensitivity = TRUE,
  local_nx = 41L,
  local_c_values = seq(0.1, 0.9, by = 0.1),
  bootstrap_local = TRUE,       # B=199: pilot only; use >=999 for final Monte Carlo review
  bootstrap_B = 199L,
  bootstrap_min_success = 0.95,  # operational screening only, NOT a coverage guarantee
  local_ridge_lambda = 0.8,     # unchanged: same scientific estimator as V10
  ridge_gradient_tol = 1e-3,   # unchanged tolerance; do NOT relax to force success
  logistic_strategy = "brglm2_AS_mean",   # bias reduction prevents infinite coefficients
  logistic_maxit = 200L,
  run_simulation = FALSE,       # TRUE only after tests and research planning
  simulation_reps = 300L,      # aim for >=1000 for meaningful tail-size verification
  simulation_ns = c(250L, 500L, 858L, 1500L),
  simulation_truth_models = c("constant", "varying"),
  run_sample_size_search = FALSE,
  theta_taus = c(0.0001, 0.0005, 0.002),  # fixed BEFORE outcome inspection
  envelope_L = 0.40,          # prespecified box [L(x),U(x)]; examples only
  envelope_U = 0.90,
  age_policy_center = 35,
  age_policy_logit_at_center = qlogis(0.65),
  age_policy_age_slope = 0.025,
  # V11 ADDITIONS: all exploratory/prespecified knobs saved in output script.
  multiplier_B = 1999L,                # >= 1999 recommended for final finite-grid band
  endpoint_interval = c(.10,.90),      # continuum exact sign property for CONSTANT c
  run_v11_simulation = FALSE,          # ENABLE before producing publication MC tables
  v11_simulation_reps = 1000L,         # publication target per MC scenario
  v11_simulation_ns = c(250L, 500L, 858L, 1500L),
  v11_simulation_scenarios = c("null","t1_strong","t2_strong","tradeoff"),
  v11_simulation_shapes = c("linear","nonlinear"),
  # Optional time-intensive calibration of the LOCAL age-by-loss bootstrap band.
  run_v11_local_grid_calibration = FALSE,
  v11_local_mc_reps = 200L,
  v11_local_mc_bootstrap_B = 199L,
  v11_local_mc_n = 858L
)

# 01 | الحزم الضرورية: لا تثبت الحزم تلقائياً لتجنب تغيّر بيئة المشروع.
required_packages <- c("readxl", "splines", "nnet", "brglm2")
not_installed <- required_packages[!vapply(required_packages,
                                           requireNamespace, logical(1), quietly = TRUE)]
if (length(not_installed)) {
  stop("Install required packages once, then source again: install.packages(c(",
       paste(sprintf("'%s'", not_installed), collapse = ", "), "))")
}
options(stringsAsFactors = FALSE, scipen = 6)
set.seed(CFG$seed)

# 02 | أدوات أساسية: clip للانحدار وكتابة CSV صالحة للقراءة والتكرار.
clip01 <- function(p, eps = 1e-7) pmax(eps, pmin(1 - eps, as.numeric(p)))
validate_binary <- function(x, label) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
      any(!(x %in% c(0, 1)))) stop(label, " must contain finite numeric 0/1 only")
  invisible(TRUE)
}
write_csv <- function(x, filename) {
  if (!dir.exists(CFG$out_dir)) {
    stop("Output directory is missing: ", CFG$out_dir,
         ". Check permissions and Excel location.")
  }
  destination <- file.path(CFG$out_dir, filename)
  tryCatch(
    utils::write.csv(x, destination, row.names = FALSE,
                     na = "", fileEncoding = "UTF-8"),
    error = function(e) {
      stop("Cannot write output '" , destination, "': ", conditionMessage(e))
    }
  )
}
logline <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = "")

# 03 | تحقق من ملف Excel الحقيقي، من الورقة الرسمية وليس من ورقة التلخيص.
read_research_data <- function(file = CFG$excel, sheet = CFG$sheet) {
  if (!file.exists(file) && interactive()) {
    message("Select the file Cervical_Cancer_858_Paired_Kappa_RStudio.xlsx")
    file <- file.choose()
  }
  if (!file.exists(file)) stop("Excel file not found: ", file,
                                ". Put TXT and XLSX in the working directory or edit CFG$excel")
  # V16 PATCH (19 Oct 2026; input handling only, estimators unchanged):
  # (a) a CSV copy of the analysis input (repository file data/*.csv) is accepted;
  # (b) extra columns such as 'source_row' no longer stop the run: the four
  #     analysis columns X, D, T1, T2 are selected in this order.
  # This removes the stop recorded in ONE_CLICK_LAST_ERROR.txt
  # ("RStudio_Input must have exactly X,D,T1,T2 in order; got: X, D, T1, T2, source_row").
  dat <- if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    as.data.frame(readxl::read_excel(file, sheet = sheet))
  }
  must <- c("X", "D", "T1", "T2")
  if (!all(must %in% names(dat))) stop("Input must contain columns X, D, T1, T2; got: ",
                                       paste(names(dat), collapse = ", "))
  if (!identical(names(dat), must)) {
    message("Using columns X, D, T1, T2 only; other columns ignored: ",
            paste(setdiff(names(dat), must), collapse = ", "))
    dat <- dat[, must, drop = FALSE]
  }
  for (v in must) dat[[v]] <- as.numeric(dat[[v]])
  if (anyNA(dat) || any(!is.finite(as.matrix(dat)))) {
    stop("RStudio_Input has missing or non-finite values; audit source workbook")
  }
  for (v in c("D", "T1", "T2")) validate_binary(dat[[v]], v)
  if (length(unique(dat$X)) < 8L || stats::sd(dat$X) <= 0)
    stop("Continuous-scale covariate X has inadequate variation")
  if (length(unique(dat$D)) < 2L) stop("D must have both disease states")
  if (nrow(dat) != 858L) warning("Expected 858 UCI records; observed ", nrow(dat),
                               ". Inspect completeness before publication")
  rownames(dat) <- NULL
  # Absolute source path survives later working-directory changes on Windows.
  attr(dat, "source_excel") <- normalizePath(file, winslash = "/", mustWork = TRUE)
  dat
}

# 04 | توثيق بيانات المرض: تقارير السجلات، والارتباط داخل الشخص.
make_qc <- function(dat) {
  n <- nrow(dat)
  vals <- data.frame(
    metric = c("n", "biopsy_positive", "biopsy_negative", "T1_positive", "T2_positive",
               "paired_discordant", "age_min", "age_max", "age_unique",
               "all_complete", "X_mean", "X_sd"),
    value = c(n, sum(dat$D), n - sum(dat$D), sum(dat$T1), sum(dat$T2),
              sum(dat$T1 != dat$T2), min(dat$X), max(dat$X),
              length(unique(dat$X)), as.integer(!anyNA(dat)),
              mean(dat$X), stats::sd(dat$X)))
  cells <- as.data.frame(table(D = dat$D, T1 = dat$T1, T2 = dat$T2))
  names(cells)[4] <- "n"
  diagacc <- do.call(rbind, lapply(c("T1","T2"), function(v) {
    y <- dat[[v]]
    data.frame(test=v, prevalence=mean(dat$D),
               sensitivity=mean(y[dat$D==1]),
               specificity=mean(1-y[dat$D==0]),
               PPV=if (sum(y==1)>0) mean(dat$D[y==1]) else NA_real_,
               NPV=if (sum(y==0)>0) mean(1-dat$D[y==0]) else NA_real_,
               false_negative=sum(dat$D==1 & y==0),
               false_positive=sum(dat$D==0 & y==1))
  }))
  list(summary = vals, paired_cells = cells, accuracy=diagacc)
}

# 05 | ثبت دوال خسارة «معلومة ومفترضة مسبقاً» بدون ملاءمة نتائج المرض.
# s(x) هو حجم الخسارة، وc(x) نسبة خسارة FN إلى إجمالي الخسارتين.
loss_scale <- function(x) rep(1, length(x))
relative_fn_loss <- function(x, name = "age_dependent") {
  if (name == "equal") return(rep(0.5, length(x)))
  if (name == "fn_075") return(rep(0.75, length(x)))
  if (name == "age_dependent")
    return(stats::plogis(CFG$age_policy_logit_at_center +
              CFG$age_policy_age_slope * (x - CFG$age_policy_center)))
  stop("Unknown clinical loss schedule: ", name)
}
box_L <- function(x) rep(CFG$envelope_L, length(x))
box_U <- function(x) rep(CFG$envelope_U, length(x))
check_loss <- function(x, s, cval, L = NULL, U = NULL) {
  if (length(s) != length(x) || length(cval) != length(x) ||
      any(!is.finite(s)) || any(s <= 0) ||
      any(!is.finite(cval)) || any(cval < 0 | cval > 1)) stop("Bad s(X) or c(X)")
  if (!is.null(L) || !is.null(U)) {
    if (length(L) != length(x) || length(U) != length(x) ||
        any(!is.finite(L)) || any(!is.finite(U)) ||
        any(L < 0 | L > U | U > 1)) stop("Invalid bounds: 0 <= L <= U <= 1")
  }
  invisible(TRUE)
}

# 06 | تدريب الانحدارات الثلاثة على جزء التدريب فقط؛ لا نفترض استقلال T1 وT2.
# لا نحتاج نموذجاً مستقلاً للاختبارين للتقدير المتكامل؛ اعتماد الاختبارين يظهر في EIF.
# تابع تقدير لوجستي مصحح لانحياز ML: brglm2 يقلل مشاكل الانفصال.
# Uses mean bias reduction, which can give finite regression coefficients even
# for complete/quasi separation. IMPORTANT: changes finite-sample estimators.
# Captured warnings are recorded in the model, NOT erased from the audit trail.
fit_safe_logistic <- function(f, train, label="unnamed model") {
  if (!identical(CFG$logistic_strategy,"brglm2_AS_mean"))
    stop("Unsupported logistic_strategy")
  caught <- character(0)
  obj <- tryCatch(
    withCallingHandlers(
      stats::glm(formula=f, data=train, family=stats::binomial(link="logit"),
                 method=brglm2::brglmFit, type="AS_mean",
                 control=stats::glm.control(maxit=CFG$logistic_maxit)),
      warning=function(w) {
        caught <<- c(caught,conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ), error=function(e) {
      stop("Bias-reduced GLM failed [",label,"]: ",conditionMessage(e))
    }
  )
  if (!isTRUE(obj$converged) || anyNA(stats::coef(obj)) ||
      any(!is.finite(stats::coef(obj))))
    stop("Bias-reduced GLM did not return a converged finite fit [",label,"]")
  attr(obj,"audit_warnings") <- unique(caught)
  attr(obj,"audit_method") <- "brglm2::brglmFit / AS_mean"
  obj
}

# لا تتغير صيغة الانحدار الأصلي؛ الذي يتغير هو طريقة التقدير.
fit_binary_model <- function(train, outcome, method, df) {
  f <- if (method == "linear") {
    stats::reformulate("X", response=outcome)
  } else if (method == "spline") {
    stats::as.formula(paste0(outcome," ~ splines::ns(X, df = ",df,")"))
  } else stop("method must be linear or spline")
  fit_safe_logistic(f,train,label=paste0(outcome," / ",method))
}

# 07 | تعيين الطيات مستقل عن الاستجابات D,T1,T2؛ لا تسرب لبيانات validation.
crossfit_probabilities <- function(dat, K = CFG$folds,
                                   method = CFG$nuisance_model,
                                   df = CFG$nuisance_df,
                                   seed = CFG$seed) {
  n <- nrow(dat)
  if (K < 2L || K > n / 20) stop("Number of cross-fitting folds is inappropriate")
  set.seed(seed)
  fold <- sample(rep(seq_len(K), length.out = n), size = n, replace = FALSE)
  pred <- matrix(NA_real_, nrow = n, ncol = 3L,
                 dimnames = list(NULL, c("p", "q1", "q2")))
  diag <- list()
  for (k in seq_len(K)) {
    id <- which(fold == k)
    tr <- dat[-id, , drop = FALSE]
    te <- dat[id, , drop = FALSE]
    info <- c(fold = k, train_n = nrow(tr), validation_n = nrow(te))
    for (j in seq_along(c("D", "T1", "T2"))) {
      response <- c("D", "T1", "T2")[j]
      m <- fit_binary_model(tr, response, method, df)
      pred[id, j] <- clip01(stats::predict(m, newdata = te, type = "response"))
      info[paste0(response, "_train_events")] <- sum(tr[[response]])
      info[paste0(response,"_fit_warnings")] <- length(attr(m,"audit_warnings"))
      info[paste0(response,"_fit_warning_details")] <-
        paste(attr(m,"audit_warnings"),collapse=" | ")
    }
    diag[[k]] <- as.data.frame(as.list(info), stringsAsFactors = FALSE)
  }
  if (anyNA(pred) || any(!is.finite(pred))) stop("Invalid out-of-fold predictions")
  dd <- do.call(rbind, diag)
  metrics <- data.frame(outcome = c("D", "T1", "T2"),
                        Brier = c(mean((dat$D-pred[, 1])^2),
                                  mean((dat$T1-pred[, 2])^2),
                                  mean((dat$T2-pred[, 3])^2)),
                        observed = c(mean(dat$D), mean(dat$T1), mean(dat$T2)),
                        predicted = colMeans(pred))
  list(p = pred[, 1], q1 = pred[, 2], q2 = pred[, 3],
       fold = fold, diagnostics = dd, calibration = metrics)
}

# 08 | أساس التقدير المتعامد: A_j = E[a_j] و C_j=E[s Cov(D,Tj |X)].
# phi_Cj=s(X)(D-p)(Tj-qj)-Cj هي دالة التأثير الكفء مع s ثابتة.
fixed_policy_inference <- function(dat, nuisance, cval, s = loss_scale(dat$X),
                                   policy_name = "specified", alpha = CFG$alpha) {
  check_loss(dat$X, s, cval)
  n <- nrow(dat)
  D <- dat$D; T1 <- dat$T1; T2 <- dat$T2
  p <- nuisance$p; q1 <- nuisance$q1; q2 <- nuisance$q2
  a1 <- s*(cval*D*(1-T1) + (1-cval)*(1-D)*T1)
  a2 <- s*(cval*D*(1-T2) + (1-cval)*(1-D)*T2)
  u1 <- s*(D-p)*(T1-q1)
  u2 <- s*(D-p)*(T2-q2)
  A1 <- mean(a1); A2 <- mean(a2); C1 <- mean(u1); C2 <- mean(u2)
  R1 <- A1+C1; R2 <- A2+C2
  if (!is.finite(R1) || !is.finite(R2) ||
      min(R1, R2) <= CFG$ref_denominator_tol)
    stop("Estimated reference loss is non-positive or near zero; no kappa inference")
  k1 <- C1/R1; k2 <- C2/R2; delta <- k1-k2
  phiA1 <- a1-A1; phiA2 <- a2-A2
  phiC1 <- u1-C1; phiC2 <- u2-C2
  phiK1 <- (A1*phiC1-C1*phiA1)/R1^2
  phiK2 <- (A2*phiC2-C2*phiA2)/R2^2
  phiDelta <- phiK1-phiK2
  std_error <- function(phi) stats::sd(phi)/sqrt(length(phi))
  se1 <- std_error(phiK1); se2 <- std_error(phiK2); seD <- std_error(phiDelta)
  # Independent algebraic computation of the SAME multivariate Delta variance.
  # Cov of subject-level components includes T1-T2 paired covariance terms.
  g_Delta <- c(-C1/R1^2, A1/R1^2, C2/R2^2, -A2/R2^2)
  phi_components <- cbind(phiA1,phiC1,phiA2,phiC2)
  Sig_theta <- stats::cov(phi_components)/n
  se_matrix_delta <- sqrt(max(0,drop(t(g_Delta) %*% Sig_theta %*% g_Delta)))
  if (abs(seD-se_matrix_delta) > 1e-8*max(1,seD))
    stop("Delta method vs IF discrepancy; abort analysis")
  loss_diff <- A1-A2
  se_loss_diff <- std_error((a1-a2)-loss_diff)
  z_loss <- if (se_loss_diff>0) loss_diff/se_loss_diff else NA_real_
  p_loss <- if (is.finite(z_loss)) 2*stats::pnorm(-abs(z_loss)) else NA_real_
  z <- if (seD > 0) delta/seD else NA_real_
  pval <- if (is.finite(z)) 2*stats::pnorm(-abs(z)) else NA_real_
  zcrit <- stats::qnorm(1-alpha/2)
  row <- data.frame(policy = policy_name, n = n, A1 = A1, A2 = A2,
                    C1 = C1, C2 = C2, R1 = R1, R2 = R2,
                    kappa1 = k1, kappa2 = k2, delta = delta,
                    se_kappa1 = se1, se_kappa2 = se2,
                    se_delta = seD, se_delta_multivariate_delta_method=se_matrix_delta,
                    delta_method_vs_eif_abs_difference=abs(seD-se_matrix_delta),
                    ci_delta_low = delta-zcrit*seD,
                    ci_delta_high = delta+zcrit*seD,
                    z = z, p_two_sided = pval,
                    mean_relative_fn_loss = mean(cval),
                    mean_loss_scale = mean(s),
                    risk_difference_A1_minus_A2 = loss_diff,
                    se_paired_loss_difference=se_loss_diff,
                    ci_loss_low=loss_diff-zcrit*se_loss_diff,
                    ci_loss_high=loss_diff+zcrit*se_loss_diff,
                    p_loss_two_sided=p_loss,
                    stringsAsFactors = FALSE)
  list(summary = row,
       influence = data.frame(phi_kappa1 = phiK1, phi_kappa2 = phiK2,
                              phi_delta = phiDelta),
       components = list(A1=A1,A2=A2,C1=C1,C2=C2,R1=R1,R2=R2),
       policy = cval)
}

# 09 | مشتقة كابا للاختبارين جزء من same-subject joint IF؛ لا نتعامل مع العينتين مستقلة.
# إجراء حساب منحنيات c الثابتة؛ CI هنا pointwise فقط وليس band على continuum.
loss_curves <- function(dat, nuisance, cgrid = seq(0.05,0.95,by=0.025)) {
  rows <- lapply(cgrid, function(cv) {
    a <- fixed_policy_inference(dat,nuisance,rep(cv,nrow(dat)),
                                policy_name="fixed_curve")$summary
    a$c_value <- cv
    a
  })
  do.call(rbind,rows)
}

# 10 | الهيمنة عبر كل دوال c(X) القابلة للقياس داخل مجال [L(X),U(X)].
# B(c)=C1*V2-C2*V1 + E[c(X)*h(X)]. إشارتها تساوي إشارة Delta إذا كانت R>0.
soft_min <- function(h, L, U, tau) {
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) stop("tau must be >0")
  a <- -L*h/tau; b <- -U*h/tau
  m <- pmax(a,b)
  -tau * (m + log(exp(a-m)+exp(b-m)))
}
soft_min_grad <- function(h, L, U, tau) {
  # weight on L: logistic((U-L)*h/tau); stable even with large arguments.
  pL <- stats::plogis((U-L)*h/tau)
  L*pL+U*(1-pL)
}

# 11 | V_j يمثل reference loss عند c(x)=0؛ نقدر V_j بمعادلة متعامدة.
prepare_envelope <- function(dat, nuisance,
                             s = loss_scale(dat$X),
                             L = box_L(dat$X), U = box_U(dat$X)) {
  check_loss(dat$X,s,rep(.5,nrow(dat)),L,U)
  D <- dat$D; T1 <- dat$T1; T2 <- dat$T2
  p <- nuisance$p; q1 <- nuisance$q1; q2 <- nuisance$q2
  eD <- D-p; e1 <- T1-q1; e2 <- T2-q2
  uC1 <- s*eD*e1; uC2 <- s*eD*e2
  C1 <- mean(uC1); C2 <- mean(uC2)
  v1 <- s*((1-p)*q1-q1*eD+(1-p)*e1)
  v2 <- s*((1-p)*q2-q2*eD+(1-p)*e2)
  V1 <- mean(v1); V2 <- mean(v2)
  phiC1 <- uC1-C1; phiC2 <- uC2-C2
  phiV1 <- v1-V1; phiV2 <- v2-V2
  h <- s*(C1*(p-q2)-C2*(p-q1))
  b0 <- C1*V2-C2*V1
  minval <- pmin(L*h,U*h)
  maxval <- pmax(L*h,U*h)
  d1 <- s*(p-q1); d2 <- s*(p-q2)
  minR1_plugin <- V1+mean(pmin(L*d1,U*d1))
  minR2_plugin <- V2+mean(pmin(L*d2,U*d2))
  if (min(minR1_plugin,minR2_plugin) <= CFG$ref_denominator_tol)
    warning("Box lower bound for Rj is small/non-positive: dominance sign uncertifiable")
  list(D=D,T1=T1,T2=T2,p=p,q1=q1,q2=q2,s=s,L=L,U=U,
       eD=eD,e1=e1,e2=e2,C1=C1,C2=C2,V1=V1,V2=V2,
       phiC1=phiC1,phiC2=phiC2,phiV1=phiV1,phiV2=phiV2,
       h=h,b0=b0,lower_plugin=b0+mean(minval),
       upper_plugin=b0+mean(maxval),
       minR1_plugin=minR1_plugin,minR2_plugin=minR2_plugin,
       xi=(C1-C2)*eD+C2*e1-C1*e2)
}

# 12 | EIF للحد الأملس theta_tau (منهجية V7): قابل للتفاضل فقط عند tau>0 ثابت.
# phi_theta = (V2+alpha1)phiC1 -(V1+alpha2)phiC2
#             + C1phiV2 - C2phiV1 + g-Eg + w*s*xi.
# لا تعتمد التقديرات المباشرة للحد غير الأملس inferentially.
smooth_dominance_test <- function(env, tau, alpha = CFG$alpha) {
  n <- length(env$p)
  g <- soft_min(env$h,env$L,env$U,tau)
  w <- soft_min_grad(env$h,env$L,env$U,tau)
  alpha1 <- mean(w*env$s*(env$p-env$q2))
  alpha2 <- mean(w*env$s*(env$p-env$q1))
  # one-step orthogonal estimate; not the raw nonsmooth plug-in.
  theta_hat <- env$b0+mean(g+w*env$s*env$xi)
  phi_theta <- (env$V2+alpha1)*env$phiC1-(env$V1+alpha2)*env$phiC2+
    env$C1*env$phiV2-env$C2*env$phiV1+
    (g-mean(g))+w*env$s*env$xi
  # estimated EIF has small nonzero sample mean; center for variance consistency.
  phi_theta <- phi_theta-mean(phi_theta)
  se <- stats::sd(phi_theta)/sqrt(n)
  zcrit <- stats::qnorm(1-alpha)
  lcb <- theta_hat-zcrit*se
  plug_theta <- env$b0+mean(g)  # optional diagnostic only
  # هذه المقامات قيود ضرورية لتفسير إشارة B على أنها إشارة Delta.
  eligible <- is.finite(env$minR1_plugin) && is.finite(env$minR2_plugin) &&
              min(env$minR1_plugin, env$minR2_plugin) > CFG$ref_denominator_tol
  row <- data.frame(tau=tau, alpha_used=alpha, theta_one_step=theta_hat,
                    theta_plugin=plug_theta,
                    se=se, lower_one_sided=lcb,
                    reject_robust_null=if (eligible) as.integer(is.finite(lcb)&&lcb>0)
                                       else NA_integer_,
                    plugin_reference_positivity=eligible,
                    lower_B_plugin=env$lower_plugin,
                    upper_B_plugin=env$upper_plugin,
                    maximum_envelope_gap=tau*log(2),
                    lower_R1_plugin=env$minR1_plugin,
                    lower_R2_plugin=env$minR2_plugin,
                    stringsAsFactors=FALSE)
  list(summary=row,influence=phi_theta)
}

# 13 | عدة قيم tau المحددة مسبقاً: Bonferroni يمنع الاختيار اللاحق غير المضبوط.
smooth_multiple_taus <- function(env, taus=CFG$theta_taus,
                                 alpha=CFG$alpha) {
  if (length(taus)<1L || any(!is.finite(taus)) || any(taus<=0))
    stop("Supply finite positive preselected smoothing parameters")
  if (anyDuplicated(taus)) stop("taus must be distinct")
  res <- lapply(taus,function(tt) smooth_dominance_test(env,tt,alpha/length(taus)))
  list(table=do.call(rbind,lapply(res,`[[`,"summary")), detail=res)
}

# 14 | الاختبارات الذاتية المستقلة؛ هذه الفحوص ليست بديلاً عن تشغيل كامل ومحاكاة.
run_self_tests_v8 <- function() {
  set.seed(661)
  x <- c(18,22,35,52,68,82)
  p <- c(.05,.10,.15,.25,.37,.48)
  se <- c(.65,.68,.75,.80,.82,.85)
  sp <- c(.96,.94,.93,.90,.88,.86)
  q <- p*se+(1-p)*(1-sp)
  for (cv in c(0,.25,.5,.75,1)) {
    F <- p*(1-se); G <- (1-p)*(1-sp)
    A <- cv*F+(1-cv)*G
    R <- cv*p*(1-q)+(1-cv)*(1-p)*q
    C <- p*(1-p)*(se+sp-1)
    stopifnot(max(abs(R-A-C))<1e-12)
  }
  for (tau in c(0.0001,.01,.1)) {
    h <- c(-.3,-.02,0,.01,.3)
    L <- rep(.25,length(h)); U <- rep(.85,length(h))
    z <- soft_min(h,L,U,tau)
    exact <- pmin(L*h,U*h)
    stopifnot(all(z<=exact+1e-12),all(z>=exact-tau*log(2)-1e-12))
    eps <- 1e-6
    num <- (soft_min(h+eps,L,U,tau)-soft_min(h-eps,L,U,tau))/(2*eps)
    stopifnot(max(abs(num-soft_min_grad(h,L,U,tau)))<1e-5)
  }
  # تحقق في عينة صغيرة من انتقال الملف إلى النماذج ومن حساب EIF دون انهيار.
  toy <- simulate_paired_diagnostic(n=180, scenario="alternative",seed=123)
  pred <- crossfit_probabilities(toy,K=3L,method="linear",df=3L,seed=888)
  est <- fixed_policy_inference(toy,pred,rep(.5,nrow(toy)),policy_name="selftest")
  stopifnot(is.finite(est$summary$delta), is.finite(est$summary$se_delta))
  env <- prepare_envelope(toy,pred)
  out <- smooth_dominance_test(env,tau=.001)
  stopifnot(is.finite(out$summary$theta_one_step),is.finite(out$summary$se))
  TRUE
}

# Additional V9 test: verify bias-reduced logistic returns finite coefficients
# even under a deliberately (quasi-)separated 0/1 response.
run_self_tests_v9 <- function() {
  stopifnot(isTRUE(run_self_tests_v8()))
  toy <- data.frame(X=seq(-2,2,length.out=40L))
  toy$D <- as.numeric(toy$X>0)
  m <- fit_safe_logistic(D ~ X,toy,"V9 separation challenge")
  pr <- as.numeric(stats::predict(m,newdata=toy,type="response"))
  stopifnot(all(is.finite(stats::coef(m))),all(is.finite(pr)),
            all(pr>0),all(pr<1))
  TRUE
}

# 15 | نموذج محلي مرن لفحص الحساسية والنوعية؛ لا يتطلب أربع خلايا مشتركة موجبة.
# ملاءمة هوامش Tj|D,X منفصلة لا تعني افتراض استقلال الاختبارين:
# الاستدلال يستبقي اقتران الأشخاص داخل paired subject-level bootstrap.
# 15 | Local ridge-penalized logistic model with a FROZEN SPLINE BASIS.
# This is a prespecified secondary conditional estimator, not the primary EIF estimator.
# Same penalty, basis and algorithm in the original sample AND every bootstrap sample.
# The ridge penalty is on standardized nonintercept regression coefficients.
# Only computational stability is targeted here, NOT empirical 95% band coverage.
# Fixing lambda avoids outcome-driven tuning inside a fragile rare-event application.
make_local_recipe <- function(dat, df_prev=CFG$local_df_prev,
                              df_pair=CFG$local_df_pair) {
  block <- function(df, with_d) {
    B <- splines::ns(dat$X, df=df, intercept=FALSE)
    knots <- attr(B,"knots"); boundary <- attr(B,"Boundary.knots")
    m <- if (with_d) cbind(D=dat$D,B, sweep(B,1,dat$D,"*")) else B
    m <- as.matrix(m)
    centers <- colMeans(m)
    scales <- apply(m,2,stats::sd)
    scales[!is.finite(scales) | scales<1e-6] <- 1
    list(knots=knots,boundary=boundary,center=centers,scale=scales,
         with_d=with_d,ncols=ncol(m))
  }
  list(prev=block(df_prev,FALSE),pair=block(df_pair,TRUE),
       x_range=range(dat$X), df_prev=df_prev, df_pair=df_pair)
}
local_design <- function(x, d=NULL, block) {
  basis <- as.matrix(splines::ns(x,knots=block$knots,
                         Boundary.knots=block$boundary,intercept=FALSE))
  if (block$with_d) {
    if (is.null(d) || length(d)!=length(x)) stop("local_design: D length mismatch")
    mm <- cbind(D=d,basis,sweep(basis,1,d,"*"))
  } else mm <- basis
  mm <- as.matrix(mm)
  if (ncol(mm)!=block$ncols) stop("Local spline matrix has unexpected dimension")
  zz <- sweep(sweep(mm,2,block$center,"-"),2,block$scale,"/")
  cbind(Intercept=1,zz)
}
# Strictly convex penalized binomial objective for lambda>0 (when intercept identifiable).
# Stable log(1+exp(eta)) avoids overflow; nlminb is bundled with base R.
# Robust convex ridge-logistic optimization (V10.1).
# Fixed objective and penalty as V10: only numerical strategy changes.
# A small gradient of the ORIGINAL penalized objective is the acceptance rule;
# optimizer status codes alone are not proof of failure or convergence.
# Damped Newton with a line search is followed, if needed, by nlminb/BFGS and
# another Newton refinement. All unresolved numerical failures remain reported.
fit_ridge_binomial <- function(xmat,y,lambda=CFG$local_ridge_lambda,label="model") {
  xmat <- as.matrix(xmat)
  storage.mode(xmat) <- "double"
  y <- as.numeric(y)
  if (length(y)!=nrow(xmat) || anyNA(y) || any(!(y %in% c(0,1))))
    stop(label,": binary outcome/design mismatch")
  if (any(!is.finite(xmat))) stop(label,": nonfinite model matrix")
  if (!is.finite(lambda) || lambda<=0) stop("lambda must be strictly positive")
  if (ncol(xmat)<1L || any(abs(xmat[,1]-1)>1e-12))
    stop(label,": first column must be a constant intercept")
  if (sum(y)==0L || sum(y)==length(y)) {
    stop(label,": only one outcome class in resample; unpenalized intercept not identified")
  }
  pcols <- ncol(xmat)
  pen <- c(0, rep(lambda,pcols-1L))
  gradient_tol <- CFG$ridge_gradient_tol
  if (!is.finite(gradient_tol) || gradient_tol<=0)
    stop("ridge_gradient_tol must be positive")
  objective <- function(beta) {
    eta <- drop(xmat %*% beta)
    sum(pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta) +
      0.5*sum(pen*beta^2)
  }
  gradient <- function(beta) {
    eta <- drop(xmat %*% beta)
    drop(crossprod(xmat, stats::plogis(eta)-y))+pen*beta
  }
  grad_norm <- function(beta) {
    if (is.null(beta) || length(beta)!=pcols || any(!is.finite(beta)))
      return(Inf)
    g <- tryCatch(gradient(beta),error=function(e) rep(Inf,pcols))
    if (any(!is.finite(g))) Inf else max(abs(g))
  }
  # Ridge penalty applies to all slope columns but not to the intercept.
  # With both 0 and 1 outcomes, the penalized optimum is finite and unique.
  newton <- function(initial,maxit=150L) {
    beta <- as.numeric(initial)
    iterations <- 0L
    status <- "iteration_limit"
    for (iter in seq_len(maxit)) {
      gn <- grad_norm(beta)
      if (gn <= gradient_tol) {
        status <- "stationary"; break
      }
      eta <- drop(xmat %*% beta)
      phat <- stats::plogis(eta)
      w <- phat*(1-phat)
      hess <- crossprod(xmat, sweep(xmat,1,w,"*")) + diag(pen,nrow=pcols)
      g <- gradient(beta)
      # Allow tiny solver damping only if needed to invert a near-singular Hessian.
      direction <- NULL
      for (jitter in c(0,1e-10,1e-8,1e-6)) {
        H <- if (jitter==0) hess else hess+diag(jitter,pcols)
        d <- tryCatch(drop(solve(H,g)),error=function(e) NULL)
        if (!is.null(d) && all(is.finite(d)) && sum(g*d)>0) {
          direction <- d; break
        }
      }
      if (is.null(direction)) {
        status <- "invalid_direction"; break
      }
      # Restrict the maximum change in any linear predictor per Newton step.
      max_eta_step <- max(abs(drop(xmat %*% direction)))
      if (!is.finite(max_eta_step)) {
        status <- "invalid_eta_step"; break
      }
      if (max_eta_step>20) direction <- direction*(20/max_eta_step)
      gd <- sum(g*direction)
      old <- objective(beta)
      accepted <- FALSE
      t <- 1
      for (ls in seq_len(45L)) {
        candidate <- beta-t*direction
        val <- objective(candidate)
        if (is.finite(val) && val<=old-1e-4*t*gd) {
          beta <- candidate; accepted <- TRUE; break
        }
        t <- t/2
      }
      if (!accepted) {
        status <- "line_search_stalled"; break
      }
      iterations <- iterations+1L
    }
    gn <- grad_norm(beta)
    list(par=beta,gradient_max=gn,objective=objective(beta),
         converged=is.finite(gn) && gn<=gradient_tol,
         iterations=iterations,status=if (gn<=gradient_tol) "stationary" else status)
  }
  initial <- rep(0,pcols)
  initial[1] <- stats::qlogis((sum(y)+.5)/(length(y)+1))
  best <- newton(initial)
  best_method <- "damped_Newton"
  optimizer_status <- best$status
  total_iterations <- best$iterations
  if (!best$converged) {
    # Use the best available solution as a warm start, not restart at zero.
    alt1 <- tryCatch(stats::nlminb(best$par,objective,gradient,
              control=list(iter.max=1200L,eval.max=2400L,rel.tol=1e-12)),
              error=function(e)e)
    if (!inherits(alt1,"error") && all(is.finite(alt1$par))) {
      a <- newton(alt1$par,maxit=150L)
      total_iterations <- total_iterations+a$iterations
      if (is.finite(a$gradient_max) && a$gradient_max<best$gradient_max) {
        best <- a; best_method <- "nlminb_then_Newton"
        optimizer_status <- paste0("nlminb_code_",alt1$convergence,"_",a$status)
      }
    }
  }
  if (!best$converged) {
    alt2 <- tryCatch(stats::optim(best$par,objective,gradient,method="BFGS",
              control=list(maxit=2500L,reltol=1e-13)),error=function(e)e)
    if (!inherits(alt2,"error") && all(is.finite(alt2$par))) {
      a <- newton(alt2$par,maxit=150L)
      total_iterations <- total_iterations+a$iterations
      if (is.finite(a$gradient_max) && a$gradient_max<best$gradient_max) {
        best <- a; best_method <- "BFGS_then_Newton"
        optimizer_status <- paste0("BFGS_code_",alt2$convergence,"_",a$status)
      }
    }
  }
  if (!best$converged || any(!is.finite(best$par))) {
    stop(label,": ridge solver failed after Newton/nlminb/BFGS; max_gradient=",
         format(best$gradient_max,digits=6),
         ", tolerance=",format(gradient_tol,digits=6),
         ", status=",optimizer_status)
  }
  list(beta=as.numeric(best$par),lambda=lambda,method=best_method,
       gradient_max=best$gradient_max, gradient_tol=gradient_tol,
       iterations=total_iterations,optimizer_status=optimizer_status,
       n_event=sum(y),objective=best$objective)
}
ridge_local_predict <- function(model,x,d=NULL,block) {
  mm <- local_design(x,d,block)
  clip01(stats::plogis(drop(mm%*%model$beta)))
}
fit_marginal_spline <- function(dat,df_prev=CFG$local_df_prev,
                               df_pair=CFG$local_df_pair,
                               reference_recipe=NULL) {
  recipe <- if (is.null(reference_recipe))
    make_local_recipe(dat,df_prev,df_pair) else reference_recipe
  xp <- local_design(dat$X,block=recipe$prev)
  xt <- local_design(dat$X,dat$D,block=recipe$pair)
  p <- fit_ridge_binomial(xp,dat$D,label="D|X")
  t1 <- fit_ridge_binomial(xt,dat$T1,label="T1|D,X")
  t2 <- fit_ridge_binomial(xt,dat$T2,label="T2|D,X")
  fits <- list(p,t1,t2)
  list(kind="ridge_marginal",prev=p,test1=t1,test2=t2,
       recipe=recipe,x_min=recipe$x_range[1],x_max=recipe$x_range[2],
       audit=data.frame(model=c("D|X","T1|D,X","T2|D,X"),
         estimator="ridge_penalized_logistic",lambda=CFG$local_ridge_lambda,
         n_events=vapply(fits,function(z)z$n_event,numeric(1)),
         optimizer=vapply(fits,function(z)z$method,character(1)),
         optimizer_status=vapply(fits,function(z)z$optimizer_status,character(1)),
         solver_iterations=vapply(fits,function(z)z$iterations,integer(1)),
         objective=vapply(fits,function(z)z$objective,numeric(1)),
         max_gradient=vapply(fits,function(z)z$gradient_max,numeric(1)),
         gradient_tolerance=CFG$ridge_gradient_tol))
}

# 15b | benchmark V5 joint multinomial P(T1,T2|D,X): secondary only.
# هذا التحليل قد يفشل لأن الخلية T1=1,T2=0 نادرة للغاية في بيانات UCI.
fit_joint_spline <- function(dat, df_prev=CFG$local_df_prev,
                             df_pair=CFG$local_df_pair) {
  if (!all(c("D","T1","T2","X") %in% names(dat))) stop("Missing columns")
  train <- dat
  train$pair <- factor(paste0(train$T1,train$T2),
                       levels=c("00","01","10","11"))
  if (any(table(train$pair)==0L)) stop("All 4 paired outcome types are required")
  fprev <- stats::as.formula(paste0("D ~ splines::ns(X, df = ",df_prev,")"))
  fpair <- stats::as.formula(paste0("pair ~ D * splines::ns(X, df = ",df_pair,")"))
  mprev <- fit_safe_logistic(fprev,train,"JOINT benchmark prevalence")
  # Penalized multinomial is secondary; it is not a drop-in inferential replacement.
  mpair <- nnet::multinom(fpair, data=train, trace=FALSE,
                          decay=0.0001, maxit=2500, MaxNWts=10000)
  if (!isTRUE(mprev$converged) ||
      any(!is.finite(stats::coef(mprev))) ||
      (!is.null(mpair$convergence) && mpair$convergence!=0))
    stop("Spline/multinomial model failed: report failure")
  list(kind="joint",prev=mprev,pair=mpair,levels=levels(train$pair),
       x_min=min(dat$X),x_max=max(dat$X))
}

joint_prob_matrix <- function(fit, x, d) {
  nd <- data.frame(X=x,D=rep(d,length(x)))
  ans <- stats::predict(fit$pair,newdata=nd,type="probs")
  if (is.vector(ans)) ans <- matrix(ans,nrow=1,dimnames=list(NULL,names(ans)))
  ans <- as.matrix(ans)
  if (!all(fit$levels %in% colnames(ans))) stop("Paired class predictions missing")
  ans[,fit$levels,drop=FALSE]
}

predict_local_components <- function(fit,x,cval=.5) {
  if (any(!is.finite(x))) stop("Invalid x grid")
  if (any(x<fit$x_min | x>fit$x_max))
    stop("Do not extrapolate outside X observed in fitted sample")
  if (length(cval)==1) cval <- rep(cval,length(x))
  if (length(cval)!=length(x) || any(cval<0|cval>1)) stop("Bad c grid")
  if (identical(fit$kind,"ridge_marginal")) {
    p <- ridge_local_predict(fit$prev,x,block=fit$recipe$prev)
    Se1 <- ridge_local_predict(fit$test1,x,rep(1,length(x)),fit$recipe$pair)
    Se2 <- ridge_local_predict(fit$test2,x,rep(1,length(x)),fit$recipe$pair)
    Sp1 <- 1-ridge_local_predict(fit$test1,x,rep(0,length(x)),fit$recipe$pair)
    Sp2 <- 1-ridge_local_predict(fit$test2,x,rep(0,length(x)),fit$recipe$pair)
  } else if (identical(fit$kind,"marginal")) {
    p <- clip01(stats::predict(fit$prev,newdata=data.frame(X=x),type="response"))
    nd0 <- data.frame(X=x,D=rep(0,length(x)))
    nd1 <- data.frame(X=x,D=rep(1,length(x)))
    Se1 <- clip01(stats::predict(fit$test1,newdata=nd1,type="response"))
    Se2 <- clip01(stats::predict(fit$test2,newdata=nd1,type="response"))
    Sp1 <- 1-clip01(stats::predict(fit$test1,newdata=nd0,type="response"))
    Sp2 <- 1-clip01(stats::predict(fit$test2,newdata=nd0,type="response"))
  } else if (identical(fit$kind,"joint")) {
    p <- clip01(stats::predict(fit$prev,newdata=data.frame(X=x),type="response"))
    d0 <- joint_prob_matrix(fit,x,0)
    d1 <- joint_prob_matrix(fit,x,1)
    Se1 <- d1[,"10"]+d1[,"11"]
    Se2 <- d1[,"01"]+d1[,"11"]
    Sp1 <- d0[,"00"]+d0[,"01"]
    Sp2 <- d0[,"00"]+d0[,"10"]
  } else stop("Unknown local fit kind")
  q1 <- p*Se1+(1-p)*(1-Sp1)
  q2 <- p*Se2+(1-p)*(1-Sp2)
  C1 <- p*(1-p)*(Se1+Sp1-1)
  C2 <- p*(1-p)*(Se2+Sp2-1)
  F1 <- p*(1-Se1); G1 <- (1-p)*(1-Sp1)
  F2 <- p*(1-Se2); G2 <- (1-p)*(1-Sp2)
  A1 <- cval*F1+(1-cval)*G1
  A2 <- cval*F2+(1-cval)*G2
  R1 <- cval*p*(1-q1)+(1-cval)*(1-p)*q1
  R2 <- cval*p*(1-q2)+(1-cval)*(1-p)*q2
  k1 <- ifelse(R1>CFG$ref_denominator_tol,C1/R1,NA_real_)
  k2 <- ifelse(R2>CFG$ref_denominator_tol,C2/R2,NA_real_)
  U1 <- p*(1-q1); V1 <- (1-p)*q1
  U2 <- p*(1-q2); V2 <- (1-p)*q2
  bb0 <- C1*V2-C2*V1
  bb1 <- C1*(U2-V2)-C2*(U1-V1)
  cr <- ifelse(abs(bb1)>1e-11, -bb0/bb1, NA_real_)
  cr <- ifelse(cr>=0 & cr<=1, cr, NA_real_)
  data.frame(X=x,c=cval,p=p,Se1=Se1,Se2=Se2,Sp1=Sp1,Sp2=Sp2,
             q1=q1,q2=q2,A1=A1,A2=A2,R1=R1,R2=R2,C1=C1,C2=C2,
             kappa1=k1,kappa2=k2,delta=k1-k2,crossover_c=cr)
}

# 16 | سطح delta(x,c)؛ نفس الموديل المشترك، بدون افتراض one-crossover عبر x.
local_surface <- function(fit,x_grid,c_grid) {
  do.call(rbind,lapply(c_grid,function(cv)
    predict_local_components(fit,x_grid,rep(cv,length(x_grid)))))
}

# 17 | benchmark integrated through standardization (not cross-fitted).
joint_integrated_policy <- function(fit,x,cval,s=loss_scale(x)) {
  z <- predict_local_components(fit,x,cval)
  A1 <- mean(s*z$A1); A2 <- mean(s*z$A2)
  C1 <- mean(s*z$C1); C2 <- mean(s*z$C2)
  data.frame(benchmark=if (identical(fit$kind,"ridge_marginal"))
                  "ridge_marginal_logistic_spline" else if
                  (identical(fit$kind,"marginal")) "marginal_logistic_spline"
                  else "joint_multinomial_spline",
             n=length(x),
             kappa1=C1/(A1+C1),kappa2=C2/(A2+C2),
             delta=C1/(A1+C1)-C2/(A2+C2),
             stringsAsFactors=FALSE)
}

# 18 | paired subject-bootstrap للسطح على finite prespecified grid.
# يعمل مع النماذج الثنائية الصحيحة فقط؛ فشل النماذج يُبلّغ ولا يُخفى.
bootstrap_joint_local <- function(dat,fit,x_grid,c_grid,
                                  B=CFG$bootstrap_B,
                                  seed=CFG$seed+25L) {
  if (B<50L) stop("At least B=50 (pilot); B>=999 recommended for publication")
  xref <- local_surface(fit,x_grid,c_grid)
  n <- nrow(dat)
  set.seed(seed)
  draws <- matrix(NA_real_,nrow=B,ncol=nrow(xref))
  attempts <- vector("list",B)
  for (b in seq_len(B)) {
    idx <- sample.int(n,n,replace=TRUE)
    one <- tryCatch({
      fb <- if (identical(fit$kind,"ridge_marginal"))
        fit_marginal_spline(dat[idx,,drop=FALSE],reference_recipe=fit$recipe)
      else if (identical(fit$kind,"marginal"))
        fit_marginal_spline(dat[idx,,drop=FALSE])
      else fit_joint_spline(dat[idx,,drop=FALSE])
      estimates <- local_surface(fb,x_grid,c_grid)$delta
      if (length(estimates)!=nrow(xref) || any(!is.finite(estimates)))
        stop("local grid produced missing/non-finite Delta")
      n_warn <- if (!is.null(fb$audit) && "warning_count" %in% names(fb$audit))
        sum(fb$audit$warning_count) else 0L
      warn_messages <- if (!is.null(fb$audit) &&
                           "warning_messages" %in% names(fb$audit))
        paste(fb$audit$warning_messages[nzchar(fb$audit$warning_messages)],
              collapse=" | ") else ""
      list(delta=estimates,warnings=n_warn,warning_messages=warn_messages,
           optimizers=paste(fb$audit$optimizer,collapse=";"),
           max_grad=max(fb$audit$max_gradient),
           iterations=sum(fb$audit$solver_iterations))
    }, error=function(e)e)
    if (inherits(one,"error")) {
      attempts[[b]] <- data.frame(replication=b,success=FALSE,
                                  warning_count=NA_integer_,
                                  reason=conditionMessage(one),
                                  warning_messages="",optimizers="",
                                  max_gradient=NA_real_,
                                  solver_iterations=NA_integer_,
                                  stringsAsFactors=FALSE)
    } else {
      draws[b,] <- one$delta
      attempts[[b]] <- data.frame(replication=b,success=TRUE,
                                  warning_count=one$warnings,
                                  reason="",warning_messages=one$warning_messages,
                                  optimizers=one$optimizers,
                                  max_gradient=one$max_grad,
                                  solver_iterations=one$iterations,
                                  stringsAsFactors=FALSE)
    }
    if (b%%50L==0L) logline("Paired bootstrap: ",b," / ",B)
  }
  attempt_table <- do.call(rbind,attempts)
  good <- which(attempt_table$success)
  n_good <- length(good)
  frac <- n_good/B
  certified <- frac>=CFG$bootstrap_min_success && n_good>=50L
  logline("Bootstrap successful: ",n_good,"/",B,"; accepted=",certified)
  failed_reasons <- attempt_table$reason[!attempt_table$success]
  freq <- if (length(failed_reasons)) {
    t <- as.data.frame(table(failed_reasons),stringsAsFactors=FALSE)
    names(t) <- c("reason","count"); t
  } else data.frame(reason=character(),count=integer())
  status <- data.frame(requested=B,success=n_good,failures=B-n_good,
                       success_rate=frac,min_required=CFG$bootstrap_min_success,
                       accepted_for_finite_grid_CI=certified,
                       method=paste0(fit$kind," / subject-level resampling"),
                       warning_bearing_successes=sum(attempt_table$success &
                           !is.na(attempt_table$warning_count) &
                           attempt_table$warning_count>0),
                       highest_successful_gradient=if (n_good>0L)
                         max(attempt_table$max_gradient[good]) else NA_real_,
                       retry_optimizer_successes=sum(grepl("then_Newton",
                         attempt_table$optimizers[good],fixed=TRUE)),
                       max_t_critical=NA_real_, stringsAsFactors=FALSE)
  # CRITICAL: preserve the attempt-by-attempt audit even if bands are rejected.
  if (!certified) {
    return(list(results=NULL,status=status,attempts=attempt_table,
                failure_messages=freq))
  }
  if (frac<0.99) warning("Some bootstrap fits failed: bands conditional on successful fits; verify by Monte Carlo")
  dd <- draws[good,,drop=FALSE]
  se <- apply(dd,2,stats::sd)
  lo <- apply(dd,2,stats::quantile,probs=CFG$alpha/2,na.rm=TRUE,names=FALSE)
  hi <- apply(dd,2,stats::quantile,probs=1-CFG$alpha/2,na.rm=TRUE,names=FALSE)
  z <- sweep(dd,2,xref$delta,FUN="-")
  positive <- which(se>1e-10)
  if (length(positive)>0L) {
    zz <- sweep(z[,positive,drop=FALSE],2,se[positive],FUN="/")
    mx <- apply(abs(zz),1,max)
    mult <- as.numeric(stats::quantile(mx,1-CFG$alpha,na.rm=TRUE,names=FALSE))
  } else mult <- NA_real_
  status$max_t_critical <- mult
  xref$bootstrap_se <- se
  xref$ci_point_low <- lo; xref$ci_point_high <- hi
  xref$simultaneous_low <- if (is.finite(mult)) xref$delta-mult*se else NA_real_
  xref$simultaneous_high <- if (is.finite(mult)) xref$delta+mult*se else NA_real_
  list(results=xref,status=status,attempts=attempt_table,
       failure_messages=freq)
}

# 19 | حفظ رسوم base R بدون اعتماد ggplot2؛ كل ملف له عنوان ودلالة واضحة.
plot_fixed_curve <- function(curve, file) {
  grDevices::png(file,width=1300,height=850,res=140)
  on.exit(grDevices::dev.off(),add=TRUE)
  plot(curve$c_value,curve$kappa1,type="l",lwd=3,col="#22577A",
       ylim=range(c(curve$kappa1,curve$kappa2),finite=TRUE),
       xlab="Prespecified constant false-negative loss weight c",
       ylab="Integrated conditional-chance kappa",
       main="Cervical cancer: loss-weighted kappa curves")
  lines(curve$c_value,curve$kappa2,lwd=3,col="#B23A48")
  legend("bottomright",legend=c("Hinselmann","Schiller"),
         col=c("#22577A","#B23A48"),lwd=3,bty="n")
}
plot_delta_curve <- function(curve,file) {
  grDevices::png(file,width=1300,height=850,res=140)
  on.exit(grDevices::dev.off(),add=TRUE)
  plot(curve$c_value,curve$delta,type="l",lwd=3,col="#22577A",
       xlab="Constant FN-loss weight c",ylab="Paired difference: kappa1 - kappa2",
       main="Cross-fitted integrated difference (pointwise estimates)")
  abline(h=0,lty=2,col="grey40")
  lines(curve$c_value,curve$ci_delta_low,lty=3,col="grey40")
  lines(curve$c_value,curve$ci_delta_high,lty=3,col="grey40")
  legend("topright",legend=c("Delta","Pointwise 95% Wald CI"),
         col=c("#22577A","grey40"),lty=c(1,3),lwd=c(3,1),bty="n")
}
plot_local_surface <- function(surface,file) {
  grDevices::png(file,width=1300,height=850,res=140)
  on.exit(grDevices::dev.off(),add=TRUE)
  xx <- sort(unique(surface$X)); cc <- sort(unique(surface$c))
  zz <- matrix(NA_real_,nrow=length(xx),ncol=length(cc))
  for (j in seq_along(cc)) {
    rr <- surface[abs(surface$c-cc[j])<1e-12,,drop=FALSE]
    rr <- rr[order(rr$X),,drop=FALSE]
    zz[,j] <- rr$delta
  }
  graphics::image(xx,cc,zz,xlab="Age (years)",ylab="FN relative loss c",
                  main="Local paired weighted-kappa difference: finite grid",
                  col=grDevices::hcl.colors(45,"Blue-Red 3"))
  graphics::contour(xx,cc,zz,levels=0,add=TRUE,lwd=2,drawlabels=FALSE)
}

# 20 | توليد DGM مرتبط داخل الشخص للمحاكاة: conditional marginals with shared U.
# under null: identical test marginals but NON-identical paired realizations.
simulate_paired_diagnostic <- function(n,scenario=c("null","alternative"),
                                       seed=NULL, dependence=.55,
                                       truth_model=c("constant","varying")) {
  scenario <- match.arg(scenario)
  truth_model <- match.arg(truth_model)
  if (!is.null(seed)) set.seed(seed)
  X <- pmin(84,pmax(13,round(stats::rnorm(n,mean=30,sd=12))))
  z <- (X-30)/12
  if (truth_model == "constant") {
    # p, Se, Sp ثابتة: true q_j constant؛ linear and spline nuisance models contain truth.
    p <- rep(stats::plogis(-1.7),n)
    Se1 <- rep(stats::plogis(.8),n)
    Sp1 <- rep(stats::plogis(2.2),n)
    if (scenario=="null") { Se2<-Se1; Sp2<-Sp1 } else {
      Se2 <- rep(stats::plogis(.55),n)
      Sp2 <- rep(stats::plogis(1.65),n)
    }
  } else {
    # هذا سيناريو حساسية قد تكون فيه نماذج q_j اللوجستية خاطئة التحديد.
    p <- stats::plogis(-1.7+.35*z)
    Se1 <- stats::plogis(.8+.22*z)
    Sp1 <- stats::plogis(2.2-.15*z)
    if (scenario=="null") { Se2 <- Se1; Sp2 <- Sp1 } else {
      Se2 <- stats::plogis(.55+.07*z)
      Sp2 <- stats::plogis(1.65+.03*z)
    }
  }
  D <- stats::rbinom(n,1,p)
  r1 <- ifelse(D==1,Se1,1-Sp1)
  r2 <- ifelse(D==1,Se2,1-Sp2)
  u1 <- stats::runif(n)
  u2 <- ifelse(stats::runif(n)<dependence,u1,stats::runif(n))
  data.frame(X=X,D=D,T1=as.integer(u1<r1),T2=as.integer(u2<r2))
}

# 21 | القيمة الحقيقية في DGM عبر تكامل قطعي عددي deterministic للـ X.
# X=round(clamp(N(30,12^2),13,84)); qnorm midpoint grid approximates distribution.
dgm_population_truth <- function(scenario=c("null","alternative"),
                                 policy="equal",m=200001L, tau=.0005,
                                 truth_model=c("constant","varying")) {
  scenario <- match.arg(scenario)
  truth_model <- match.arg(truth_model)
  u <- (seq_len(m)-0.5)/m
  X <- pmin(84,pmax(13,round(stats::qnorm(u,mean=30,sd=12))))
  z <- (X-30)/12
  if (truth_model=="constant") {
    p <- rep(stats::plogis(-1.7),m)
    Se1 <- rep(stats::plogis(.8),m); Sp1 <- rep(stats::plogis(2.2),m)
    if (scenario=="null") { Se2<-Se1; Sp2<-Sp1 } else {
      Se2<-rep(stats::plogis(.55),m); Sp2<-rep(stats::plogis(1.65),m)
    }
  } else {
    p <- stats::plogis(-1.7+.35*z)
    Se1 <- stats::plogis(.8+.22*z); Sp1 <- stats::plogis(2.2-.15*z)
    if (scenario=="null") { Se2<-Se1; Sp2<-Sp1 } else {
      Se2<-stats::plogis(.55+.07*z); Sp2<-stats::plogis(1.65+.03*z)
    }
  }
  q1 <- p*Se1+(1-p)*(1-Sp1)
  q2 <- p*Se2+(1-p)*(1-Sp2)
  s <- loss_scale(X); cval <- relative_fn_loss(X,policy)
  C1 <- mean(s*p*(1-p)*(Se1+Sp1-1))
  C2 <- mean(s*p*(1-p)*(Se2+Sp2-1))
  A1 <- mean(s*(cval*p*(1-Se1)+(1-cval)*(1-p)*(1-Sp1)))
  A2 <- mean(s*(cval*p*(1-Se2)+(1-cval)*(1-p)*(1-Sp2)))
  delta <- C1/(A1+C1)-C2/(A2+C2)
  V1 <- mean(s*(1-p)*q1); V2 <- mean(s*(1-p)*q2)
  h <- s*(C1*(p-q2)-C2*(p-q1))
  L <- box_L(X); U <- box_U(X)
  b0 <- C1*V2-C2*V1
  data.frame(scenario=scenario, truth_model=truth_model,policy=policy,delta=delta,
             lower_B=b0+mean(pmin(L*h,U*h)),
             theta_tau=b0+mean(soft_min(h,L,U,tau)),tau=tau)
}

# 22 | محاكاة مستقلة للنماذج والتحليل الأساسي + conservative dominance.
# الفشل جزء من التقرير؛ لا نعطي SE/coverage مصطنعاً للنسخ الفاشلة.
run_monte_carlo <- function(n, replications=CFG$simulation_reps,
                            scenario=c("null","alternative"),
                            policy="equal",tau=.0005, seed=1699,
                            truth_model=c("constant","varying")) {
  scenario <- match.arg(scenario)
  truth_model <- match.arg(truth_model)
  truth <- dgm_population_truth(scenario,policy,tau=tau,truth_model=truth_model)
  rows <- vector("list",replications)
  failures <- character(0)
  for (r in seq_len(replications)) {
    set.seed(seed+100000L*r+n)
    dat <- simulate_paired_diagnostic(n,scenario,truth_model=truth_model)
    result <- tryCatch({
      nuisance <- crossfit_probabilities(dat,method=CFG$nuisance_model,
                                         df=CFG$nuisance_df,K=CFG$folds,
                                         seed=seed+r)
      cval <- relative_fn_loss(dat$X,policy)
      f <- fixed_policy_inference(dat,nuisance,cval,policy_name=policy)
      e <- prepare_envelope(dat,nuisance)
      d <- smooth_dominance_test(e,tau=tau)
      data.frame(replication=r,delta_hat=f$summary$delta,
                 se_delta=f$summary$se_delta,
                 ci_low=f$summary$ci_delta_low,
                 ci_high=f$summary$ci_delta_high,
                 reject_delta=as.integer(is.finite(f$summary$p_two_sided) &&
                                          f$summary$p_two_sided<CFG$alpha),
                 theta_hat=d$summary$theta_one_step,
                 lower_dominance=d$summary$lower_one_sided,
                 reject_dominance=d$summary$reject_robust_null)
    },error=function(err) err)
    if (inherits(result,"error")) failures <- c(failures,conditionMessage(result))
    else rows[[r]] <- result
    if (r %% 50L == 0L) logline("MC ",scenario," n=",n,": ",r,"/",replications)
  }
  good <- Filter(Negate(is.null),rows)
  raw <- if (length(good)) do.call(rbind,good) else data.frame()
  if (!nrow(raw)) stop("Every Monte Carlo fit failed")
  ss <- data.frame(
    n=n,scenario=scenario,truth_model=truth_model,policy=policy,
    nuisance_model=CFG$nuisance_model,
    logistic_method=CFG$logistic_strategy,replications=replications,
    successes=nrow(raw),failure_rate=1-nrow(raw)/replications,
    true_delta=truth$delta,
    delta_bias=mean(raw$delta_hat-truth$delta),
    delta_RMSE=sqrt(mean((raw$delta_hat-truth$delta)^2)),
    delta_empirical_SD=stats::sd(raw$delta_hat),
    delta_mean_SE=mean(raw$se_delta),
    delta_coverage=mean(raw$ci_low<=truth$delta & raw$ci_high>=truth$delta),
    delta_mean_CI_width=mean(raw$ci_high-raw$ci_low),
    delta_reject_fraction=mean(raw$reject_delta),
    lower_B_truth=truth$lower_B,
    theta_tau_truth=truth$theta_tau,
    dominance_reject_fraction=mean(raw$reject_dominance))
  ss$coverage_MCSE <- sqrt(ss$delta_coverage*(1-ss$delta_coverage)/nrow(raw))
  list(summary=ss,raw=raw,failed_messages=as.data.frame(table(failures)))
}

# 23 | بحث حجم العينة بالمحاكاة؛ لا توجد صيغة مغلقة تدّعي الدقة لهذه الإجراءات.
run_sample_size_search <- function(ns=CFG$simulation_ns,
                                   reps=CFG$simulation_reps,seed=500,
                                   target_power=0.80) {
  if (target_power<=0 || target_power>=1)
    stop("target_power must be strictly between 0 and 1")
  out <- list()
  for (nn in sort(unique(ns))) {
    obj <- run_monte_carlo(nn,reps,scenario="alternative",
                           policy="equal",seed=seed+nn,truth_model="constant")
    out[[as.character(nn)]] <- obj$summary
  }
  ans <- do.call(rbind,out)
  ans$target_power <- target_power
  ans$meets_target <- ans$delta_reject_fraction >= target_power
  eligible <- ans$n[ans$meets_target & is.finite(ans$delta_reject_fraction)]
  ans$smallest_eligible_n <- if (length(eligible)) min(eligible) else NA_integer_
  ans$power_MCSE <- sqrt(ans$delta_reject_fraction *
                          (1-ans$delta_reject_fraction)/ans$successes)
  # Candidate-grid exploratory planning only; power MCSE / model misspecification matter.
  ans
}

# 24 | تقرير README تلقائي يحدد معنى كل نتيجة وحدود الاستدلال.
write_readme <- function(dat,boot_status=NULL) {
  txt <- c(
    "JBS V10 | STABLE LOCAL BOOTSTRAP + JOURNAL-READY EXPORTS | NOT PEER REVIEWED",
    paste("Execution timestamp:",as.character(Sys.time())),
    paste("R version:",R.version.string),
    paste("Sample size:",nrow(dat),"Biopsy positives:",sum(dat$D)),
    "Source: UCI Cervical Cancer (Risk Factors), 858 patient-level rows.",
    "X=Age, D=Biopsy, T1=Hinselmann, T2=Schiller.",
    "Biopsy is a recorded reference, not proven error-free or clinically independent.",
    "All 858 rows preserved in RStudio_Input; this program assumes complete D/T/X.",
    "LOSS POLICIES ARE ILLUSTRATIVE, NOT OBTAINED FROM CLINICAL COST STUDIES.",
    "PRIMARY: cross-fitted orthogonal integrated kappa and paired Delta.",
    "Paired mean clinical loss difference and its influence-based SE are supplementary targets.",
    "LIMITATION: statistical validity assumes nuisance L2 consistency/product rates.",
    "Losses s(X),c(X), and [L(X),U(X)] must be fixed independently of outcomes.",
    "ROBUST DOMINANCE: preselected FIXED tau soft lower envelope and Bonferroni.",
    "The test is asymptotically POINTWISE conservative, not uniformly calibrated.",
    "With fixed tau>0, the asymptotic rejection probability under each fixed null is ZERO.",
    "The test cannot certify dominance when its smooth target remains negative.",
    "Plug-in lower_B is DESCRIPTIVE ONLY. It has no claimed Wald/ordinary bootstrap CI.",
    "LOCAL regression uses FIXED ridge-penalized logistic spline (not brglm2).",
    paste("Prespecified ridge penalty:",CFG$local_ridge_lambda),
    paste("Ridge absolute gradient tolerance (unchanged):",CFG$ridge_gradient_tol),
    "Local V10.1 solver: damped Newton + line search; warm-start nlminb/BFGS if needed.",
    "The V10.1 numerical solver does not change the V10 ridge objective/penalty; conditional-on-success bootstrap bands still need empirical coverage checks.",
    "Cross-fit nuisances ALSO use brglm2; this does not guarantee L2 rates in finite samples.",
    "Multivariate Delta SE is checked against paired EIF SE to numerical precision.",
    "LOCAL marginal logistic/spline analysis is SECONDARY and model-dependent.",
    "The paired bootstrap re-samples WHOLE subjects, retaining T1-T2 dependence.",
    "Joint multinomial benchmark is exploratory due to one T1+T2- joint record.",
    "Local simultaneous bootstrap confidence regions refer to a FINITE grid only.",
    "Near-zero positivity, sparse Biopsy-positive outcomes and separation can destabilize fits.",
    "Monte Carlo results are not final without ample replications and independent checking.",
    "MC constant DGM includes correctly specified regressions; varying DGM stresses misspecification.",
    "Sample-size grid is exploratory: smallest candidate meeting target power, not a closed-form guaranteed minimum.",
    "Self-tests are not Monte Carlo evidence; local ridge and main brglm2 fitting alter finite-sample estimators.",
    "Acceptance by Journal of Biopharmaceutical Statistics CANNOT be guaranteed.",
    paste("Cross-fitting K:",CFG$folds,"model:",CFG$nuisance_model,
          "df:",CFG$nuisance_df),
    paste("Smooth taus:",paste(CFG$theta_taus,collapse=", ")),
    paste("Primary clinical age-policy: invlogit(",CFG$age_policy_logit_at_center,
          " + ",CFG$age_policy_age_slope,"*(age-",CFG$age_policy_center,"))",sep=""),
    paste("Envelope bounds:",CFG$envelope_L,"to",CFG$envelope_U),
    paste("Local bootstrap enabled:",CFG$bootstrap_local,
          "B:",CFG$bootstrap_B),
    paste("R package versions:",paste(paste(required_packages,
       vapply(required_packages,function(p) as.character(utils::packageVersion(p)),""),sep="="),collapse="; "))
  )
  if (!is.null(boot_status)) txt <- c(txt,
    paste("Local bootstrap fit success rate:",boot_status$success_rate),
    paste("Local bootstrap bands accepted:",boot_status$accepted_for_finite_grid_CI),
    paste("Local bootstrap failures:",boot_status$failures))
  writeLines(txt,file.path(CFG$out_dir,"README_Results_and_Limitations.txt"),useBytes=TRUE)
  invisible(txt)
}

# 25 | نقطة التشغيل الموحدة: Source is enough if Excel lies beside script in working dir.
# V10 reproducibility check for the original/bootstrap identical ridge recipe.
run_self_tests_v10 <- function() {
  stopifnot(isTRUE(run_self_tests_v9()))
  toy <- simulate_paired_diagnostic(220,scenario="alternative",seed=738,
                                    truth_model="varying")
  fit <- fit_marginal_spline(toy)
  ages <- as.numeric(stats::quantile(toy$X,c(.1,.5,.9)))
  first <- predict_local_components(fit,ages,0.5)
  second <- fit_marginal_spline(toy[sample.int(nrow(toy),nrow(toy),TRUE),],
                                reference_recipe=fit$recipe)
  pred <- predict_local_components(second,ages,0.5)
  stopifnot(all(is.finite(first$delta)),all(is.finite(pred$delta)),
            all(first$R1>0),all(first$R2>0))
  # V10.1 numerical optimizer stress tests, without altering the estimand:
  # strong separation, rare events, and duplicate predictors.
  zz <- seq(-2,2,length.out=120L)
  Xm <- cbind(Intercept=1, X=zz, duplicate_X=zz, nonlinear=zz^2)
  yy <- as.numeric(zz>1.25)
  testfit <- fit_ridge_binomial(Xm,yy,label="stress: near separation")
  stopifnot(is.finite(testfit$objective),
            testfit$gradient_max<=CFG$ridge_gradient_tol,
            all(is.finite(testfit$beta)))
  yrare <- rep(0,120L)
  yrare[c(15L,71L,110L)] <- 1
  testfit2 <- fit_ridge_binomial(Xm,yrare,label="stress: rare 3/120")
  stopifnot(is.finite(testfit2$objective),
            testfit2$gradient_max<=CFG$ridge_gradient_tol,
            all(is.finite(testfit2$beta)))
  invisible(TRUE)
}

# 24b | Publication tables and 600-dpi figures (BASE R, no extra dependencies).
# Journal submission policy: captions, disclosures and real-data validation remain
# author responsibilities. No graphics or code can guarantee acceptance.
publication_outputs <- function(dat,qc,primary,curve,dominance,nuisance,
                                local=NULL,bs=NULL,benchmark=NULL) {
  tabdir <- file.path(CFG$out_dir,"JBS_Tables")
  figdir <- file.path(CFG$out_dir,"JBS_Figures")
  dir.create(tabdir,recursive=TRUE,showWarnings=FALSE)
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
  if (!dir.exists(tabdir)||!dir.exists(figdir)) stop("Cannot create publication folders")
  html_escape <- function(x) {
    x <- gsub("&","&amp;",as.character(x),fixed=TRUE)
    x <- gsub("<","&lt;",x,fixed=TRUE)
    gsub(">","&gt;",x,fixed=TRUE)
  }
  # V11.1: explicitly enumerate integer fields. NEVER match a trailing "n":
  # that inadvertently rounded Kappa_Hinselmann and lower_B_plugin to 0/1.
  # All values in the numeric CSV remain unmodified; HTML is presentation only.
  value_format <- function(x, colname="") {
    if (is.logical(x)) return(ifelse(is.na(x), "NA", ifelse(x,"TRUE","FALSE")))
    if (!is.numeric(x)) return(as.character(x))
    is_integer_field <- grepl(
      "^(n|.*_n|Participants|False_negatives|False_positives|requested|success|failures|events|warning_bearing_successes|retry_optimizer_successes|count|counts|replications)$",
      colname, ignore.case=TRUE)
    is_p_field <- grepl(
      "^(p|P_two_sided|p_two_sided|p_loss_two_sided|p_at_L|p_at_U|p_IUT|cf_p|reg_p|.*_pvalue|p_value)$",
      colname, ignore.case=TRUE)
    if (is_integer_field) return(ifelse(is.na(x),"NA",sprintf("%.0f",x)))
    if (is_p_field) return(ifelse(is.na(x),"NA",ifelse(x<.001,"<0.001",sprintf("%.4f",x))))
    # Precision-preserving output: no artificial zeros for small nonzero effects.
    out <- ifelse(abs(x)>0 & abs(x)<1e-5, sprintf("%.4e",x),sprintf("%.6f",x))
    ifelse(is.na(x),"NA",out)
  }
  # Executable display regression tests. No statistical estimation is changed.
  stopifnot(identical(value_format(.535466,"Kappa_Hinselmann"), "0.535466"),
            identical(value_format(-.00135143,"lower_B_plugin"), "-0.001351"),
            identical(value_format(.00099808,"highest_successful_gradient"), "0.000998"),
            identical(value_format(858,"n"), "858"))
  write_publication_table <- function(name,dat,caption,note="") {
    utils::write.csv(dat,file.path(tabdir,paste0(name,".csv")),
                     row.names=FALSE,fileEncoding="UTF-8",na="")
    contents <- lapply(seq_along(dat),function(j) value_format(dat[[j]],names(dat)[j]))
    pretty <- as.data.frame(contents,stringsAsFactors=FALSE)
    names(pretty) <- names(dat)
    headers <- paste0("<th>",html_escape(names(pretty)),"</th>",collapse="")
    rows <- apply(pretty,1,function(x)
      paste0("<tr>",paste0("<td>",html_escape(x),"</td>",collapse=""),"</tr>"))
    html <- c("<!DOCTYPE html><html><head><meta charset='UTF-8'>",
      "<style>body{font-family:Arial,'Times New Roman',sans-serif;color:#222;}",
      "table{border-collapse:collapse;width:100%;font-size:10pt;}",
      "th{text-align:left;border-top:1.5px solid #000;border-bottom:1px solid #000;padding:7px;}",
      "td{border-bottom:1px solid #ddd;padding:6px;}",
      "caption{text-align:left;font-weight:bold;font-size:11pt;padding-bottom:8px;}",
      ".note{font-size:9pt;margin-top:8px;}</style></head><body>",
      "<table><caption>",html_escape(caption),"</caption><thead><tr>",headers,
      "</tr></thead><tbody>",paste(rows,collapse="\n"),
      "</tbody></table><div class='note'>",html_escape(note),"</div></body></html>")
    writeLines(html,file.path(tabdir,paste0(name,".html")),useBytes=TRUE)
    invisible(dat)
  }
  fmt_ci <- function(v,lo,hi) sprintf("%.3f (%.3f, %.3f)",v,lo,hi)
  age_stats <- function(y) sprintf("%.1f (%.1f)",mean(y),stats::sd(y))
  n0 <- sum(dat$D==0); n1 <- sum(dat$D==1)
  tb1 <- data.frame(Characteristic=c("Participants, n","Age, mean (SD), y",
                                   "Age, median [IQR], y","Hinselmann positive, n (%)",
                                   "Schiller positive, n (%)"),
                    Overall=NA_character_,`Biopsy negative`=NA_character_,
                    `Biopsy positive`=NA_character_,check.names=FALSE)
  ss <- list(dat,dat[dat$D==0,,drop=FALSE],dat[dat$D==1,,drop=FALSE])
  for (i in seq_along(ss)) {
    a <- ss[[i]]
    tb1[[i+1]] <- c(as.character(nrow(a)),age_stats(a$X),
      sprintf("%.0f [%.0f–%.0f]",stats::median(a$X),
              stats::quantile(a$X,.25),stats::quantile(a$X,.75)),
      sprintf("%d (%.1f%%)",sum(a$T1),100*mean(a$T1)),
      sprintf("%d (%.1f%%)",sum(a$T2),100*mean(a$T2)))
  }
  write_publication_table("Table_1_Study_Characteristics",tb1,
    "Table 1. Characteristics of the cervical-cancer analysis sample",
    "UCI risk factors dataset; Age in years. Biopsy is the working reference status. Verification and sampling mechanisms are not independently established.")
  exact_prop <- function(success,total) {
    if (!total) return("Not estimable")
    ci <- stats::binom.test(success,total)$conf.int
    fmt_ci(success/total,ci[1],ci[2])
  }
  tt <- lapply(c("T1","T2"),function(v) {
    y <- dat[[v]]; d <- dat$D
    data.frame(Test=if (v=="T1") "Hinselmann" else "Schiller",
      `Sensitivity (exact 95% CI)`=exact_prop(sum(y==1 & d==1),sum(d==1)),
      `Specificity (exact 95% CI)`=exact_prop(sum(y==0 & d==0),sum(d==0)),
      `PPV (exact 95% CI)`=exact_prop(sum(y==1 & d==1),sum(y==1)),
      `NPV (exact 95% CI)`=exact_prop(sum(y==0 & d==0),sum(y==0)),
      False_negatives=sum(y==0 & d==1),False_positives=sum(y==1 & d==0),
      check.names=FALSE)
  })
  write_publication_table("Table_2_Unadjusted_Accuracy",do.call(rbind,tt),
    "Table 2. Unadjusted binary test accuracy against biopsy",
    "Intervals are exact binomial intervals for each metric individually; they are not paired comparative confidence intervals.")
  p3 <- primary[,c("policy","kappa1","kappa2","delta","se_delta",
             "ci_delta_low","ci_delta_high","p_two_sided",
             "risk_difference_A1_minus_A2")]
  names(p3) <- c("Assumed_loss_policy","Kappa_Hinselmann","Kappa_Schiller",
                 "Paired_difference","SE_difference","CI_low","CI_high",
                 "P_two_sided","Direct_loss_difference")
  write_publication_table("Table_3_Primary_Paired_Kappa",p3,
    "Table 3. Cross-fitted integrated loss-weighted kappas and paired contrasts",
    "CI: pointwise normal approximation using influence functions (multivariate Delta Method). Loss weights are hypothetical sensitivity schedules, not clinically elicited costs. Holm-adjusted tests are available in the raw results.")
  dd <- dominance
  cols <- intersect(c("tau","alpha_used","theta_one_step","se",
                  "lower_one_sided","reject_robust_null",
                  "plugin_reference_positivity","lower_B_plugin"),names(dd))
  write_publication_table("Table_4_Exploratory_Dominance",dd[,cols,drop=FALSE],
    "Table 4. Exploratory conservative robust-dominance sensitivity analysis",
    "Fixed tau soft-min test is conservative (rejection probability tends to zero at fixed null distributions). Plugin positivity is not a confidence bound. Do not claim certified dominance without validated inference and clinical loss justification.")
  write_publication_table("Table_5_Nuisance_Diagnostics",nuisance$calibration,
    "Table 5. Out-of-fold probability-model diagnostics",
    "Brier score, observed event frequency and mean cross-fitted prediction. These do not establish nuisance-model consistency.")
  if (!is.null(bs)) {
    write_publication_table("Table_6_Local_Bootstrap_Audit",bs$status,
      "Table 6. Subject-level local bootstrap computational audit",
      "The 95% successful-fit threshold is an operational guardrail, not proof of 95% confidence coverage. Bands are finite-grid and exploratory until Monte Carlo calibration.")
    if (nrow(bs$failure_messages))
      write_publication_table("Table_6b_Bootstrap_Failure_Reasons",bs$failure_messages,
        "Supplementary table. Reasons for bootstrap replication failures")
  }
  if (!is.null(benchmark))
    write_publication_table("Table_7_Model_Based_Benchmark",benchmark,
      "Supplementary table. Regression-standardized conditional benchmark",
      "The penalized local regression estimator is secondary; regularization changes finite-sample estimates.")
  report <- character()
  plot_export <- function(stem,draw,width=7.1,height=4.7) {
    for (type in c("pdf","tiff")) {
      target <- file.path(figdir,paste0(stem,".",type))
      succeeded <- tryCatch({
        if (type=="pdf") grDevices::pdf(target,width=width,height=height,
                                         useDingbats=FALSE)
        else grDevices::tiff(target,width=width,height=height,units="in",
                             res=600,compression="lzw")
        tryCatch(draw(),finally=grDevices::dev.off())
        TRUE
      },error=function(e) {
        report <<- c(report,paste(stem,type,conditionMessage(e),sep=" | "))
        FALSE
      })
      if (!succeeded && file.exists(target)) unlink(target)
    }
  }
  # Publication graphics intentionally omit decorative titles; captions below carry methods.
  plot_export("Figure_1_Age_Distribution",function() {
    graphics::hist(dat$X,breaks=seq(12,88,by=4),col="grey85",border="white",
                   main="",xlab="Age (years)",ylab="Participants",las=1)
    graphics::box(bty="l")
  })
  plot_export("Figure_2_Hypothetical_Loss_Policies",function() {
    xx <- seq(min(dat$X),max(dat$X),length.out=250)
    graphics::plot(xx,relative_fn_loss(xx,"age_dependent"),type="l",lwd=2,
      ylim=c(.3,.95),xlab="Age (years)",ylab="Relative false-negative loss c(x)",las=1)
    graphics::lines(xx,relative_fn_loss(xx,"equal"),lty=2,lwd=2)
    graphics::lines(xx,relative_fn_loss(xx,"fn_075"),lty=3,lwd=2)
    graphics::legend("topleft",c("Age-dependent (assumed)","Equal FN/FP", "FN weight 0.75"),
      lty=1:3,lwd=2,bty="n",cex=.82)
  })
  plot_export("Figure_3_Integrated_Kappa_Curves",function() {
    yr <- range(c(curve$kappa1,curve$kappa2),finite=TRUE)
    graphics::plot(curve$c_value,curve$kappa1,type="l",lwd=2,lty=1,
      ylim=yr,xlab="Constant false-negative loss weight c",
      ylab="Integrated loss-weighted kappa",las=1)
    graphics::lines(curve$c_value,curve$kappa2,lty=2,lwd=2)
    graphics::legend("bottomright",c("Hinselmann","Schiller"),lty=1:2,lwd=2,bty="n")
  })
  plot_export("Figure_4_Paired_Difference_Pointwise",function() {
    y <- range(c(curve$delta,curve$ci_delta_low,curve$ci_delta_high),finite=TRUE)
    graphics::plot(curve$c_value,curve$delta,type="n",ylim=y,
      xlab="Constant false-negative loss weight c",ylab="Paired kappa difference",las=1)
    graphics::polygon(c(curve$c_value,rev(curve$c_value)),
      c(curve$ci_delta_low,rev(curve$ci_delta_high)),
      col="grey90",border=NA)
    graphics::abline(h=0,lty=3)
    graphics::lines(curve$c_value,curve$delta,lwd=2)
    graphics::lines(curve$c_value,curve$ci_delta_low,lty=2)
    graphics::lines(curve$c_value,curve$ci_delta_high,lty=2)
    graphics::legend("topright",c("Paired contrast","Pointwise 95% CI"),
        lty=c(1,2),lwd=c(2,1),bty="n",cex=.85)
  })
  if (!is.null(local)) {
    plot_export("Figure_5_Conditional_Difference_Surface",function() {
      xx <- sort(unique(local$X)); cc <- sort(unique(local$c))
      # construct grid matrix by sorting within each loss weight.
      zz <- do.call(cbind,lapply(cc,function(z) {
        a <- local[abs(local$c-z)<1e-10,,drop=FALSE]
        a$delta[order(a$X)]
      }))
      graphics::image(xx,cc,zz,col=grDevices::hcl.colors(35,"Blue-Red 3"),
        xlab="Age (years)",ylab="Relative false-negative loss c",las=1)
      if (min(zz)<0 && max(zz)>0) graphics::contour(xx,cc,zz,levels=0,
                                               add=TRUE,drawlabels=FALSE,lwd=2)
    },width=6.8,height=5.0)
    plot_export("Figure_6_Age_Conditional_Contrast",function() {
      key <- CFG$local_c_values[which.min(abs(CFG$local_c_values-0.5))]
      a <- local[abs(local$c-key)<1e-10,,drop=FALSE]
      a <- a[order(a$X),,drop=FALSE]
      lim <- range(a$delta,finite=TRUE)
      bands <- NULL
      if (!is.null(bs) && !is.null(bs$results)) {
        bands <- bs$results[abs(bs$results$c-key)<1e-10,,drop=FALSE]
        bands <- bands[order(bands$X),,drop=FALSE]
        if (nrow(bands)==nrow(a))
          lim <- range(c(lim,bands$simultaneous_low,bands$simultaneous_high),finite=TRUE)
      }
      graphics::plot(a$X,a$delta,type="n",ylim=lim,xlab="Age (years)",
        ylab="Conditional paired kappa difference",las=1)
      graphics::abline(h=0,lty=3)
      if (!is.null(bands) && nrow(bands)==nrow(a)) {
        graphics::polygon(c(bands$X,rev(bands$X)),
          c(bands$simultaneous_low,rev(bands$simultaneous_high)),
          col="grey90",border=NA)
      }
      graphics::lines(a$X,a$delta,lwd=2)
      graphics::mtext(paste0("Fixed relative FN-loss c=",format(key),
        if (is.null(bands)) "; descriptive only" else "; exploratory finite-grid band"),
        side=3,cex=.8)
    })
  }
  if (length(report)) writeLines(report,file.path(figdir,"FIGURE_EXPORT_ERRORS.txt"))
  captions <- c(
    "Figure 1. Age distribution in the analytic sample (n=858).",
    "Figure 2. Prespecified hypothetical loss-weight schedules; not measured clinical preferences.",
    "Figure 3. Cross-fitted integrated weighted-kappa curves for two paired tests.",
    "Figure 4. Paired integrated kappa contrast with POINTWISE IF/Delta-method 95% Wald intervals; not a simultaneous band.",
    "Figure 5. Ridge-penalized descriptive conditional contrast on a finite age-by-loss grid.",
    "Figure 6. Age-conditional paired contrast at a fixed loss weight; bands only if bootstrap fit-screen passes; exploratory pending Monte Carlo.")
  writeLines(c("PUBLICATION DRAFT FIGURE CAPTIONS",captions,
     "Figures are not submitted as evidence of clinical loss elicitation or exact continuum coverage.",
     "The study must provide model diagnostics, ethically valid source description and MC calibration."),
     file.path(figdir,"Figure_Captions_and_Restrictions.txt"))
  writeLines(c("JBS TABLE PACKAGE: CSV + Word-openable HTML. All results auto-generated, never fabricated.",
    "Recheck values and citations before copying into a manuscript.",
    "The cross-fitted estimator requires nuisance convergence/rate assumptions not verified by this single dataset.",
    "Finite bootstrap fit success rate >95% is a computational guardrail only.",
    "If bootstrap bands were not accepted, do not display or report them as inferential results.",
    "No journal acceptance can be guaranteed by tables or plots."),
    file.path(tabdir,"TABLE_NOTES_READ_FIRST.txt"))
  invisible(list(tables=tabdir,figures=figdir))
}

# ==============================================================================
# V11 SUPPLEMENT: independent regression Delta/sandwich, bidirectional robust
# comparisons, finite-grid multiplier band, and independently generated Monte Carlo.
# This module is embedded in the single-file V11 code, not sourced externally.
# ==============================================================================

# A. Full joint sandwich covariance for THREE separate parametric mean regressions.
# The scores are STACKED PER SUBJECT. Their meat has the cross-test covariance.
# Unlike the algebraic IF-vs-Delta check, this is an INDEPENDENT regression estimator.
regression_delta_v11 <- function(dat, policies = c("equal", "fn_075", "age_dependent"),
                                 alpha = CFG$alpha) {
  n <- nrow(dat)
  Z <- (dat$X - mean(dat$X))/stats::sd(dat$X)
  XD <- cbind("Intercept"=1, "Z"=Z)
  XT <- cbind("Intercept"=1, "D"=dat$D, "Z"=Z, "D_Z"=dat$D*Z)
  designs <- list(XD, XT, XT)
  outcomes <- list(dat$D, dat$T1, dat$T2)
  models <- lapply(1:3, function(j) {
    fit <- tryCatch(suppressWarnings(stats::glm.fit(
      x = designs[[j]], y = outcomes[[j]], family=stats::binomial(),
      control=stats::glm.control(maxit=150L,epsilon=1e-9))),
      error=function(e)e)
    if (inherits(fit,"error") || !isTRUE(fit$converged) ||
        fit$rank != ncol(designs[[j]]) || any(!is.finite(fit$coefficients))) {
      stop("Regression Delta benchmark not estimable: separation, rank, or convergence in ",
           c("D|X","T1|D,X","T2|D,X")[j], ". This benchmark is not silently replaced.")
    }
    fit
  })
  beta <- unlist(lapply(models, `[[`, "coefficients"),use.names=FALSE)
  pred <- lapply(seq_len(3), function(j) stats::plogis(drop(designs[[j]] %*% models[[j]]$coefficients)))
  S <- do.call(cbind, Map(function(X,y,p) X*as.numeric(y-p), designs,outcomes,pred))
  # H = block diagonal mean information, normalized per observation.
  dims <- vapply(designs,ncol,integer(1))
  cuts <- cumsum(c(0L,dims))
  H <- matrix(0,length(beta),length(beta))
  for (j in 1:3) {
    idx <- (cuts[j]+1L):cuts[j+1L]
    pp <- pred[[j]]
    H[idx,idx] <- crossprod(designs[[j]], designs[[j]]*as.numeric(pp*(1-pp)))/n
  }
  if (!is.finite(kappa(H)) || kappa(H)>1e10) stop("Regression Delta benchmark: unstable information matrix")
  Hinv <- solve(H)
  # Predict under D=0 and D=1 from each Tj|D,X regression.
  X0 <- cbind(1,0,Z,0)
  X1 <- cbind(1,1,Z,Z)
  compute_components <- function(par, cv, sc) {
    bD <- par[1:2]; b1 <- par[3:6]; b2 <- par[7:10]
    p <- stats::plogis(drop(XD %*% bD))
    t10 <- stats::plogis(drop(X0 %*% b1)); t11 <- stats::plogis(drop(X1 %*% b1))
    t20 <- stats::plogis(drop(X0 %*% b2)); t21 <- stats::plogis(drop(X1 %*% b2))
    a1 <- sc*(cv*p*(1-t11)+(1-cv)*(1-p)*t10)
    a2 <- sc*(cv*p*(1-t21)+(1-cv)*(1-p)*t20)
    c1 <- sc*p*(1-p)*(t11-t10)
    c2 <- sc*p*(1-p)*(t21-t20)
    cbind(A1=a1,C1=c1,A2=a2,C2=c2)
  }
  out <- list()
  for (policy in policies) {
    cv <- relative_fn_loss(dat$X,policy); sc <- loss_scale(dat$X)
    per_person <- compute_components(beta,cv,sc)
    theta <- colMeans(per_person)
    R1 <- theta["A1"]+theta["C1"]; R2 <- theta["A2"]+theta["C2"]
    if (min(R1,R2)<=CFG$ref_denominator_tol) stop("Regression Delta reference loss <= 0")
    delta <- theta["C1"]/R1-theta["C2"]/R2
    # Numerical differentiation is only used for the smooth parametric mapping.
    G <- matrix(NA_real_,4,length(beta))
    for (m in seq_along(beta)) {
      step <- 1e-5*max(1,abs(beta[m]))
      plus <- minus <- beta; plus[m] <- plus[m]+step; minus[m] <- minus[m]-step
      G[,m] <- (colMeans(compute_components(plus,cv,sc))-
                  colMeans(compute_components(minus,cv,sc)))/(2*step)
    }
    # IF of integration over empirical X + IF of stacked regression coefficients.
    # Sandwich meat uses all 3 scores on the SAME subject; paired covariance retained.
    phi_t <- sweep(per_person,2,theta,"-") + S %*% t(Hinv) %*% t(G)
    grad <- c(-theta["C1"]/R1^2, theta["A1"]/R1^2,
              theta["C2"]/R2^2,-theta["A2"]/R2^2)
    phi_d <- drop(phi_t %*% grad)
    se <- stats::sd(phi_d)/sqrt(n)
    z <- if (se>0) delta/se else NA_real_
    zz <- stats::qnorm(1-alpha/2)
    out[[policy]] <- data.frame(policy=policy,method="regression_stacked_sandwich_delta",
      n=n,kappa1=theta["C1"]/R1,kappa2=theta["C2"]/R2,
      delta=unname(delta),se_delta=se,ci_low=delta-zz*se,
      ci_high=delta+zz*se,p_two_sided=if(is.finite(z)) 2*stats::pnorm(-abs(z)) else NA_real_,
      numerator_condition_number=kappa(H),n_parameters=length(beta),
      model_label="D~age linear; Tj~D*age linear; joint subject-level sandwich",
      stringsAsFactors=FALSE)
  }
  do.call(rbind,out)
}

# B. EXACT interval-sign property for CONSTANT loss weight c.
# IUT tests the two endpoint one-sided nulls; no multiple-c penalty needed for
# intersection-union reasoning, since rejecting requires BOTH rejections.
endpoint_iut_v11 <- function(dat, nuisance, cL=0.10, cU=0.90,
                             alpha=CFG$alpha) {
  if (length(cL)!=1L || length(cU)!=1L || !is.finite(cL+cU) ||
      cL<0 || cU>1 || cL>=cU) stop("Choose 0 <= cL < cU <= 1")
  ends <- lapply(c(cL,cU),function(cc)
    fixed_policy_inference(dat,nuisance,rep(cc,nrow(dat)),policy_name="endpoint")$summary)
  ss <- do.call(rbind,ends)
  z <- ss$delta/ss$se_delta
  if (any(!is.finite(z))) stop("Endpoint IUT has degenerate zero standard error")
  # A positive Delta supports T1, a negative Delta supports T2.
  p_t1 <- stats::pnorm(z,lower.tail=FALSE)
  p_t2 <- stats::pnorm(z)
  p1 <- max(p_t1); p2 <- max(p_t2)
  data.frame(direction=c("T1_over_T2","T2_over_T1"),
    range_L=cL,range_U=cU,delta_at_L=ss$delta[1],delta_at_U=ss$delta[2],
    se_at_L=ss$se_delta[1],se_at_U=ss$se_delta[2],
    p_at_L=c(p_t1[1],p_t2[1]),p_at_U=c(p_t1[2],p_t2[2]),
    p_IUT=c(p1,p2),reject_at_alpha=c(p1<alpha,p2<alpha),
    applies_to="ALL constant c between interval endpoints (positive reference losses)",
    stringsAsFactors=FALSE)
}

# C. Smooth upper envelope and its EIF, parallel to smooth lower envelope.
# b0 + E[max(Lh,Uh)] < 0 CERTIFIES T2 dominance across all functions in box.
# tau fixed: upper smooth bound >= exact upper bound, making the test conservative.
soft_max_v11 <- function(h,L,U,tau) {
  if (length(tau)!=1L || !is.finite(tau) || tau<=0) stop("tau must be positive")
  a <- L*h/tau; b <- U*h/tau; m <- pmax(a,b)
  tau*(m+log(exp(a-m)+exp(b-m)))
}
soft_max_gradient_v11 <- function(h,L,U,tau) {
  # weight on U is exp(Uh/tau)/[exp(Lh/tau)+exp(Uh/tau)]
  prob_U <- stats::plogis((U-L)*h/tau)
  L*(1-prob_U)+U*prob_U
}
upper_dominance_test_v11 <- function(env,tau,alpha=CFG$alpha) {
  n <- length(env$p)
  g <- soft_max_v11(env$h,env$L,env$U,tau)
  w <- soft_max_gradient_v11(env$h,env$L,env$U,tau)
  aa1 <- mean(w*env$s*(env$p-env$q2))
  aa2 <- mean(w*env$s*(env$p-env$q1))
  theta <- env$b0+mean(g+w*env$s*env$xi)
  phi <- (env$V2+aa1)*env$phiC1-(env$V1+aa2)*env$phiC2+
    env$C1*env$phiV2-env$C2*env$phiV1+g-mean(g)+w*env$s*env$xi
  phi <- phi-mean(phi)
  se <- stats::sd(phi)/sqrt(n)
  ub <- theta+stats::qnorm(1-alpha)*se
  eligible <- all(is.finite(c(env$minR1_plugin,env$minR2_plugin))) &&
    min(env$minR1_plugin,env$minR2_plugin)>CFG$ref_denominator_tol
  list(summary=data.frame(direction="T2_over_T1",tau=tau,
    upper_B_plugin=env$upper_plugin,theta_upper_smooth=theta,
    se=se,upper_one_sided=ub,alpha_used=alpha,
    reject_robust_null=if(eligible) as.integer(is.finite(ub)&&ub<0) else NA_integer_,
    reference_positivity_plugin=eligible,
    smoothing_bound_gap=tau*log(2),
    qualification="pointwise fixed-tau conservative; needs nuisance rate and MC verification",
    stringsAsFactors=FALSE),influence=phi)
}
robust_two_direction_v11 <- function(env,taus=CFG$theta_taus,alpha=CFG$alpha) {
  if (any(!is.finite(taus)) || any(taus<=0) || anyDuplicated(taus))
    stop("taus must be distinct, finite, positive and prespecified")
  # Family-wise Bonferroni over tau choices AND two directions.
  a_each <- alpha/(2*length(taus))
  ret <- lapply(taus,function(tau) {
    low <- smooth_dominance_test(env,tau,alpha=a_each)$summary
    hi <- upper_dominance_test_v11(env,tau,alpha=a_each)$summary
    a <- data.frame(direction="T1_over_T2",tau=tau,
      bound_plugin=low$lower_B_plugin,theta_smooth=low$theta_one_step,
      se=low$se,one_sided_bound=low$lower_one_sided,
      alpha_used=a_each,reject=low$reject_robust_null,
      bound_type="lower",pointwise_test="LCB(theta_lower_tau) > 0",
      stringsAsFactors=FALSE)
    b <- data.frame(direction="T2_over_T1",tau=tau,
      bound_plugin=hi$upper_B_plugin,theta_smooth=hi$theta_upper_smooth,
      se=hi$se,one_sided_bound=hi$upper_one_sided,
      alpha_used=a_each,reject=hi$reject_robust_null,
      bound_type="upper",pointwise_test="UCB(theta_upper_tau) < 0",
      stringsAsFactors=FALSE)
    rbind(a,b)
  })
  ans <- do.call(rbind,ret); rownames(ans) <- NULL
  ans
}

# D. FINITE-GRID simultaneous bands for the INTEGRATED constant-c contrast.
# Multiplier bootstrap of influence functions reuses frozen cross-fitted nuisance
# predictions. This is a fast first-order method; NOT exact continuum coverage.
multiplier_integrated_band_v11 <- function(dat,nuisance,
        grid=seq(0.05,0.95,by=.025),B=1999L,
        seed=CFG$seed+111L,alpha=CFG$alpha) {
  if (B<199L || length(grid)<3L || any(diff(grid)<=0))
    stop("Prespecify an increasing c-grid and at least 199 multiplier draws")
  n <- nrow(dat)
  obj <- lapply(grid,function(cv) fixed_policy_inference(
    dat,nuisance,rep(cv,n),policy_name="fixed_weight"))
  mat <- do.call(cbind,lapply(obj,function(o)o$influence$phi_delta))
  se_cols <- apply(mat,2,stats::sd)
  if (any(!is.finite(se_cols)) || any(se_cols<1e-10))
    stop("Multiplier band: nonpositive influence-function variance on grid")
  set.seed(seed)
  # Matrix of standard normal multipliers (one vector per draw).
  crit_max <- numeric(B)
  block <- 250L
  for (start in seq.int(1L,B,by=block)) {
    m <- min(block,B-start+1L)
    noise <- matrix(stats::rnorm(n*m),nrow=n,ncol=m)
    numer <- crossprod(noise,mat)
    standardized <- sweep(numer,2,sqrt(n)*se_cols,"/")
    crit_max[start:(start+m-1L)] <- apply(abs(standardized),1,max)
  }
  critical <- unname(stats::quantile(crit_max,1-alpha,names=FALSE))
  vv <- vapply(obj,function(o)o$summary$delta,numeric(1))
  se <- se_cols/sqrt(n)
  data.frame(c=grid,delta=vv,se=se,
     pointwise_low=vv-stats::qnorm(1-alpha/2)*se,
     pointwise_high=vv+stats::qnorm(1-alpha/2)*se,
     simultaneous_grid_low=vv-critical*se,
     simultaneous_grid_high=vv+critical*se,
     critical_max_t=critical,multiplier_draws=B,
     interpretation="finite prespecified c grid; asymptotic, subject-level multiplier IF",
     stringsAsFactors=FALSE)
}

# E. True continuous X (truncated normal), paired test outcomes via common random U.
# Null has equal test sensitivity/specificity functions, not identical realizations.
# Strong alternatives and a nonlinear misspecification scenario are included.
mc_truth_functions_v11 <- function(x,scenario,shape) {
  z <- (x-30)/12
  nonlinear <- if(shape=="nonlinear") 0.40*sin(z*1.4) else rep(0,length(z))
  p <- stats::plogis(-1.7+.35*z+nonlinear)
  se1 <- stats::plogis(.8+.20*z + .4*nonlinear)
  sp1 <- stats::plogis(2.2-.13*z -.2*nonlinear)
  if (scenario=="null") {se2<-se1;sp2<-sp1}
  else if (scenario=="t2_local") {
    # V11.1 moderate alternative: stress power near the null, instead of relying
    # solely on almost-certain-power strong alternatives. Not a null scenario.
    se2 <- stats::plogis(.95+.20*z+.4*nonlinear)
    sp2 <- stats::plogis(2.35-.13*z-.2*nonlinear)
  }
  else if (scenario=="t2_strong") {
    se2 <- stats::plogis(2.2+.20*z+.4*nonlinear)
    sp2 <- stats::plogis(4.0-.13*z-.2*nonlinear)
  } else if(scenario=="t1_strong") {
    se2 <- stats::plogis(-.1+.20*z+.4*nonlinear)
    sp2 <- stats::plogis(1.0-.13*z-.2*nonlinear)
  } else if(scenario=="tradeoff") {
    se2 <- stats::plogis(1.9+.20*z+.4*nonlinear)
    sp2 <- stats::plogis(1.5-.13*z-.2*nonlinear)
  } else stop("Unknown simulation scenario")
  list(p=p,se1=se1,sp1=sp1,se2=se2,sp2=sp2)
}
mc_generate_v11 <- function(n,scenario="null",shape="linear",dependence=.55) {
  low <- stats::pnorm(13,mean=30,sd=12)
  high <- stats::pnorm(84,mean=30,sd=12)
  X <- stats::qnorm(stats::runif(n,low,high),mean=30,sd=12)
  f <- mc_truth_functions_v11(X,scenario,shape)
  D <- stats::rbinom(n,1,f$p)
  rate1 <- ifelse(D==1,f$se1,1-f$sp1)
  rate2 <- ifelse(D==1,f$se2,1-f$sp2)
  u1 <- stats::runif(n)
  u2 <- ifelse(stats::runif(n)<dependence,u1,stats::runif(n))
  data.frame(X=X,D=D,T1=as.integer(u1<rate1),T2=as.integer(u2<rate2))
}
mc_truth_v11 <- function(scenario="null",shape="linear",policy="age_dependent",
                         m=50001L) {
  low <- stats::pnorm(13,mean=30,sd=12)
  high <- stats::pnorm(84,mean=30,sd=12)
  X <- stats::qnorm(low+(high-low)*((seq_len(m)-0.5)/m),mean=30,sd=12)
  ff <- mc_truth_functions_v11(X,scenario,shape)
  p <- ff$p
  q1 <- p*ff$se1+(1-p)*(1-ff$sp1)
  q2 <- p*ff$se2+(1-p)*(1-ff$sp2)
  s <- loss_scale(X); cval <- relative_fn_loss(X,policy)
  C1 <- mean(s*p*(1-p)*(ff$se1+ff$sp1-1))
  C2 <- mean(s*p*(1-p)*(ff$se2+ff$sp2-1))
  A1 <- mean(s*(cval*p*(1-ff$se1)+(1-cval)*(1-p)*(1-ff$sp1)))
  A2 <- mean(s*(cval*p*(1-ff$se2)+(1-cval)*(1-p)*(1-ff$sp2)))
  V1 <- mean(s*(1-p)*q1); V2 <- mean(s*(1-p)*q2)
  h <- s*(C1*(p-q2)-C2*(p-q1))
  b0 <- C1*V2-C2*V1
  L <- box_L(X); U <- box_U(X)
  list(delta=C1/(A1+C1)-C2/(A2+C2),
       lower=b0+mean(pmin(L*h,U*h)),upper=b0+mean(pmax(L*h,U*h)),
       ratio1=C1/(A1+C1),ratio2=C2/(A2+C2))
}

# F. Independent MC comparisons of Cross-fitting vs REGRESSION DELTA and both
# dominance directions; report estimator failures, MCSE, and whether the null
# is TRUE at the chosen data-generating mechanism.
run_mc_v11 <- function(n=858L,reps=500L,scenario="null",shape="linear",
                       policy="age_dependent",tau=.0005,
                       seed=20261201L,dependence=.55) {
  if(!is.finite(dependence) || dependence<0 || dependence>1)
    stop("dependence must belong to [0,1]")
  if (reps<20L || n<100L) stop("Use >=20 MC replicates, n >=100")
  tr <- mc_truth_v11(scenario,shape,policy)
  ans <- vector("list",reps)
  for (r in seq_len(reps)) {
    set.seed(seed+as.integer(r)*1009L+as.integer(n)*17L)
    dat <- mc_generate_v11(n,scenario,shape,dependence=dependence)
    one <- tryCatch({
      pred <- crossfit_probabilities(dat,K=CFG$folds,method=CFG$nuisance_model,
                                    df=CFG$nuisance_df,seed=seed+r)
      f <- fixed_policy_inference(dat,pred,relative_fn_loss(dat$X,policy),policy_name=policy)
      env <- prepare_envelope(dat,pred)
      lo <- smooth_dominance_test(env,tau,alpha=CFG$alpha/2)$summary
      hi <- upper_dominance_test_v11(env,tau,alpha=CFG$alpha/2)$summary
      # Regression comparator failure is recorded INDEPENDENTLY from primary.
      bm <- tryCatch(regression_delta_v11(dat,policies=policy),error=function(e)e)
      data.frame(replication=r,primary_ok=TRUE,regression_ok=!inherits(bm,"error"),
        cf_delta=f$summary$delta,cf_se=f$summary$se_delta,
        cf_ci_low=f$summary$ci_delta_low,cf_ci_high=f$summary$ci_delta_high,
        cf_reject=as.integer(f$summary$p_two_sided<CFG$alpha),
        reg_delta=if(!inherits(bm,"error")) bm$delta else NA_real_,
        reg_se=if(!inherits(bm,"error")) bm$se_delta else NA_real_,
        reg_low=if(!inherits(bm,"error")) bm$ci_low else NA_real_,
        reg_high=if(!inherits(bm,"error")) bm$ci_high else NA_real_,
        reg_reject=if(!inherits(bm,"error")) as.integer(bm$p_two_sided<CFG$alpha) else NA_integer_,
        theta_lower_hat=lo$theta_one_step,theta_upper_hat=hi$theta_upper_smooth,
        lower_bound_hat=lo$lower_one_sided,upper_bound_hat=hi$upper_one_sided,
        t1_certified=lo$reject_robust_null,t2_certified=hi$reject_robust_null,
        error_message=if(inherits(bm,"error")) paste("regression:",conditionMessage(bm)) else "",
        stringsAsFactors=FALSE)
    },error=function(e)e)
    if (inherits(one,"error"))
      ans[[r]] <- data.frame(replication=r,primary_ok=FALSE,regression_ok=FALSE,
        cf_delta=NA_real_,cf_se=NA_real_,cf_ci_low=NA_real_,cf_ci_high=NA_real_,
        cf_reject=NA_integer_,reg_delta=NA_real_,reg_se=NA_real_,
        reg_low=NA_real_,reg_high=NA_real_,reg_reject=NA_integer_,
        theta_lower_hat=NA_real_,theta_upper_hat=NA_real_,
        lower_bound_hat=NA_real_,upper_bound_hat=NA_real_,
        t1_certified=NA_integer_,t2_certified=NA_integer_,
        error_message=conditionMessage(one),stringsAsFactors=FALSE)
    else ans[[r]] <- one
    if(r%%50L==0L) logline("V11 MC: ",scenario,"/",shape," n=",n," rep ",r,"/",reps)
  }
  raw <- do.call(rbind,ans)
  good <- raw[raw$primary_ok,,drop=FALSE]
  reg <- raw[raw$regression_ok,,drop=FALSE]
  summarize_method <- function(one,which) {
    if(!nrow(one))return(data.frame(method=which,successful=0L,bias=NA_real_,
      rmse=NA_real_,empirical_sd=NA_real_,mean_se=NA_real_,
      coverage=NA_real_,coverage_mcse=NA_real_,rejection_rate=NA_real_,
      rejection_mcse=NA_real_,mean_width=NA_real_))
    if(which=="crossfit") {delta<-one$cf_delta;se<-one$cf_se;
      lower<-one$cf_ci_low;upper<-one$cf_ci_high;reject<-one$cf_reject}
    else {delta<-one$reg_delta;se<-one$reg_se;
      lower<-one$reg_low;upper<-one$reg_high;reject<-one$reg_reject}
    cover <- mean(lower<=tr$delta&upper>=tr$delta)
    rr <- mean(reject)
    data.frame(method=which,successful=nrow(one),bias=mean(delta)-tr$delta,
      rmse=sqrt(mean((delta-tr$delta)^2)),empirical_sd=stats::sd(delta),
      mean_se=mean(se),coverage=cover,
      coverage_mcse=sqrt(cover*(1-cover)/nrow(one)),
      rejection_rate=rr,rejection_mcse=sqrt(rr*(1-rr)/nrow(one)),
      mean_width=mean(upper-lower))
  }
  tab <- rbind(summarize_method(good,"crossfit"),summarize_method(reg,"regression_stacked_delta"))
  tab$n <- n;tab$scenario <- scenario;tab$shape <- shape;tab$policy <- policy
  tab$dependence <- dependence
  tab$truth_delta <- tr$delta;tab$delta_null_true <- abs(tr$delta)<1e-9
  tab$total_reps <- reps;tab$primary_failure_fraction <- 1-nrow(good)/reps
  tab$regression_failure_fraction <- 1-nrow(reg)/reps
  dom <- data.frame(scenario=scenario,shape=shape,n=n,policy=policy,
    dependence=dependence,replications=reps,valid_reps=nrow(good),
    truth_lower_B=tr$lower,truth_upper_B=tr$upper,
    null_T1_no_dominance=(tr$lower<=0),
    null_T2_no_dominance=(tr$upper>=0),
    T1_certification_fraction=if(nrow(good))mean(good$t1_certified,na.rm=TRUE) else NA_real_,
    T2_certification_fraction=if(nrow(good))mean(good$t2_certified,na.rm=TRUE) else NA_real_,
    inference="fixed tau; 2-direction Bonferroni; not uniform in distribution")
  list(summary=tab,dominance=dom,raw=raw,truth=tr)
}

# G0. INDEPENDENT LOCAL GRID COVERAGE CALIBRATION (computationally expensive).
# This tests the exact local ridge/bootstrap workflow under a KNOWN continuous-X
# DGM. A failed fit or rejected bootstrap band is NOT counted as a success.
mc_true_local_delta_v11 <- function(x_grid,c_grid,scenario,shape) {
  ff <- mc_truth_functions_v11(rep(x_grid,each=length(c_grid)),scenario,shape)
  xx <- rep(x_grid,each=length(c_grid)); cc<-rep(c_grid,times=length(x_grid))
  p<-ff$p
  C1<-p*(1-p)*(ff$se1+ff$sp1-1)
  C2<-p*(1-p)*(ff$se2+ff$sp2-1)
  A1<-cc*p*(1-ff$se1)+(1-cc)*(1-p)*(1-ff$sp1)
  A2<-cc*p*(1-ff$se2)+(1-cc)*(1-p)*(1-ff$sp2)
  data.frame(X=xx,c=cc,delta=C1/(A1+C1)-C2/(A2+C2))
}
run_local_grid_calibration_v11 <- function(n=CFG$v11_local_mc_n,
    reps=CFG$v11_local_mc_reps,B=CFG$v11_local_mc_bootstrap_B,
    scenario="null",shape="linear",seed=20265001L,
    x_grid=c(20,25,30,35,40),c_grid=c(.25,.50,.75)) {
  if (reps<20L || B<50L || n<200L) stop("Local calibration: n >= 200, reps >=20, B >=50")
  truth<-mc_true_local_delta_v11(x_grid,c_grid,scenario,shape)
  log<-vector("list",reps)
  for(k in seq_len(reps)) {
    set.seed(seed+1019L*k+n)
    dat<-mc_generate_v11(n,scenario,shape)
    one<-tryCatch({
      fit<-fit_marginal_spline(dat)
      boot<-bootstrap_joint_local(dat,fit,x_grid,c_grid,B=B,seed=seed+k)
      if(!isTRUE(boot$status$accepted_for_finite_grid_CI) || is.null(boot$results)) {
        data.frame(replication=k,accepted=FALSE,coverage=NA,
          bootstrap_success_rate=boot$status$success_rate,
          reason="bootstrap bands rejected by numerical fit screen")
      } else {
        r<-boot$results
        loc <- match(paste(truth$X,truth$c),paste(r$X,r$c))
        if(anyNA(loc))stop("Local grid coordinates mismatched the truth grid")
        covered <- all(r$simultaneous_low[loc]<=truth$delta &
                       r$simultaneous_high[loc]>=truth$delta)
        data.frame(replication=k,accepted=TRUE,coverage=covered,
          bootstrap_success_rate=boot$status$success_rate,reason="")
      }
    },error=function(e)data.frame(replication=k,accepted=FALSE,coverage=NA,
      bootstrap_success_rate=NA_real_,reason=conditionMessage(e)))
    log[[k]]<-one
    if(k%%20L==0L)logline("Local finite-grid MC ",scenario,"/",shape,": ",k,"/",reps)
  }
  raw<-do.call(rbind,log)
  num_ok<-sum(raw$accepted)
  success_cov <- if(num_ok)mean(raw$coverage[raw$accepted])else NA_real_
  summary<-data.frame(n=n,reps=reps,B=B,scenario=scenario,shape=shape,
    x_min=min(x_grid),x_max=max(x_grid),grid_points=nrow(truth),
    accepted_fraction=num_ok/reps,
    simultaneous_coverage_given_accepted=success_cov,
    mcse_coverage_given_accepted=if(num_ok) sqrt(success_cov*(1-success_cov)/num_ok) else NA_real_,
    unconditional_coverage_and_accepted=sum(raw$accepted & (raw$coverage %in% TRUE))/reps,
    note="finite prespecified grid only; rejected fits count as no inference")
  list(summary=summary,raw=raw,truth=truth)
}

# G. Tables, diagnostics, captions, new figures and code self-tests.
plot_v11_pdf_tiff <- function(base,fun,width=6.8,height=4.8) {
  for(ext in c("pdf","tiff")){
    filepath<-paste0(base,".",ext)
    if(ext=="pdf") grDevices::pdf(filepath,width=width,height=height,
                                   family="Helvetica",useDingbats=FALSE)
    else grDevices::tiff(filepath,width=width,height=height,units="in",
                         res=600,compression="lzw")
    tryCatch(fun(),finally=grDevices::dev.off())
  }
}
run_self_tests_v11 <- function() {
  stopifnot(isTRUE(run_self_tests_v10()))
  hh <- c(-3,-.2,0,.4,2); LL<-rep(.4,5);UU<-rep(.9,5);tau<-.03
  rawmin<-pmin(LL*hh,UU*hh);rawmax<-pmax(LL*hh,UU*hh)
  stopifnot(all(soft_min(hh,LL,UU,tau) <= rawmin+1e-12),
            all(rawmin-soft_min(hh,LL,UU,tau)<=tau*log(2)+1e-12),
            all(soft_max_v11(hh,LL,UU,tau)>=rawmax-1e-12),
            all(soft_max_v11(hh,LL,UU,tau)-rawmax<=tau*log(2)+1e-12))
  # Numerical derivatives of both smooth operators.
  eps<-1e-6
  check_upper <- (soft_max_v11(hh+eps,LL,UU,tau)-soft_max_v11(hh-eps,LL,UU,tau))/(2*eps)
  check_lower <- (soft_min(hh+eps,LL,UU,tau)-soft_min(hh-eps,LL,UU,tau))/(2*eps)
  stopifnot(max(abs(check_upper-soft_max_gradient_v11(hh,LL,UU,tau)))<1e-6,
            max(abs(check_lower-soft_min_grad(hh,LL,UU,tau)))<1e-6)
  # Finite synthetic checks for separately fitted stacked GLM benchmark.
  set.seed(44342L)
  xd <- stats::runif(450,18,55)
  d <- stats::rbinom(450,1,stats::plogis(-1+.03*(xd-35)))
  tt1<-stats::rbinom(450,1,stats::plogis(-2+2*d+.015*(xd-35)))
  tt2<-stats::rbinom(450,1,stats::plogis(-1.5+1.7*d+.01*(xd-35)))
  tmp <- regression_delta_v11(data.frame(X=xd,D=d,T1=tt1,T2=tt2),policies="equal")
  stopifnot(nrow(tmp)==1,all(is.finite(tmp$delta)),tmp$se_delta>0)
  invisible(TRUE)
}

v11_extensions <- function(dat,pred,main_primary,curve,env) {
  logline("V11 add-on START: regression sandwich Delta; two-direction robust tests; multiplier band")
  if (!isTRUE(run_self_tests_v11())) stop("V11 unit tests failed")
  output <- CFG$out_dir
  table_dir<-file.path(output,"JBS_Tables")
  figure_dir<-file.path(output,"JBS_Figures")
  dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
  dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
  save_both <- function(df,name) {
    write_csv(df,paste0("V11_",name,".csv"))
    utils::write.csv(df,file.path(table_dir,paste0("V11_",name,".csv")),row.names=FALSE)
    # Standalone HTML table for import into Word; numeric CSV remains the source of truth.
    escape_html <- function(z) {
      z <- gsub("&","&amp;",as.character(z),fixed=TRUE)
      z <- gsub("<","&lt;",z,fixed=TRUE)
      gsub(">","&gt;",z,fixed=TRUE)
    }
    head <- paste0("<th>",escape_html(names(df)),"</th>",collapse="")
    rows <- if(nrow(df)) vapply(seq_len(nrow(df)),function(i){
      vals <- vapply(df[i,,drop=FALSE],function(z){
        if(is.numeric(z) && !is.na(z)) formatC(z,format="fg",digits=5)
        else as.character(z)
      },character(1))
      paste0("<tr>",paste0("<td>",escape_html(vals),"</td>",collapse=""),"</tr>")
    },character(1)) else character(0)
    lines <- c("<!DOCTYPE html><html><head><meta charset='utf-8'>",
      "<style>body{font:11pt Arial;color:#222}table{border-collapse:collapse}",
      "th{border-top:2px solid black;border-bottom:1px solid black;padding:7px}",
      "td{border-bottom:1px solid #ddd;padding:6px}</style></head><body>",
      paste0("<h3>",escape_html(name),"</h3><table><thead><tr>",head,"</tr></thead><tbody>"),
      rows,"</tbody></table>",
      "<p><em>Preliminary analytic output; verify and annotate before submission.</em></p>",
      "</body></html>")
    writeLines(lines,file.path(table_dir,paste0("V11_",name,".html")),useBytes=TRUE)
  }
  regression <- tryCatch(regression_delta_v11(dat),error=function(e)e)
  if(inherits(regression,"error")){
    writeLines(conditionMessage(regression),file.path(output,"V11_REGRESSION_DELTA_FAILURE.txt"))
    warning("Independent regression Delta benchmark failed; see V11_REGRESSION_DELTA_FAILURE.txt")
  } else save_both(regression,"Table8_Regression_Delta_Sandwich")
  iut <- endpoint_iut_v11(dat,pred,cL=CFG$endpoint_interval[1],
                            cU=CFG$endpoint_interval[2])
  save_both(iut,"Table9_Constant_Loss_Endpoint_IUT")
  both <- robust_two_direction_v11(env,taus=CFG$theta_taus)
  save_both(both,"Table10_Bidirectional_Functional_Dominance")
  mult <- multiplier_integrated_band_v11(dat,pred,
      grid=seq(.05,.95,by=.025),B=CFG$multiplier_B)
  save_both(mult,"Table11_Integrated_Contrast_Multiplier_Grid_Band")
  compare <- data.frame(policy=main_primary$policy,
       cf_delta=main_primary$delta,cf_se=main_primary$se_delta,
       cf_lower=main_primary$ci_delta_low,cf_upper=main_primary$ci_delta_high,
       cf_p=main_primary$p_two_sided)
  if(!inherits(regression,"error")){
    ord <- match(compare$policy,regression$policy)
    compare$reg_delta <- regression$delta[ord]
    compare$reg_se <- regression$se_delta[ord]
    compare$reg_lower<-regression$ci_low[ord]
    compare$reg_upper<-regression$ci_high[ord]
    compare$reg_p<-regression$p_two_sided[ord]
    compare$delta_estimator_difference <- compare$cf_delta-compare$reg_delta
  }
  save_both(compare,"Table12_Primary_Versus_Regression_Benchmark")
  plot_v11_pdf_tiff(file.path(figure_dir,"V11_Figure7_Simultaneous_Integrated_Grid_Band"),function(){
    graphics::plot(mult$c,mult$delta,type="n",ylim=range(mult$simultaneous_grid_low,
       mult$simultaneous_grid_high),xlab="Relative false-negative loss c",
       ylab="Hinselmann minus Schiller: weighted kappa")
    graphics::polygon(c(mult$c,rev(mult$c)),
      c(mult$simultaneous_grid_low,rev(mult$simultaneous_grid_high)),
      col="grey85",border=NA)
    graphics::abline(h=0,lty=3)
    graphics::lines(mult$c,mult$delta,lwd=2)
    graphics::lines(mult$c,mult$pointwise_low,lty=2)
    graphics::lines(mult$c,mult$pointwise_high,lty=2)
    graphics::legend("topright",c("Paired contrast","95% finite-grid simultaneous band",
      "Pointwise 95% Wald"),lty=c(1,1,2),lwd=c(2,8,1),
      col=c("black","grey85","black"),bty="n",cex=.75)
  })
  if(!inherits(regression,"error")) plot_v11_pdf_tiff(
      file.path(figure_dir,"V11_Figure8_Estimator_Comparison"),function(){
    xx<-seq_along(compare$policy)
    rng<-range(c(compare$cf_lower,compare$cf_upper,
                 compare$reg_lower,compare$reg_upper),finite=TRUE)
    graphics::plot(xx,compare$cf_delta,pch=16,ylim=rng,xlim=c(.5,length(xx)+.5),
      xaxt="n",xlab="Prespecified loss policy",ylab="Paired kappa contrast")
    graphics::axis(1,at=xx,labels=compare$policy)
    graphics::abline(h=0,lty=3)
    graphics::arrows(xx-.07,compare$cf_lower,xx-.07,compare$cf_upper,
                     angle=90,code=3,length=.05)
    graphics::points(xx-.07,compare$cf_delta,pch=16)
    graphics::arrows(xx+.07,compare$reg_lower,xx+.07,compare$reg_upper,
                     angle=90,code=3,length=.05,lty=2)
    graphics::points(xx+.07,compare$reg_delta,pch=1)
    graphics::legend("topright",c("Cross-fitted EIF/Delta","Regression sandwich Delta"),
                     pch=c(16,1),lty=c(1,2),bty="n",cex=.8)
  })
  notes <- c(
    "V11 JOURNAL-REVIEW INTERPRETATION RULES",
    "1. No code, table, simulation, or image guarantees acceptance in Journal of Biopharmaceutical Statistics.",
    "2. Data are from UCI observational collection; selection/verification mechanism must be justified; Biopsy is treated as recorded reference.",
    "3. Age-dependent c(X) is a hypothetical sensitivity-analysis loss preference; it was NOT elicited from clinicians.",
    "3b. Recorded age is integer-valued in years: numeric-scale smoothing is reasonable but literal absolutely continuous-X assumptions do NOT automatically hold.",
    "4. Regression Delta benchmark fits D|age and each Tj|D,age with a joint subject-level score sandwich; dependence between tests is retained.",
    "5. Its regression benchmark is a distinct estimator, and can be biased if conditional logistic mean models are misspecified.",
    "6. Cross-fitting with fixed spline df does NOT prove nuisance consistency over arbitrary data-generating laws.",
    "7. Pointwise inference relies on nuisance convergence and crossfit product-rate conditions; nominal coverage needs MC validation.",
    "8. Functional loss-box robust tests use fixed tau smooth LOWER and UPPER bounds; each test is conservative and may have low power.",
    "9. Bonferroni for two directions x all prespecified taus; positivity checks are PLUG-IN, not an inferential positivity certificate.",
    "10. Endpoint IUT is for ALL CONSTANT loss weights in [cL,cU] with positive R; it is NOT inference across arbitrary functions c(X).",
    "11. Multiplier band is simultaneous over a prespecified FINITE GRID of constant c only, not exact continuum coverage.",
    "12. Ridge local bootstrap is exploratory: convergence does not establish 95% simultaneous coverage; validate via independent MC.",
    "13. Missing data/verification bias, diagnostics in rare diseased strata, policy elicitation, reproducibility and MCSE must be disclosed.",
    "14. Final manuscript: benchmark against unadjusted approaches, show paired clinical risk contrast, prespecification and simulation failures.",
    paste0("R session: ",R.version.string),
    paste0("Analysis date UTC: ",format(Sys.time(),tz="UTC",usetz=TRUE)))
  writeLines(notes,file.path(output,"V11_STATISTICAL_DISCLOSURES_READ_FIRST.txt"))
  writeLines(c("Additional V11 manuscript captions:",
    "Figure 7. Cross-fitted integrated paired weighted-kappa difference under constant loss: asymptotic finite-grid multiplier simultaneous confidence band (gray) and pointwise Wald boundaries (dashed). Bands require Monte Carlo coverage validation.",
    "Figure 8. Prespecified loss policies: comparison of cross-fitted influence-function/Delta and separate stacked-score regression sandwich-Delta estimates with pointwise 95% intervals.",
    "Table 8. Independent regression Delta benchmark retaining paired score covariance.",
    "Table 9. Endpoint intersection-union tests for dominance across all constant loss weights in prespecified interval.",
    "Table 10. Both-direction conservative functional-loss dominance tests with fixed-tau smoothing and multiplicity correction.",
    "Table 11. Finite-grid integrated-contrast multiplier band; continuum uniform coverage not claimed.",
    "Table 12. Empirical comparison of estimation frameworks."),
    file.path(output,"V11_MANUSCRIPT_TABLE_FIGURE_NOTES.txt"))
  if (isTRUE(CFG$run_v11_simulation)) {
    if (CFG$v11_simulation_reps < 500L) warning(
       "MC <500 per cell is PILOT ONLY; do not infer nominal 5% rejection accurately")
    sims <- list(); doms <- list()
    for(nn in CFG$v11_simulation_ns)
      for(sc in CFG$v11_simulation_scenarios)
        for(shape in CFG$v11_simulation_shapes) {
          name<-paste(sc,shape,nn,sep="_")
          vv <- tryCatch(run_mc_v11(n=nn,reps=CFG$v11_simulation_reps,
            scenario=sc,shape=shape,policy="age_dependent",
            tau=CFG$theta_taus[2],seed=CFG$seed+nn),error=function(e)e)
          if(inherits(vv,"error")) {
            writeLines(conditionMessage(vv),file.path(output,paste0("V11_MC_ERROR_",name,".txt")))
          } else {
            sims[[name]]<-vv$summary; doms[[name]]<-vv$dominance
            write_csv(vv$raw,paste0("V11_MC_RAW_",name,".csv"))
          }
        }
    if(length(sims)) {
      mc_sum <- do.call(rbind,sims)
      save_both(mc_sum,"Table13_MC_Bias_RMSE_Coverage_TypeI_Power")
      save_both(do.call(rbind,doms),"Table14_MC_Both_Direction_Dominance")
      # Grid search using a conservative estimated lower bound for power, not
      # a claimed closed-form n formula. Missing cells yield no certification.
      aa <- mc_sum[mc_sum$scenario %in% c("t1_strong","t2_strong") &
                   mc_sum$method=="crossfit",,drop=FALSE]
      if(nrow(aa)) {
        aa$power_lower_normal <- pmax(0,aa$rejection_rate-1.96*aa$rejection_mcse)
        aa$target_power <- .80
        aa$qualifies <- (aa$power_lower_normal>=aa$target_power) &
          (aa$primary_failure_fraction<=.05)
        aa$min_qualifying_n <- NA_integer_
        for(z in unique(paste(aa$scenario,aa$shape,sep="|"))) {
          ii <- which(paste(aa$scenario,aa$shape,sep="|")==z)
          eligible_n <- aa$n[ii][aa$qualifies[ii]]
          aa$min_qualifying_n[ii] <- if(length(eligible_n)) min(eligible_n) else NA_integer_
        }
        save_both(aa,"Table15_MC_Power_Based_Candidate_N_Planning")
      }
      plot_v11_pdf_tiff(file.path(figure_dir,"V11_Figure9_MC_Coverage_Diagnostics"),function(){
        dd <- mc_sum[mc_sum$scenario=="null",,drop=FALSE]
        if(!nrow(dd)) {graphics::plot.new();graphics::text(.5,.5,"No null simulations");return(invisible(NULL))}
        graphics::plot(range(dd$n),c(0,1),type="n",xlab="Sample size",
          ylab="Coverage of paired contrast CI")
        graphics::abline(h=1-CFG$alpha,lty=3)
        for(k in seq_along(unique(paste(dd$shape,dd$method)))) {
          grp <- unique(paste(dd$shape,dd$method))[k]
          w <- dd[paste(dd$shape,dd$method)==grp,,drop=FALSE]
          w <- w[order(w$n),,drop=FALSE]
          graphics::lines(w$n,w$coverage,lty=k,lwd=2)
          graphics::points(w$n,w$coverage,pch=k)
        }
        graphics::legend("bottomright",unique(paste(dd$shape,dd$method)),
          lty=seq_along(unique(paste(dd$shape,dd$method))),bty="n",cex=.75)
      })
    }
  }
  # Optional separate validation of age-specific local bootstrap band coverage.
  # Extremely costly: repetitions times B repeated ridge fits; disabled by default.
  if (isTRUE(CFG$run_v11_local_grid_calibration)) {
    local_cal <- tryCatch(run_local_grid_calibration_v11(
      n=CFG$v11_local_mc_n,reps=CFG$v11_local_mc_reps,
      B=CFG$v11_local_mc_bootstrap_B,scenario="null",shape="linear"),
      error=function(e)e)
    if(inherits(local_cal,"error"))
      writeLines(conditionMessage(local_cal),
                 file.path(output,"V11_LOCAL_GRID_CALIBRATION_ERROR.txt"))
    else {
      save_both(local_cal$summary,"Table16_MC_Local_Grid_Band_Coverage")
      write_csv(local_cal$raw,"V11_MC_Local_Grid_Replications.csv")
    }
  }
  logline("V11 add-on FINISHED: new Tables 8-12, Figures 7-8, audit and optional MC")
  invisible(list(regression=regression,endpoint=iut,directions=both,
                 grid_band=mult,compare=compare))
}


main_v10 <- function() {
  # V9.1: request/select workbook FIRST, then establish a stable output folder.
  logline("START — read verified Excel (sheet RStudio_Input)")
  dat <- read_research_data()
  source_excel <- attr(dat, "source_excel")
  if (is.null(source_excel) || !file.exists(source_excel))
    stop("Could not resolve the absolute input Excel path")
  CFG$excel <<- source_excel  # Store exact source for manifests and checksums.
  configured_out <- CFG$out_dir
  is_absolute <- startsWith(configured_out, "/") ||
    startsWith(configured_out, "\\") ||
    (nchar(configured_out) >= 3L &&
       grepl("^[A-Za-z]:", configured_out) &&
       substr(configured_out, 3L, 3L) %in% c("/", "\\"))
  target_dir <- if (is_absolute) configured_out else
    file.path(dirname(source_excel), configured_out)
  if (!dir.exists(target_dir)) {
    created <- dir.create(target_dir, showWarnings = TRUE, recursive = TRUE)
    if (!isTRUE(created) && !dir.exists(target_dir))
      stop("Cannot create results directory: ", target_dir,
           ". Verify that the Excel folder is writable, or set CFG$out_dir ",
           "to a writable absolute path.")
  }
  if (!dir.exists(target_dir))
    stop("Results folder does not exist: ", target_dir)
  # Update the parent configuration: writer/plot functions read CFG globally.
  CFG$out_dir <<- normalizePath(target_dir, winslash = "/", mustWork = TRUE)
  probe <- file.path(CFG$out_dir, paste0(".v9_write_test_", Sys.getpid(), ".tmp"))
  can_write <- tryCatch(file.create(probe),
                        warning = function(w) FALSE,
                        error = function(e) FALSE)
  if (!isTRUE(can_write))
    stop("Output directory is not writable: ", CFG$out_dir,
         ". Move Excel to Documents or choose another writable directory.")
  unlink(probe)
  logline("Verified Excel: ", source_excel)
  logline("Results directory: ", CFG$out_dir)
  qc <- make_qc(dat)
  write_csv(qc$summary,"01_Data_QC.csv")
  write_csv(qc$paired_cells,"02_Paired_D_T1_T2_Cells.csv")
  write_csv(qc$accuracy,"02b_Unadjusted_Diagnostic_Accuracy.csv")
  logline("n=",nrow(dat)," biopsy positive=",sum(dat$D),
          " paired discordance=",sum(dat$T1!=dat$T2))

  # Define truth-free algorithmic checks before analysis.
  stopifnot(isTRUE(run_self_tests_v10())) # algebraic + ridge-bootstrap and logistic separation tests
  logline("Algebra, smooth envelope, and pilot numerical self-tests: PASSED")

  # Freeze the actual fold assignment for this analysis and record it.
  pred <- crossfit_probabilities(dat)
  write_csv(data.frame(X=dat$X,D=dat$D,T1=dat$T1,T2=dat$T2,
                       fold=pred$fold,p_hat=pred$p,q1_hat=pred$q1,q2_hat=pred$q2),
            "03_Out_Of_Fold_Predictions.csv")
  write_csv(pred$diagnostics,"04_Folds_Training_Diagnostics.csv")
  write_csv(pred$calibration,"05_Probability_Prediction_Diagnostics.csv")

  # Primary fixed externally declared loss schedules.
  names_p <- c("equal","fn_075","age_dependent")
  primary <- lapply(names_p,function(nm)
    fixed_policy_inference(dat,pred,relative_fn_loss(dat$X,nm),policy_name=nm))
  prim_tab <- do.call(rbind,lapply(primary,`[[`,"summary"))
  prim_tab$p_holm_exploratory <- stats::p.adjust(prim_tab$p_two_sided,method="holm")
  write_csv(prim_tab,"06_Integrated_Weighted_Kappa_and_Paired_Comparison.csv")
  write_csv(primary[[3]]$influence,"07_EIF_Influences_AgeDependent.csv")
  write_csv(data.frame(X=dat$X,
                       c_age=relative_fn_loss(dat$X,"age_dependent"),
                       s=loss_scale(dat$X),L=box_L(dat$X),U=box_U(dat$X)),
            "08_Prespecified_Loss_Functions.csv")
  logline("Integrated efficient-influence comparisons complete")

  curve <- loss_curves(dat,pred)
  write_csv(curve,"09_Constant_Loss_Kappa_Curve_Pointwise.csv")
  plot_fixed_curve(curve,file.path(CFG$out_dir,"09_Kappa_Curves.png"))
  plot_delta_curve(curve,file.path(CFG$out_dir,"10_Delta_Curve_Pointwise.png"))

  # Smooth conservative dominance: one pre-specified finite tau family.
  env <- prepare_envelope(dat,pred)
  dom <- smooth_multiple_taus(env)
  write_csv(dom$table,"11_Smooth_Conservative_Dominance_Bonferroni.csv")
  logline("Dominance certification computed; fixed tau is conservative at the null")
  if (min(env$minR1_plugin,env$minR2_plugin) <= CFG$ref_denominator_tol)
    warning("Robust dominance sign comparison lacks plugin positivity: do not certify")

  # Secondary local joint model, no forced extrapolation.
  local_ok <- TRUE
  joint <- tryCatch(fit_marginal_spline(dat),error=function(e) {
    local_ok <<- FALSE
    logline("Local marginal logistic model failure: ",conditionMessage(e))
    NULL
  })
  bs_status <- NULL
  bs <- NULL
  local <- NULL
  joint0 <- NULL
  if (local_ok) {
    xg <- as.numeric(stats::quantile(dat$X,c(.05,.95),names=FALSE))
    xg <- seq(xg[1],xg[2],length.out=CFG$local_nx)
    write_csv(joint$audit,"12a_Local_Logistic_Fit_Diagnostics.csv")
    local <- local_surface(joint,xg,CFG$local_c_values)
    write_csv(local,"12_Local_Continuous_Kappa_Surface.csv")
    plot_local_surface(local,file.path(CFG$out_dir,"12_Delta_Age_C_Surface.png"))
    joint0 <- joint_integrated_policy(joint,dat$X,
                        relative_fn_loss(dat$X,"age_dependent"))
    write_csv(joint0,"13_ModelBased_Integrated_Benchmark.csv")
    # joint multinomial is optional and exploratory because 1 record is T1+ T2-.
    if (isTRUE(CFG$joint_multinomial_sensitivity)) {
      zjoint <- tryCatch({
        joint_sens <- fit_joint_spline(dat)
        joint_integrated_policy(joint_sens,dat$X,
                                relative_fn_loss(dat$X,"age_dependent"))
      }, error=function(e)e)
      if (inherits(zjoint,"error")) {
        writeLines(conditionMessage(zjoint),
                   file.path(CFG$out_dir,"13b_JointMultinomial_Failure.txt"))
      } else {
        write_csv(zjoint,"13b_JointMultinomial_Sensitivity.csv")
      }
    }
    # Subject-level NONPARAMETRIC bootstrap, not parametric bootstrapping.
    if (isTRUE(CFG$bootstrap_local)) {
      bs <- tryCatch(bootstrap_joint_local(dat,joint,xg,CFG$local_c_values),
                     error=function(e)e)
      if (inherits(bs,"error")) {
        warning("Bootstrap run had a structural error: ",conditionMessage(bs))
        writeLines(conditionMessage(bs),
                   file.path(CFG$out_dir,"14_Local_Bootstrap_Structural_Error.txt"))
      } else {
        # Always save attempts and reasons, even when band is rejected.
        write_csv(bs$status,"15_Local_Bootstrap_Status.csv")
        write_csv(bs$failure_messages,"16_Local_Bootstrap_Failures.csv")
        write_csv(bs$attempts,"16b_Local_Bootstrap_All_Attempts.csv")
        bs_status <- bs$status
        if (isTRUE(bs$status$accepted_for_finite_grid_CI)) {
          write_csv(bs$results,"14_Local_FiniteGrid_Bootstrap_Bands.csv")
          logline("Exploratory local finite-grid bands created; coverage not certified without MC")
        } else {
          msg <- paste0("Bootstrap success ",bs$status$success,
                        "/",bs$status$requested," below prespecified ",
                        CFG$bootstrap_min_success,
                        "; no confidence bands issued. Inspect attempts CSV.")
          warning(msg)
          writeLines(msg,file.path(CFG$out_dir,"14_Local_Bootstrap_Not_Accepted.txt"))
        }
      }
    }
  }

  # Independent Monte Carlo, disabled by default because publication-grade reps take time.
  # For documented pilot/full MC execution use the separate V11_1_MC_RUNNER.R.
  if (isTRUE(CFG$run_simulation)) {
    for (truth_model in CFG$simulation_truth_models)
      for (sc in c("null","alternative")) for (nn in CFG$simulation_ns) {
      z <- tryCatch(run_monte_carlo(nn,CFG$simulation_reps,sc,
                                    policy="equal",tau=CFG$theta_taus[2],
                                    truth_model=truth_model),
                    error=function(e)e)
      tag <- paste0(truth_model,"_",sc,"_n",nn)
      if (inherits(z,"error")) {
        writeLines(conditionMessage(z),file.path(CFG$out_dir,
                   paste0("MC_FAILED_",tag,".txt")))
      } else {
        write_csv(z$summary,paste0("MC_Summary_",tag,".csv"))
        write_csv(z$raw,paste0("MC_Replications_",tag,".csv"))
        write_csv(z$failed_messages,paste0("MC_Failures_",tag,".csv"))
      }
    }
  }
  if (isTRUE(CFG$run_sample_size_search)) {
    plan <- run_sample_size_search()
    write_csv(plan,"MC_SimulationBased_SampleSize_Search.csv")
  }
  pub_result <- tryCatch(publication_outputs(dat,qc,prim_tab,curve,dom$table,pred,
                               local=local,bs=bs,benchmark=joint0),error=function(e)e)
  if (inherits(pub_result,"error")) {
    msg <- paste("Formatting/export of supplementary publication material failed:",
                 conditionMessage(pub_result))
    warning(msg)
    writeLines(msg,file.path(CFG$out_dir,"PUBLICATION_EXPORT_ERROR.txt"))
  } else {
    logline("Tables: ",file.path(CFG$out_dir,"JBS_Tables"),
            "; figures: ",file.path(CFG$out_dir,"JBS_Figures"))
  }
  # V11's independent inferential and reproducibility package.
  # Any failure is recorded clearly rather than silently altering the primary analysis.
  additional <- tryCatch(v11_extensions(dat,pred,prim_tab,curve,env),
                         error=function(e)e)
  if (inherits(additional,"error")) {
    writeLines(conditionMessage(additional),
               file.path(CFG$out_dir,"V11_EXTENSION_FATAL_ERROR.txt"))
    stop("V11 inferential extension failed: ",conditionMessage(additional))
  }
  write_readme(dat,bs_status)
  utils::capture.output(sessionInfo(),
             file=file.path(CFG$out_dir,"R_SESSION_INFORMATION.txt"))
  logline("FINISHED. Inspect ",CFG$out_dir,"/ and the README before interpretation.")
  invisible(list(primary=prim_tab,dominance=dom$table,diagnostics=qc$summary))
}

# V12: Auto-run moved to explicit staged driver at the end of this file.


# ==============================================================================
# V12 | INTEGRATED RESEARCH WORKFLOW / إدارة مراحل البحث
# ==============================================================================
# Reuses the tested algebraic definitions and numerical routines above but does
# NOT conflate computational convergence with correct statistical coverage.
#
# STAGES (choose BEFORE source()):
#   Sys.setenv(JBS_STAGE="selftest")           -- algebra/numerical checks only
#   Sys.setenv(JBS_STAGE="analysis")           -- real-data analysis + plots/tables
#   Sys.setenv(JBS_STAGE="pilot")              -- analysis + pilot Monte Carlo
#   Sys.setenv(JBS_STAGE="full")               -- analysis + full Monte Carlo
#   Sys.setenv(JBS_STAGE="mc_pilot")           -- simulation only (no Excel)
#   Sys.setenv(JBS_STAGE="mc_full")            -- simulation only (no Excel)
#   Sys.setenv(JBS_STAGE="local_calibration")  -- expensive local-band coverage MC
#
# The MC driver records every cell separately. If the process stops, previously
# completed cells remain on disk; MC output is not silently invented.
# Override run settings with environment variables JBS_MC_REPS, JBS_BOOTSTRAP_B,
# JBS_MULTIPLIER_B, JBS_EXCEL_PATH and JBS_OUTPUT_DIR before sourcing.

v12_env_int <- function(key, default, min_value=1L) {
  v <- Sys.getenv(key, unset="")
  if (!nzchar(v)) return(as.integer(default))
  z <- suppressWarnings(as.integer(v))
  if (is.na(z) || z < min_value) stop(key, " must be an integer >= ", min_value)
  z
}

v12_configure <- function(stage) {
  if (!stage %in% c("selftest", "analysis", "pilot", "full",
                    "mc_pilot", "mc_full", "local_calibration"))
    stop("Unknown JBS_STAGE = ", stage)
  excel <- Sys.getenv("JBS_EXCEL_PATH",unset="")
  if (nzchar(excel)) CFG$excel <<- excel
  out <- Sys.getenv("JBS_OUTPUT_DIR",unset="")
  if (nzchar(out)) CFG$out_dir <<- out
  CFG$run_simulation <<- FALSE                 # legacy simulation off
  CFG$run_sample_size_search <<- FALSE
  CFG$run_v11_simulation <<- FALSE            # unified V12 MC below
  CFG$run_v11_local_grid_calibration <<- FALSE
  CFG$bootstrap_B <<- v12_env_int("JBS_BOOTSTRAP_B",
    if(stage=="full") 999L else 199L, 50L)
  CFG$multiplier_B <<- v12_env_int("JBS_MULTIPLIER_B",1999L,199L)
  invisible(stage)
}

v12_prepare_mc_dir <- function() {
  if (!dir.exists(CFG$out_dir)) {
    if (!dir.create(CFG$out_dir,recursive=TRUE,showWarnings=TRUE) &&
        !dir.exists(CFG$out_dir)) stop("Could not create output folder")
  }
  CFG$out_dir <<- normalizePath(CFG$out_dir,winslash="/",mustWork=TRUE)
  probe <- file.path(CFG$out_dir,".v12_write_test")
  if (!isTRUE(file.create(probe))) stop("MC output folder is not writable: ",CFG$out_dir)
  unlink(probe)
}

# Extra exact-algebra self-tests independent of the regression fit.
# This tests the *definition* of chance loss, signed-kappa numerator, envelope,
# as well as the error of smooth bounds. It does not prove statistical novelty.
v12_algebra_audit <- function() {
  x <- seq(15,75,length.out=151)
  s <- rep(1,length(x))
  p <- stats::plogis(-1.8+.04*(x-30))
  se1 <- stats::plogis(.9+.018*(x-30)); se2 <- stats::plogis(1.2+.012*(x-30))
  sp1 <- stats::plogis(2.9-.005*(x-30)); sp2 <- stats::plogis(2.0-.005*(x-30))
  q1 <- p*se1+(1-p)*(1-sp1); q2 <- p*se2+(1-p)*(1-sp2)
  cc <- .4+.45*stats::plogis((x-35)/12)
  a1 <- s*(cc*p*(1-se1)+(1-cc)*(1-p)*(1-sp1))
  a2 <- s*(cc*p*(1-se2)+(1-cc)*(1-p)*(1-sp2))
  r1 <- s*(cc*p*(1-q1)+(1-cc)*(1-p)*q1)
  r2 <- s*(cc*p*(1-q2)+(1-cc)*(1-p)*q2)
  cov1 <- s*p*(1-p)*(se1+sp1-1)
  cov2 <- s*p*(1-p)*(se2+sp2-1)
  err_cov <- max(abs(r1-a1-cov1),abs(r2-a2-cov2))
  A1<-mean(a1);A2<-mean(a2);C1<-mean(cov1);C2<-mean(cov2)
  R1<-mean(r1);R2<-mean(r2)
  b_num <- C1*R2-C2*R1
  delta <- C1/R1-C2/R2
  err_sign <- abs(delta-b_num/(R1*R2))
  L<-rep(.4,length(x));U<-rep(.9,length(x))
  V1<-mean(s*(1-p)*q1); V2<-mean(s*(1-p)*q2)
  h<-s*(C1*(p-q2)-C2*(p-q1))
  b0<-C1*V2-C2*V1
  blo<-b0+mean(pmin(L*h,U*h)); bup<-b0+mean(pmax(L*h,U*h))
  # Compute B(c) in two different ways at a random loss function.
  Bref <- C1*R2-C2*R1
  Bexplicit <- b0+mean(cc*h)
  err_B <- abs(Bref-Bexplicit)
  err_box <- max(blo-Bexplicit,Bexplicit-bup,0)
  taus<-c(.0001,.0005,.002)
  err_smooth<-0
  for(tau in taus) {
    low<-soft_min(h,L,U,tau)
    up<-soft_max_v11(h,L,U,tau)
    err_smooth<-max(err_smooth,
       max(low-pmin(L*h,U*h)), max(pmax(L*h,U*h)-up),
       max(pmin(L*h,U*h)-low-tau*log(2)),
       max(up-pmax(L*h,U*h)-tau*log(2)))
  }
  audit<-data.frame(
    identity=c("R-A=weighted covariance","Delta=B/(R1*R2)",
       "B=b0+E[c(X)h(X)]","box extrema bound random c", "smooth lower-upper errors"),
    maximum_absolute_violation=c(err_cov,err_sign,err_B,err_box,err_smooth),
    permitted_tolerance=rep(1e-9,5L))
  audit$passed <- audit$maximum_absolute_violation <= audit$permitted_tolerance
  if (!all(audit$passed)) stop("V12 identity test failed")
  audit
}

# Non-numerical manuscript logic. These categories are checked by humans,
# never converted into automatic claims of publication suitability.
v12_claim_register <- function() {
  data.frame(
    claim=c("Loss-weighted conditional covariance identity",
      "Finite-box functional loss extrema", "Cross-fitted IF / Delta asymptotics",
      "Independent stacked-sandwich regression benchmark",
      "Finite-grid multiplier band", "Local ridge bootstrap band",
      "Fixed-tau functional dominance test", "Uniform validity near ties",
      "Age-dependent clinical loss elicitation", "Biopsy verification design",
      "Monte Carlo performance", "Novelty & literature assessment"),
    current_status=c("ALGEBRAIC", "ALGEBRAIC", "THEOREM_WITH_ASSUMPTIONS",
      "MODEL_BASED", "ASYMPTOTIC_FINITE_GRID", "EXPLORATORY",
      "POINTWISE_CONSERVATIVE", "NOT_ESTABLISHED",
      "HYPOTHETICAL", "NOT_VERIFIED", "MUST_RUN_REAL_R_MC",
      "REQUIRES_INDEPENDENT_REVIEW"),
    author_action=c("Retain positivity conditions", "State measurable pointwise box constraints",
      "Check nuisance product rates, bounded s, R positive, variance nonzero",
      "Evaluate misspecification and finite-sample separation",
      "Validate coverage in MC; predefine c-grid", 
      "Report rejected resamples and MC conditional/unconditional coverage",
      "Do not claim exact size or uniform validity",
      "Derive new method or avoid claim", "Elicit losses with domain experts",
      "Audit clinical dataset source, selection and verification bias",
      "Inspect Tables 13-16 and failure diagnostics", 
      "Systematically compare with prior literature"),stringsAsFactors=FALSE)
}

v12_mc_plan <- function(mode="pilot") {
  if(mode=="pilot") {
    list(reps=v12_env_int("JBS_MC_REPS",50L,20L),
      ns=c(250L,858L),
      scenarios=c("null","t2_local","t2_strong"),
      shapes=c("linear","nonlinear"), stage="PILOT_NOT_PUBLICATION")
  } else {
    list(reps=v12_env_int("JBS_MC_REPS",1000L,20L),
      ns=c(250L,500L,858L,1500L),
      scenarios=c("null","t2_local","t1_strong","t2_strong","tradeoff"),
      shapes=c("linear","nonlinear"),
      stage="FULL_DESIGN_REQUIRES_REVIEW")
  }
}

v12_exact_lower_binomial <- function(events, n, alpha=.05) {
  if(n<1L || is.na(events)) return(NA_real_)
  if(events<=0L) return(0)
  stats::qbeta(alpha,events,n-events+1L)
}

v12_add_mc_fields <- function(mc) {
  raw<-mc$raw; truth<-mc$truth$delta
  n<-nrow(raw)
  for(met in c("cf","reg")) {
    good <- if(met=="cf") raw$primary_ok else raw$regression_ok
    lo <- if(met=="cf") raw$cf_ci_low else raw$reg_low
    hi <- if(met=="cf") raw$cf_ci_high else raw$reg_high
    rej <- if(met=="cf") raw$cf_reject else raw$reg_reject
    valid <- as.logical(good) & is.finite(lo) & is.finite(hi)
    coverage_all <- sum(valid & lo<=truth & hi>=truth,na.rm=TRUE)/n
    reject_all <- sum(valid & rej==1,na.rm=TRUE)/n
    rows <- which(mc$summary$method==if(met=="cf") "crossfit" else "regression_stacked_delta")
    if(length(rows)) {
      mc$summary$coverage_and_success[rows] <- coverage_all
      mc$summary$rejection_and_success[rows] <- reject_all
      mc$summary$nominal_ci_valid_fraction[rows] <- mean(valid)
      mc$summary$coverage_monte_carlo_se_unconditional[rows] <- sqrt(coverage_all*(1-coverage_all)/n)
      mc$summary$rejection_monte_carlo_se_unconditional[rows] <- sqrt(reject_all*(1-reject_all)/n)
    }
  }
  mc
}

v12_generate_mc_artifacts <- function(base, rows, doms, mode) {
  if (!length(rows)) stop("No successful Monte Carlo cells")
  tb13<-do.call(rbind,rows); tb14<-do.call(rbind,doms)
  utils::write.csv(tb13,file.path(base,"Table13_MC_Bias_RMSE_Coverage_Rejection.csv"),row.names=FALSE)
  utils::write.csv(tb14,file.path(base,"Table14_MC_Robust_Dominance.csv"),row.names=FALSE)
  # Type-I error is measurable ONLY in true-null scenarios; otherwise rejection=power.
  tb13$empirical_typeI <- ifelse(tb13$delta_null_true,tb13$rejection_rate,NA_real_)
  tb13$empirical_power <- ifelse(!tb13$delta_null_true,tb13$rejection_rate,NA_real_)
  tb13$rejection_failure_adjusted <- tb13$rejection_and_success
  utils::write.csv(tb13,file.path(base,"Table13b_MC_TypeI_Power_Failure_Adjusted.csv"),row.names=FALSE)
  # Candidate size: success-aware power with a one-sided EXACT binomial lower bound.
  choices <- tb13[tb13$method=="crossfit" & !tb13$delta_null_true,,drop=FALSE]
  if(nrow(choices)) {
    choices$success_and_reject_n <- round(choices$rejection_and_success*choices$total_reps)
    choices$lower95_power <- mapply(v12_exact_lower_binomial,
      choices$success_and_reject_n,choices$total_reps)
    choices$meets_target <- choices$lower95_power>=.80 & choices$primary_failure_fraction<=.05
    choices$chosen_candidate_n <- NA_integer_
    group <- paste(choices$scenario,choices$shape,choices$policy,
      if("dependence" %in% names(choices)) choices$dependence else 0.55,sep="::")
    for(k in unique(group)) {
      ix<-which(group==k)
      candidates<-choices$n[ix][choices$meets_target[ix]]
      if(length(candidates)) choices$chosen_candidate_n[ix]<-min(candidates)
    }
    utils::write.csv(choices,file.path(base,"Table15_MC_Sample_Size_Candidates.csv"),row.names=FALSE)
  }
  # Figures per manuscript methodological comparison, without altering estimates.
  for(kind in c("coverage","rejection")) {
    plotter <- function() {
      yy <- if(kind=="coverage") "coverage" else "rejection_rate"
      x <- tb13$n; y <- tb13[[yy]]
      good <- is.finite(x)&is.finite(y)
      if(!any(good)) {graphics::plot.new();graphics::text(.5,.5,"No complete Monte Carlo results");return()}
      graphics::plot(range(x[good]),c(0,1),type="n",xlab="Sample size",
         ylab=if(kind=="coverage") "95% CI coverage (conditional on successful fits)" else
           "Rejection probability",main=paste("Monte Carlo",kind,"-",mode))
      if(kind=="coverage") graphics::abline(h=.95,lty=3)
      if(kind=="rejection") graphics::abline(h=.05,lty=3)
      nullrows <- tb13[tb13$scenario=="null",,drop=FALSE]
      shapes <- unique(nullrows$shape)
      id<-0L
      for(shape in shapes) for(method in unique(nullrows$method)) {
        id<-id+1L
        sub<-nullrows[nullrows$shape==shape&nullrows$method==method,,drop=FALSE]
        sub<-sub[order(sub$n),,drop=FALSE]
        if(nrow(sub)) {graphics::lines(sub$n,sub[[yy]],lty=id,lwd=2);
                      graphics::points(sub$n,sub[[yy]],pch=id)}
      }
      graphics::legend("topright",legend=as.vector(outer(shapes,unique(nullrows$method),paste)),
                       lty=seq_len(id),bty="n",cex=.65)
    }
    plot_v11_pdf_tiff(file.path(base,paste0("Figure_MC_",kind)),plotter)
  }
  invisible(tb13)
}

v12_run_mc <- function(mode="pilot") {
  plan<-v12_mc_plan(mode)
  root<-file.path(CFG$out_dir,paste0("MONTE_CARLO_",toupper(mode)))
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(root))stop("Cannot make Monte Carlo folder: ",root)
  spec<-c(paste("Mode:",plan$stage),paste("Reps per cell:",plan$reps),
      paste("N:",paste(plan$ns,collapse=", ")),
      paste("Scenarios:",paste(plan$scenarios,collapse=", ")),
      paste("Shapes:",paste(plan$shapes,collapse=", ")),
      "MC summary coverage and rejection conditional on estimability",
      "Failure-adjusted rates count failed fits as no interval and no rejection",
      "All rows are simulated; real observed UCI dataset is NOT the simulation truth",
      "No sample size is certified without external clinical effect assumptions",
      "Pilot simulation cannot validate Type-I rate or confidence coverage",
      "Fixed-tau dominance conservative; regularity and near-tie validity still require study")
  writeLines(spec,file.path(root,"MC_DESIGN_AND_CAUTION.txt"))
  rows<-list();doms<-list(); errors<-list();idx<-0L
  for(nn in plan$ns) for(sc in plan$scenarios) for(shape in plan$shapes){
    idx<-idx+1L
    tag<-paste0("n",nn,"_",sc,"_",shape)
    logline("MC START ",tag," reps=",plan$reps)
    ans<-tryCatch(run_mc_v11(n=nn,reps=plan$reps,scenario=sc,shape=shape,
          policy="age_dependent",tau=.0005,seed=20261301L+idx*41003L),
          error=function(e)e)
    if(inherits(ans,"error")) {
      msg<-conditionMessage(ans)
      errors[[tag]]<-data.frame(scenario_cell=tag,error=msg)
      writeLines(msg,file.path(root,paste0("ERROR_",tag,".txt")))
      logline("MC ERROR ",tag,": ",msg)
    } else {
      ans<-v12_add_mc_fields(ans)
      rows[[tag]]<-ans$summary; doms[[tag]]<-ans$dominance
      utils::write.csv(ans$raw,file.path(root,paste0("RAW_",tag,".csv")),row.names=FALSE)
      # CHECKPOINT tables: each successfully completed cell is permanent on disk.
      utils::write.csv(do.call(rbind,rows),file.path(root,"CHECKPOINT_MC_SUMMARY.csv"),row.names=FALSE)
      utils::write.csv(do.call(rbind,doms),file.path(root,"CHECKPOINT_MC_DOMINANCE.csv"),row.names=FALSE)
      logline("MC DONE ",tag)
    }
  }
  if(length(errors))utils::write.csv(do.call(rbind,errors),file.path(root,"MC_SCENARIO_ERRORS.csv"),row.names=FALSE)
  if(length(rows))v12_generate_mc_artifacts(root,rows,doms,mode)
  utils::capture.output(sessionInfo(),file=file.path(root,"R_SESSION_INFO.txt"))
  missing <- length(plan$ns)*length(plan$scenarios)*length(plan$shapes)-length(rows)
  writeLines(c(paste("Completed scenario cells:",length(rows)),
               paste("Failed scenario cells:",missing),
               paste("Run stage:",plan$stage),
               "Do not treat conditional-on-success performance as unconditional reliability."),
             file.path(root,"MC_COMPLETENESS_AUDIT.txt"))
  invisible(root)
}

v12_write_research_status <- function() {
  if(!dir.exists(CFG$out_dir))return(invisible(NULL))
  utils::write.csv(v12_claim_register(),
      file.path(CFG$out_dir,"Table17_Claim_and_Validation_Register.csv"),row.names=FALSE)
  utils::write.csv(v12_algebra_audit(),
      file.path(CFG$out_dir,"V12_Algebra_Independent_SelfTests.csv"),row.names=FALSE)
  utils::capture.output(str(CFG),file=file.path(CFG$out_dir,"V12_Analysis_Configuration.txt"))
  utils::capture.output(sessionInfo(),file=file.path(CFG$out_dir,"V12_R_Session_Info.txt"))
  if(!is.null(CFG$excel) && file.exists(CFG$excel))
    utils::write.csv(data.frame(file=normalizePath(CFG$excel,winslash="/"),
       md5=unname(tools::md5sum(CFG$excel))),
       file.path(CFG$out_dir,"V12_Input_File_Checksum.csv"),row.names=FALSE)
  ctn<-list.files(CFG$out_dir,recursive=TRUE,full.names=FALSE)
  utils::write.csv(data.frame(relative_path=ctn),
       file.path(CFG$out_dir,"V12_Results_File_Manifest.csv"),row.names=FALSE)
  invisible(TRUE)
}

v12_run <- function(stage=tolower(Sys.getenv("JBS_STAGE",unset="analysis"))) {
  v12_configure(stage)
  logline("V12 PIPELINE stage=",stage)
  # Mathematical and implementation checks are always run before data analysis.
  print(v12_algebra_audit(),row.names=FALSE)
  if (!isTRUE(run_self_tests_v11())) stop("V11 engine self-tests failed")
  if (stage=="selftest") {
    logline("SELFTEST FINISHED: no patient data analysis or MC run")
    return(invisible(TRUE))
  }
  if (stage %in% c("analysis","pilot","full")) {
    # Existing engine creates JBS_Tables & JBS_Figures, Delta benchmark,
    # paired EIF, constant-loss IUT, both directions of loss-functional dominance,
    # bootstrap finite-grid bands, QC, and reproducibility files.
    result <- main_v10()
    v12_write_research_status()
  } else {
    v12_prepare_mc_dir()
  }
  if(stage %in% c("pilot","full","mc_pilot","mc_full")) {
    mode <- if(stage %in% c("pilot","mc_pilot")) "pilot" else "full"
    v12_run_mc(mode)
  }
  if(stage=="local_calibration") {
    # Computationally expensive. Smaller pilot configuration by default.
    nrep<-v12_env_int("JBS_LOCAL_REPS",60L,20L)
    B<-v12_env_int("JBS_LOCAL_B",199L,50L)
    local<-run_local_grid_calibration_v11(n=858L,reps=nrep,B=B,
                            scenario="null",shape="nonlinear")
    utils::write.csv(local$summary,file.path(CFG$out_dir,"Table16_Local_Band_MC_Coverage.csv"),row.names=FALSE)
    utils::write.csv(local$raw,file.path(CFG$out_dir,"RAW_Local_Band_MC_Replications.csv"),row.names=FALSE)
  }
  v12_write_research_status()
  logline("V12 COMPLETE. Results: ",CFG$out_dir)
  invisible(TRUE)
}


# ============================================================================
# V13.2 INTEGRATED RESEARCH: PRIMARY KAPPA VS EXTERNAL PV MANUSCRIPT
# ============================================================================
# Distinct statistical questions - never pool their p values or effect measures.
# PRIMARY: integrated covariate-dependent loss-weighted kappa for two PAIRED BDTs,
#          X continuous; orthogonal CF estimator, regression-sandwich Delta,
#          functional dominance, finite-grid multiplier/local-bootstrap inference.
# CONTEXT: paired PPV & NPV descriptive/global Wald, X grouped only for a
#          deliberately separate REFERENCE analysis; DOES NOT claim to reproduce
#          the published regression/WLS estimators or global tests exactly.
# Published comparator paper: Roldan-Nofuentes and Regad, DOI
#          10.1080/10543406.2025.2547589 (predictive values, categorical X).
# All loss policies in this implementation are HYPOTHETICAL sensitivity analyses.
# No data, simulation, code or decision rule guarantees journal acceptance.

v13_log <- function(...) cat(format(Sys.time(), '%H:%M:%S'), ' | V13 | ',
                              paste0(...,collapse=''), '\n', sep='')
v13_env_int <- function(key, def, lo=1L) {
  val <- Sys.getenv(key, unset='')
  if (!nzchar(val)) return(as.integer(def))
  x <- suppressWarnings(as.integer(val))
  if (!is.finite(x) || x < lo) stop(key,' must be an integer >= ',lo)
  x
}
v13_safe_output <- function(path) {
  ok <- dir.exists(path) || dir.create(path,recursive=TRUE,showWarnings=FALSE)
  if (!isTRUE(ok)) stop('Cannot create directory: ',path)
  p <- normalizePath(path,mustWork=TRUE,winslash='/')
  f <- tempfile('.test_',tmpdir=p)
  if (!isTRUE(file.create(f))) stop('Output path not writable: ',p)
  unlink(f);p
}
v13_write <- function(z, filename, dir=CFG$out_dir) {
  utils::write.csv(z,file.path(v13_safe_output(dir),filename),row.names=FALSE,
                   na='',fileEncoding='UTF-8')
}
v13_status <- function(mode, pass, reason, evidence='',type='CHECK') {
  data.frame(mode=mode,status=if(isTRUE(pass)) 'PASS' else 'BLOCKED',
             reason=reason,evidence=evidence,type=type,stringsAsFactors=FALSE)
}

# -------- Context ONLY: paired PPV and NPV, NOT loss-weighted kappa -----------
v13_pv_global <- function(dat, alpha=CFG$alpha) {
  # P(D=1|Tj=1) and P(D=0|Tj=0), both tests measured in SAME subjects.
  # Ratio EIF accounts for paired covariance and is asymptotic Wald (2 df).
  stopifnot(all(c('D','T1','T2') %in% names(dat)))
  d <- as.numeric(dat$D); t1 <- as.numeric(dat$T1); t2 <- as.numeric(dat$T2)
  n <- length(d)
  pv <- function(t) {
    prp <- mean(t); prn <- mean(1-t)
    if (min(prp,prn) < .005) stop('Predictive-value denominator near zero')
    ppv <- mean(d*t)/prp
    npv <- mean((1-d)*(1-t))/prn
    list(ppv=ppv,npv=npv,
         ifppv=t*(d-ppv)/prp,
         ifnpv=(1-t)*((1-d)-npv)/prn)
  }
  a<-pv(t1);b<-pv(t2)
  v <- c(ppv_difference=a$ppv-b$ppv,npv_difference=a$npv-b$npv)
  phi <- cbind(a$ifppv-b$ifppv,a$ifnpv-b$ifnpv)
  vc <- stats::cov(phi)/n
  eig <- eigen(vc,symmetric=TRUE,only.values=TRUE)$values
  if (!all(is.finite(eig)) || min(eig)<=1e-12 || max(eig)/min(eig)>1e12)
    stop('PV contrast covariance singular; global Wald test not estimable')
  stat <- drop(t(v) %*% solve(vc,v))
  data.frame(n=n,PPV_T1=a$ppv,PPV_T2=b$ppv,NPV_T1=a$npv,NPV_T2=b$npv,
    delta_PPV=unname(v[1]),delta_NPV=unname(v[2]),
    se_delta_PPV=sqrt(vc[1,1]),se_delta_NPV=sqrt(vc[2,2]),
    Q_global=stat,df=2L,p_global=stats::pchisq(stat,2,lower.tail=FALSE),
    reject_global=as.integer(stats::pchisq(stat,2,lower.tail=FALSE)<alpha),
    estimand='Paired marginal PPV/NPV global equality (NOT weighted kappa)',
    stringsAsFactors=FALSE)
}

v13_pv_categorical_context <- function(dat, category_count=2L, alpha=CFG$alpha) {
  # Contextual stratified DISPLAY only; binning of X is not primary modeling.
  # Prespecified cut for 2 categories, tertile quantiles for 3 categories.
  stopifnot(category_count %in% c(2L,3L))
  z <- dat
  if (category_count==2L) {
    # clinically illustrative only; not a discovered clinical threshold
    z$category <- ifelse(z$X<30L,'age_below_30','age_30_plus')
    design <- 'age < 30 versus >= 30; illustrative fixed threshold'
  } else {
    knots <- unique(as.numeric(stats::quantile(z$X,c(0,1/3,2/3,1),na.rm=TRUE)))
    if(length(knots)<4L) stop('Age tertile cutpoints not distinct')
    knots[1]<- -Inf; knots[4]<-Inf
    z$category<-as.character(cut(z$X,breaks=knots,include.lowest=TRUE))
    design <- 'empirical age tertiles; sample-derived exploratory groups'
  }
  global <- tryCatch(v13_pv_global(z,alpha),error=function(e)e)
  if(inherits(global,'error')) stop(conditionMessage(global))
  sub <- lapply(split(z,z$category),function(v){
    ans <- tryCatch(v13_pv_global(v,alpha),error=function(e)e)
    if(inherits(ans,'error')) data.frame(n=nrow(v),p_global=NA_real_,
                     delta_PPV=NA_real_,delta_NPV=NA_real_,issue=conditionMessage(ans))
    else {ans$issue<-''; ans}
  })
  # Keep a harmonized key-field view: no combining marginal and conditional tests.
  strat <- do.call(rbind,lapply(names(sub),function(k) {
    q<-sub[[k]]
    data.frame(group=k,n=q$n[1],delta_PPV=q$delta_PPV[1],
       delta_NPV=q$delta_NPV[1],p_global=q$p_global[1],
       issue=q$issue[1],stratification=design,stringsAsFactors=FALSE)
  }))
  list(marginal=global,by_stratum=strat,design=design)
}

# Dedicated categorical-data generator. Distinct from continuous-X kappa DGM.
v13_generate_pv <- function(n=500L, categories=2L, scenario='null', dependence=.60,
                           seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  if(!categories %in% c(2L,3L))stop('categories must be 2 or 3')
  if(!scenario %in% c('null','alternative'))stop('Unknown PV scenario')
  probs <- if(categories==2L)c(.6,.4)else c(.45,.35,.2)
  K<-sample.int(categories,n,replace=TRUE,prob=probs)
  p <- c(.12,.32,.44)[K]
  D <- stats::rbinom(n,1,p)
  se1<-c(.72,.78,.83)[K]; sp1<-c(.93,.90,.87)[K]
  if(scenario=='null') {se2<-se1;sp2<-sp1}
  else {se2<-pmin(.99,se1+.13);sp2<-pmax(.5,sp1-.045)}
  r1<-ifelse(D==1,se1,1-sp1);r2<-ifelse(D==1,se2,1-sp2)
  u<-stats::runif(n);v<-ifelse(stats::runif(n)<dependence,u,stats::runif(n))
  data.frame(X_category=K,D=D,T1=as.integer(u<r1),T2=as.integer(v<r2))
}

# Deterministic truth for categorical PV DGM; integrated across strata.
v13_truth_pv <- function(categories=2L,scenario='null') {
  weights<-if(categories==2L)c(.6,.4)else c(.45,.35,.2)
  p<-c(.12,.32,.44)[seq_len(categories)]
  s1<-c(.72,.78,.83)[seq_len(categories)]
  t1<-c(.93,.90,.87)[seq_len(categories)]
  s2<-if(scenario=='null')s1 else pmin(.99,s1+.13)
  t2<-if(scenario=='null')t1 else pmax(.5,t1-.045)
  comp<-function(se,sp) {
    pos<-sum(weights*(p*se+(1-p)*(1-sp)))
    neg<-sum(weights*(p*(1-se)+(1-p)*sp))
    c(ppv=sum(weights*p*se)/pos,npv=sum(weights*(1-p)*sp)/neg)
  }
  a<-comp(s1,t1);b<-comp(s2,t2)
  c(ppv_difference=a[1]-b[1],npv_difference=a[2]-b[2])
}
v13_mc_pv <- function(mode='pilot', reps=50L, outdir=CFG$out_dir) {
  outdir<-v13_safe_output(file.path(outdir,paste0('REFERENCE_PV_MC_',toupper(mode))))
  ns<-if(mode=='pilot')c(250L,858L)else c(100L,250L,500L,858L)
  groups<-c(2L,3L);scenarios<-c('null','alternative')
  results<-list();caseid<-0L
  for(K in groups)for(sc in scenarios)for(nn in ns){
    caseid<-caseid+1L
    tag<-paste('K',K,'scenario',sc,'n',nn,sep='_')
    target<-file.path(outdir,paste0('RAW_',tag,'.csv'))
    if(file.exists(target)){
      tab<-utils::read.csv(target,stringsAsFactors=FALSE)
      if(nrow(tab)!=reps)stop('Existing PV MC has a different replication count: ',tag)
    } else {
      rows<-vector('list',reps)
      for(i in seq_len(reps)){
        dd<-v13_generate_pv(nn,K,sc,seed=20261017L+caseid*100000L+i)
        pv<-tryCatch(v13_pv_global(dd),error=function(e)e)
        if(inherits(pv,'error')) rows[[i]]<-data.frame(rep=i,ok=FALSE,
          p=NA_real_,reject=NA_integer_,message=conditionMessage(pv))
        else rows[[i]]<-data.frame(rep=i,ok=TRUE,p=pv$p_global,
          reject=pv$reject_global,message='')
      }
      tab<-do.call(rbind,rows);v13_write(tab,basename(target),outdir)
    }
    good<-tab$ok & !is.na(tab$reject)
    q<-if(sum(good))mean(tab$reject[good]) else NA_real_
    q_all<-sum(tab$reject[good]==1L)/reps
    tr<-v13_truth_pv(K,sc)
    results[[tag]]<-data.frame(groups=K,scenario=sc,n=nn,reps=reps,
      completed=sum(good),failed=reps-sum(good),
      estimated_typeI=if(sc=='null')q else NA_real_,
      estimated_power=if(sc!='null')q else NA_real_,
      coverage_note='PV global test only; NOT loss-weighted kappa',
      rejection_and_success=q_all,
      mcse_conditional=if(is.finite(q))sqrt(q*(1-q)/sum(good))else NA_real_,
      truth_delta_ppv=unname(tr[1]),truth_delta_npv=unname(tr[2]),
      stringsAsFactors=FALSE)
  }
  out<-do.call(rbind,results)
  v13_write(out,'TablePV_Context_MC_TypeI_Power_Failures.csv',outdir)
  out
}

# Simulated individual-level datasets: expressly SYNTHETIC, never used as UCI.
v13_generate_demo <- function(n=858L, scenario='null', shape='nonlinear',
                              dir=CFG$out_dir,seed=20261017L) {
  set.seed(seed)
  dat<-mc_generate_v11(n=n,scenario=scenario,shape=shape,dependence=.55)
  truth<-mc_truth_v11(scenario=scenario,shape=shape,policy='age_dependent')
  outdir<-v13_safe_output(file.path(dir,'SYNTHETIC_DATA'))
  name<-paste0('SYNTHETIC_',scenario,'_',shape,'_n',n,'.csv')
  v13_write(dat,name,outdir)
  v13_write(data.frame(scenario=scenario,shape=shape,n=n,
    truth_delta=truth$delta,truth_lower_B=truth$lower,truth_upper_B=truth$upper,
    NOT_REAL_CLINICAL_DATA=TRUE),paste0('TRUTH_',name),outdir)
  v13_log('Synthetic data and independent true estimands saved: ',outdir)
  invisible(dat)
}

# NO-X comparator: same clinical data, DIFFERENT target (unconditional chance).
# It is NOT a valid same-estimand comparator for the conditional chance-adjusted
# Kappa in the main paper, and it is NOT the categorical reference manuscript WLS.
v13_no_x_kappa_context <- function(dat, weights=c(0.5,0.75), alpha=CFG$alpha){
  d<-dat$D; t1<-dat$T1; t2<-dat$T2; n<-length(d)
  p<-mean(d); q1<-mean(t1); q2<-mean(t2)
  one<-function(t,cv){
    q<-mean(t)
    a<-cv*d*(1-t)+(1-cv)*(1-d)*t
    A<-mean(a)
    u<-(d-p)*(t-q)
    C<-mean(u); R<-A+C
    if(!is.finite(R)||R<=CFG$ref_denominator_tol)
      stop('No-X reference chance loss is nonpositive')
    phiA<-a-A;phiC<-u-C
    list(kappa=C/R,phi=(A*phiC-C*phiA)/R^2)
  }
  do.call(rbind,lapply(weights,function(cv){
    a<-one(t1,cv);b<-one(t2,cv)
    diff<-a$kappa-b$kappa
    se<-stats::sd(a$phi-b$phi)/sqrt(n)
    z<-stats::qnorm(1-alpha/2)
    data.frame(constant_loss_weight=cv,kappa1_noX=a$kappa,
      kappa2_noX=b$kappa,paired_delta_noX=diff,se_delta_noX=se,
      lower=diff-z*se,upper=diff+z*se,
      target='MARGINAL chance reference; distinct from conditional chance-adjusted Kappa',
      stringsAsFactors=FALSE)
  }))
}

# Paper-to-research methodological correspondence: REF paper is a DIFFERENT target.
v13_methodology_map <- function(){
  data.frame(
    phase=c("Data structure","Target parameter","Null hypotheses","Paired dependence",
      "Covariate modelling","Estimation comparator","Monte Carlo generation",
      "Null size and power","Estimator operating characteristics","Sample size",
      "Software reproducibility","Applied examples","Scientific interpretation"),
    reference_predictive_values=c(
      "D,T1,T2 and categorical covariate patterns",
      "Two PPVs + two NPVs, marginal/global",
      "Joint equality of both PPV contrasts and NPV contrasts; 2-df global Wald",
      "Paired joint multinomial tables within categories",
      "Categorical regression and weighted least squares",
      "Categorical RM vs WLS vs covariate-ignored Wu method",
      "2/3-category multinomial simulations with dependence settings",
      "Simulation Type-I and power at prespecified alpha",
      "Relative bias and numerical estimation difficulties",
      "Study n=50,100,...,500; not a direct loss-kappa formula",
      "Published R functions for reference PPV problem",
      "Two clinical applications of predictive-value comparison",
      "PVs are different parameters from weighted Kappa"),
    present_loss_weighted_kappa=c(
      "Individual paired D,T1,T2 and numeric-scale X; same-subject dependence",
      "Kappa1, Kappa2, Delta under predeclared c(X) and s(X)",
      "Delta=0; directional constant-c IUT; functional-envelope sign",
      "Shared-subject IF/Sandwich; simulated correlated paired Bernoulli",
      "Cross-fitted logistic splines, model-based regression Delta; X not discretized in primary",
      "CF/EIF vs independent regression Delta; direct clinical loss as distinct outcome",
      "Continuous-X null/near-null/strong/tradeoff, line/nonlinear; optional dependence sensitivity",
      "Monte Carlo nominal size, power and MCSE for each n and scenario",
      "Bias RMSE empirical SD mean SE coverage failure-adjusted coverage and model failures",
      "Power-based candidate n with binomial lower-bound safeguard",
      "Staged RStudio runs, fixed seed, CSV checkpoints, provenance and diagnostics",
      "UCI cervical study + clearly labeled SYNTHETIC data",
      "Do not claim equality with PPV paper estimand, WLS method or its novelty"),
    stringsAsFactors=FALSE)
}

# Sensitivity to the nuisance specification: estimator unchanged, nuisance altered.
# DO NOT interpret as data-driven selection of the best significant p-value.
v13_nuisance_sensitivity <- function(dat, outdir=CFG$out_dir){
  settings<-data.frame(model=c('linear','spline','spline','spline'),
                       df=c(1L,2L,3L,4L),stringsAsFactors=FALSE)
  rows<-vector('list',nrow(settings))
  for(k in seq_len(nrow(settings))){
    method<-settings$model[k]; df<-settings$df[k]
    res<-tryCatch({
      pr<-crossfit_probabilities(dat,K=CFG$folds,method=method,
                                  df=df,seed=CFG$seed)
      est<-fixed_policy_inference(dat,pr,relative_fn_loss(dat$X,'age_dependent'),
                                  policy_name='age_dependent')$summary
      data.frame(model=method,df=df,success=TRUE,
         delta=est$delta,se_delta=est$se_delta,
         ci_low=est$ci_delta_low,ci_high=est$ci_delta_high,
         p_two_sided=est$p_two_sided,
         Brier_D=pr$calibration$Brier[1],
         Brier_T1=pr$calibration$Brier[2],
         Brier_T2=pr$calibration$Brier[3],error='',stringsAsFactors=FALSE)
    },error=function(e) data.frame(model=method,df=df,success=FALSE,
         delta=NA_real_,se_delta=NA_real_,ci_low=NA_real_,ci_high=NA_real_,
         p_two_sided=NA_real_,Brier_D=NA_real_,Brier_T1=NA_real_,
         Brier_T2=NA_real_,error=conditionMessage(e),stringsAsFactors=FALSE))
    rows[[k]]<-res
  }
  ans<-do.call(rbind,rows)
  v13_write(ans,'V13_Nuisance_Model_Sensitivity.csv',outdir)
  ans
}

# Dependence settings are optional, distinct Monte Carlo scenario cells.
v13_dependence_levels <- function(){
  raw<-Sys.getenv('JBS_DEPENDENCE_LEVELS',unset='0.55')
  zz<-suppressWarnings(as.numeric(trimws(strsplit(raw,',',fixed=TRUE)[[1]])))
  if(!length(zz)||any(!is.finite(zz))||any(zz<0|zz>1)||anyDuplicated(zz)||
     anyDuplicated(sprintf('%.6f',zz)))
     stop('JBS_DEPENDENCE_LEVELS: distinct comma-separated values in [0,1] required')
  zz
}

# Rigorous MC for PRIMARY kappa: restartable, independent true DGM, failures kept.
v13_mc_reps <- function(mode){
  stage<-tolower(Sys.getenv('JBS_STAGE',unset=''))
  specific<-if(mode=='pilot') 'JBS_MC_PILOT_REPS' else 'JBS_MC_FULL_REPS'
  default<-if(mode=='pilot')50L else 1000L
  if(nzchar(Sys.getenv(specific,unset='')))
    return(v13_env_int(specific,default,20L))
  if(stage %in% c('mc_sequence','all_stages'))return(default)
  v13_env_int('JBS_MC_REPS',default,20L)
}
v13_mc_plan <- function(mode='pilot'){
  if(mode=='pilot') list(ns=c(250L,858L),
      scenarios=c('null','t2_local','t2_strong'), shapes=c('linear','nonlinear'),
      dependence=v13_dependence_levels(),
      reps=v13_mc_reps('pilot'))
  else list(ns=c(250L,500L,858L,1500L),
      scenarios=c('null','t2_local','t1_strong','t2_strong','tradeoff'),
      shapes=c('linear','nonlinear'),dependence=v13_dependence_levels(),
      reps=v13_mc_reps('full'))
}

v13_mc <- function(mode='pilot'){
  pl<-v13_mc_plan(mode)
  out<-v13_safe_output(file.path(CFG$out_dir,paste0('MONTE_CARLO_',toupper(mode))))
  writeLines(c(paste('RESEARCH MONTE CARLO:',mode),paste('Replications per cell:',pl$reps),
    paste('Sample sizes:',paste(pl$ns,collapse=',')),
    paste('Scenarios:',paste(pl$scenarios,collapse=',')),
    paste('Shapes:',paste(pl$shapes,collapse=',')),
    paste('Dependence values:',paste(pl$dependence,collapse=',')),
    'Truth is computed from KNOWN DGM; UCI sample not the truth.',
    'Failure-adjusted and conditional performance are separate.',
    'PV context simulations are separate and test a different estimand.',
    'Pilot outputs are not evidence of valid nominal coverage.',
    'No automatic journal-acceptance or global-theorem claim.'),
    file.path(out,'MC_SPECIFICATION.txt'))
  idx<-0L; summaries<-list();doms<-list();errors<-list()
  for(nn in pl$ns)for(sc in pl$scenarios)for(shape in pl$shapes)for(rho in pl$dependence){
    idx<-idx+1L;tag<-paste0('n',nn,'_',sc,'_',shape,'_rho',gsub('.', 'p',sprintf('%.6f',rho),fixed=TRUE))
    f1<-file.path(out,paste0('SUMMARY_',tag,'.csv'))
    f2<-file.path(out,paste0('DOMINANCE_',tag,'.csv'))
    f3<-file.path(out,paste0('RAW_',tag,'.csv'))
    if(all(file.exists(c(f1,f2,f3)))) {
      saved<-utils::read.csv(f3,stringsAsFactors=FALSE)
      if(nrow(saved)!=pl$reps)stop('Previous scenario has wrong rep count: ',tag,
                            '. Move/archive old scenario files before running')
      summaries[[tag]]<-utils::read.csv(f1);doms[[tag]]<-utils::read.csv(f2)
      v13_log('RESUME completed cell: ',tag)
      next
    }
    v13_log('MC START ',tag,' reps=',pl$reps)
    z<-tryCatch(run_mc_v11(n=nn,reps=pl$reps,scenario=sc,shape=shape,
       policy='age_dependent',tau=.0005,seed=20261301L+idx*41003L,
       dependence=rho),
       error=function(e)e)
    if(inherits(z,'error')) {
      msg<-conditionMessage(z)
      errors[[tag]]<-data.frame(cell=tag,error=msg)
      writeLines(msg,file.path(out,paste0('ERROR_',tag,'.txt')))
      v13_log('MC FAILED ',tag,': ',msg)
    }else{
      z<-v12_add_mc_fields(z)
      v13_write(z$raw,basename(f3),out)
      v13_write(z$summary,basename(f1),out)
      v13_write(z$dominance,basename(f2),out)
      summaries[[tag]]<-z$summary;doms[[tag]]<-z$dominance
      v13_write(do.call(rbind,summaries),'CHECKPOINT_MC_SUMMARY.csv',out)
      v13_log('MC SAVED ',tag)
    }
  }
  if(length(errors))v13_write(do.call(rbind,errors),'MC_CELL_ERRORS.csv',out)
  if(length(summaries)) v12_generate_mc_artifacts(out,summaries,doms,mode)
  complete<-length(summaries);expected<-length(pl$ns)*length(pl$scenarios)*length(pl$shapes)*length(pl$dependence)
  writeLines(c(paste('Successful scenario cells:',complete),
      paste('Expected scenario cells:',expected),
      paste('Skipped or failed:',expected-complete),
      'Check coverage MCSE, type-I, failures, nonregular null boundary before publication.'),
      file.path(out,'MC_COMPLETENESS_AUDIT.txt'))
  utils::capture.output(sessionInfo(),file=file.path(out,'R_SESSION_INFO.txt'))
  invisible(out)
}

# Explicit conservatism: reference PV paper and new Kappa cannot share H0.
v13_scope <- function(){
  data.frame(feature=c('Estimand','Covariate','Null hypothesis','Inference methods',
       'Data generating design','Software comparison','Published work relationship'),
   attached_manuscript=c(
      'Two tests: PPV and NPV (4-vector) and a global two-contrast test',
      'Categorical patterns (binary, 3-category simulations)',
      'PPV1=PPV2 AND NPV1=NPV2',
      'Regression models and categorical WLS; paper-specific global Wald',
      'Multinomial paired diagnostic cells within covariate patterns',
      'RM versus WLS versus no-covariate Wu method',
      'Existing published methodology; not the proposed novelty'),
   current_research=c(
      'Covariate-dependent loss-weighted kappa1, kappa2, Delta and B envelope',
      'X quantitative age; c(X) allowed to vary across individuals',
      'Delta=0 at policy; or robust B lower>0 / upper<0 across loss box',
      'Cross-fitted EIF, regression sandwich Delta, multiplier/local bootstrap',
      'Continuous-X paired Bernoulli with correlated test results',
      'Cross-fit versus regression Delta and sensitivity/bootstrap',
      'Separate extension; no claim to reproduce categorical PV paper'),
   stringsAsFactors=FALSE)
}

# Scientific readiness gate. Check actual MC outputs; computational success is insufficient.
v13_decision <- function(dir=CFG$out_dir,mode='full') {
  out<-v13_safe_output(dir)
  rows<-list()
  add<-function(area,passed,why,evidence='') {
    rows[[length(rows)+1L]]<<-data.frame(area=area,
      status=if(isTRUE(passed)) 'CHECK_PASSED' else 'NOT_YET_DEMONSTRATED',
      rationale=why,evidence=evidence,stringsAsFactors=FALSE)
  }
  baseline<-file.path(out,'06_Integrated_Weighted_Kappa_and_Paired_Comparison.csv')
  add('Real-data integrated kappa',file.exists(baseline),
      'Need validated real-data paired estimates',baseline)
  rec<-file.path(out,'V11_Table8_Regression_Delta_Sandwich.csv')
  # Some engine versions store this inside root; both paths checked.
  alt<-file.path(out,'JBS_Tables','V11_Table8_Regression_Delta_Sandwich.csv')
  add('Independent regression Delta',file.exists(rec)||file.exists(alt),
      'Independent model-based covariance comparison should be present',rec)
  mdir<-file.path(out,paste0('MONTE_CARLO_',toupper(mode)))
  sumfile<-file.path(mdir,'Table13_MC_Bias_RMSE_Coverage_Rejection.csv')
  completefile<-file.path(mdir,'MC_COMPLETENESS_AUDIT.txt')
  complete_ok<-FALSE
  if(file.exists(completefile)) {
    zz<-readLines(completefile,warn=FALSE)
    succ<-suppressWarnings(as.integer(sub('.*: *','',grep('^Successful scenario cells:',zz,value=TRUE))))
    expc<-suppressWarnings(as.integer(sub('.*: *','',grep('^Expected scenario cells:',zz,value=TRUE))))
    complete_ok<-length(succ)==1L && length(expc)==1L &&
       is.finite(succ) && is.finite(expc) && expc>0L && succ==expc
  }
  add('Monte Carlo scenario completeness',complete_ok,
    'No skipped simulation cells; check failure logs even when complete',completefile)
  if(!file.exists(sumfile)) {
    add('Monte Carlo: inferential operating characteristics',FALSE,
        'Missing real R MC output; do not report coverage or empirical type-I',sumfile)
  }else{
    tab<-utils::read.csv(sumfile)
    nrep<-if(nrow(tab))min(tab$total_reps,na.rm=TRUE)else 0L
    add('Monte Carlo replication precision',is.finite(nrep)&&nrep>=1000L,
        'Target >=1000 replications per scenario; MCSE must be displayed',sumfile)
    null<-tab[tab$delta_null_true & tab$method=='crossfit',,drop=FALSE]
    # This numeric threshold is a NON-BINDING inspection flag, not proof of asymptotics.
    ok_null<-nrow(null)>0 && all(is.finite(null$rejection_and_success)) &&
      all(null$rejection_and_success <= .08)
    add('Null type-I screen',ok_null,
        'Descriptive screen only: inspect MC confidence bounds and near-null/ties separately',sumfile)
    ok_cov<-nrow(tab)>0 && all(is.finite(tab$coverage_and_success)) &&
      all(tab$coverage_and_success >= .90)
    add('Coverage and fit-success screen',ok_cov,
        'Descriptive screen only: target 0.95, verify MCSE and conditional coverage',sumfile)
    failure<-nrow(tab)>0 && all(is.finite(tab$nominal_ci_valid_fraction)) &&
      all(tab$nominal_ci_valid_fraction >= .95)
    add('Estimator stability across DGM',failure,
        'No replacement of failed replications; check each failure reason',sumfile)
  }
  localfile<-file.path(out,'Table16_Local_Band_MC_Coverage.csv')
  add('Local conditional band calibration',file.exists(localfile),
      'Finite-grid local bootstrap requires its own repeated-sample coverage MC',localfile)
  add('Nonregular functional dominance uniformity',FALSE,
      'Fixed-tau conservative tests do not establish uniform size near ties',
      'New theoretical or calibrated method required')
  add('Clinician-elicited age-varying costs',FALSE,
      'Current c(X) scenarios are hypothetical; cannot be learned from UCI data')
  add('Biopsy verification and sampling validation',FALSE,
      'UCI dataset does not independently document unbiased reference verification')
  add('Literature novelty and coauthor approval',FALSE,
      'Requires documented prior-art review and authorship/ethics approval')
  z<-do.call(rbind,rows)
  v13_write(z,'V13_MANUSCRIPT_READINESS_CHECKLIST.csv',out)
  decision<-if(any(z$status!='CHECK_PASSED'))
       'NOT YET READY FOR SCIENTIFIC SUBMISSION: unmet and unverified criteria' else
       'ALL COMPUTATIONAL/PROTOCOL CHECKS PASSED; EDITORIAL ACCEPTANCE NOT GUARANTEED'
  writeLines(c(decision,'No software can verify novelty, assumptions, ethics or guarantee acceptance.',
          paste('Evaluation time:',format(Sys.time(),tz='UTC',usetz=TRUE))),
          file.path(out,'V13_RESEARCH_DECISION.txt'))
  print(z,row.names=FALSE)
  v13_log(decision)
  invisible(z)
}

v13_unit_tests <- function(){
  stopifnot(all(v12_algebra_audit()$passed),isTRUE(run_self_tests_v11()))
  # A strong paired PV null with non-identical observations must test H0 globally.
  dd<-v13_generate_pv(n=800L,categories=3L,scenario='alternative',seed=3381L)
  gg<-v13_pv_global(dd)
  stopifnot(is.finite(gg$p_global),gg$df==2L,gg$PPV_T1>=0,gg$PPV_T1<=1)
  fg<-v13_pv_categorical_context(data.frame(X=rep(20:59,length.out=800),
                             D=dd$D,T1=dd$T1,T2=dd$T2),2L)
  stopifnot(nrow(fg$by_stratum)==2L)
  tr<-v13_truth_pv(2L,'null')
  stopifnot(max(abs(tr))<1e-12)
  v13_log('V13 PASS: Kappa identities, paired PV global Wald, categorical DGM truth')
  TRUE
}

v13_run <- function(stage=tolower(Sys.getenv('JBS_STAGE',unset='selftest'))) {
  options(stringsAsFactors=FALSE,scipen=8)
  accepted<-c('selftest','protocol','analysis','generate_data','sensitivity',
             'reference_pv','mc_pilot','mc_full','pv_mc_pilot','pv_mc_full',
             'local_calibration','decision','all_pilot','all_full','all_stages',
             'mc_sequence')
  if(!stage %in% accepted)stop('Stage invalid: ',stage)
  # No Excel is requested in simulation-only modes.
  analysis_modes<-c('analysis','reference_pv','sensitivity',
                    'all_pilot','all_full','all_stages')
  v12_configure(if(stage %in% c('all_full','all_stages')) 'full' else
                if(stage %in% analysis_modes) 'analysis' else
                if(stage=='mc_full') 'mc_full' else if(stage=='local_calibration')
                'local_calibration' else 'mc_pilot')
  # Stable absolute result path; if analysis requires Excel, choose it FIRST.
  excel_env <- Sys.getenv('JBS_EXCEL_PATH',unset='')
  if(nzchar(excel_env)) CFG$excel <<- excel_env
  if(stage %in% analysis_modes){
    if(!file.exists(CFG$excel) && interactive()) {
      message('Select the verified cervical-cancer XLSX workbook')
      CFG$excel <<- file.choose()
    }
    if(!file.exists(CFG$excel))stop('Missing workbook; set JBS_EXCEL_PATH')
    CFG$excel <<- normalizePath(CFG$excel,winslash='/',mustWork=TRUE)
  }
  out_env <- Sys.getenv('JBS_OUTPUT_DIR',unset='')
  if(nzchar(out_env)) CFG$out_dir <<- out_env
  else if(stage %in% analysis_modes || nzchar(excel_env)) CFG$out_dir <<-
     file.path(dirname(if(stage %in% analysis_modes) CFG$excel else
          normalizePath(excel_env,winslash='/',mustWork=TRUE)),
          'JBS_V13_2_RESULTS')
  else CFG$out_dir <<- file.path(getwd(),'JBS_V13_2_RESULTS')
  CFG$out_dir <<- v13_safe_output(CFG$out_dir)
  # Write before checks to capture exact research distinctions.
  v13_write(v13_scope(),'V13_REFERENCE_VS_NEW_RESEARCH_SCOPE.csv')
  v13_write(v13_methodology_map(),'V13_METHODS_MAP_REFERENCE_VS_KAPPA.csv')
  v13_log('START stage=',stage,'  output=',CFG$out_dir)
  if(stage=='protocol') {
    v13_write(v12_claim_register(),'V13_CLAIMS_AND_PROOF_STATUS.csv')
    writeLines(c('CURRENT QUESTION: paired loss-weighted Kappa with c(X), not the published PPV global hypothesis.',
      'Inference with fixed-tau dominance smoothing is conservative, not uniform near ties.',
      'Clinical cost policy is hypothetical; group cutoffs used only for context.',
      'This is a protocol map, not proof of novelty or suitability for acceptance.'),
      file.path(CFG$out_dir,'READ_FIRST_DIFFERENT_RESEARCH_QUESTIONS.txt'))
    return(invisible(TRUE))
  }
  if(stage=='selftest'){
    v13_unit_tests();return(invisible(TRUE))
  }
  if(stage %in% analysis_modes && stage!='sensitivity'){
    # main_v10 performs actual analysis, figures 1-8, Tables 1-12, audit.
    main_v10()
    v12_write_research_status()
    # separate PV contextual analysis of same subjects, distinct question.
    df<-read_research_data(CFG$excel)
    pv<-v13_pv_global(df)
    v13_write(pv,'V13_Context_Paired_PV_Global_Wald.csv')
    v13_write(v13_no_x_kappa_context(df),
      'V13_Context_Unadjusted_Kappa_Different_Estimand.csv')
    for(K in 2:3){
      pp<-v13_pv_categorical_context(df,as.integer(K))
      v13_write(pp$by_stratum,paste0('V13_Context_PV_AgeGroups_',K,'.csv'))
    }
  }
  if(stage %in% c('sensitivity','all_pilot','all_full','all_stages')) {
    df_sens<-read_research_data(CFG$excel)
    v13_nuisance_sensitivity(df_sens,CFG$out_dir)
    v13_log('Nuisance linear/spline sensitivity audit saved')
  }
  if(stage %in% c('generate_data','all_stages')){
    for(sc in c('null','t2_local','t2_strong','tradeoff'))
      for(sh in c('linear','nonlinear'))
        v13_generate_demo(858L,sc,sh,CFG$out_dir,seed=33000L+match(sc,c('null','t2_local','t2_strong','tradeoff'))*37L+match(sh,c('linear','nonlinear')))
    for(K in 2:3)for(sc in c('null','alternative')){
      d<-v13_generate_pv(858L,as.integer(K),sc,seed=202610L+K+as.integer(sc=='alternative'))
      v13_write(d,paste0('SYNTHETIC_PV_cat',K,'_',sc,'.csv'),
        v13_safe_output(file.path(CFG$out_dir,'SYNTHETIC_DATA')))
    }
  }
  if(stage=='reference_pv'){
    # Stage just runs real-data analysis plus contextual PV summaries as requested.
  }
  if(stage %in% c('mc_pilot','mc_full','all_pilot','all_full',
                  'all_stages','mc_sequence')){
    modes <- if(stage %in% c('mc_sequence','all_stages')) c('pilot','full') else
       if(stage %in% c('mc_pilot','all_pilot')) 'pilot' else 'full'
    for(mode in modes){
      v13_mc(mode)
      pv_specific <- if(mode=='pilot') 'JBS_PV_PILOT_REPS' else 'JBS_PV_FULL_REPS'
      pv_def <- if(mode=='pilot')50L else 1000L
      pv_n <- if(nzchar(Sys.getenv(pv_specific,unset='')))
          v13_env_int(pv_specific,pv_def,20L) else
          if(stage %in% c('mc_sequence','all_stages')) pv_def else
          v13_env_int('JBS_PV_MC_REPS',pv_def,20L)
      v13_mc_pv(mode,pv_n,CFG$out_dir)
    }
  }
  if(stage %in% c('pv_mc_pilot','pv_mc_full')){
    mode<-if(stage=='pv_mc_pilot')'pilot' else 'full'
    v13_mc_pv(mode,v13_env_int('JBS_PV_MC_REPS',if(mode=='pilot')50L else 1000L,20L),CFG$out_dir)
  }
  if(stage=='local_calibration' ||
      (stage=='all_stages' && identical(Sys.getenv('JBS_RUN_LOCAL_CALIBRATION',unset='0'),'1'))){
    nrep<-v13_env_int('JBS_LOCAL_REPS',60L,20L)
    B<-v13_env_int('JBS_LOCAL_B',199L,50L)
    ll<-run_local_grid_calibration_v11(n=858L,reps=nrep,B=B,
         scenario='null',shape='nonlinear')
    v13_write(ll$summary,'Table16_Local_Band_MC_Coverage.csv')
    v13_write(ll$raw,'RAW_Local_Band_MC_Replications.csv')
  }
  if(stage %in% c('decision','mc_full','all_full','all_pilot',
                 'all_stages','mc_sequence'))
    v13_decision(CFG$out_dir,'full')
  v13_write(v13_scope(),'V13_REFERENCE_VS_NEW_RESEARCH_SCOPE.csv')
  utils::capture.output(sessionInfo(),file=file.path(CFG$out_dir,'V13_SESSION_INFO.txt'))
  v13_log('FINISHED stage=',stage,'; output=',CFG$out_dir)
  invisible(TRUE)
}

# ==============================================================================
# V13.3 ONE-CLICK FULL RESEARCH ORCHESTRATOR (PATCH OVER UNCHANGED V13.2)
# Default: source(this_file) runs ALL planned phases in one RStudio execution.
# Select the verified Excel ONCE when prompted. Requires long uninterrupted runtime.
# Main kappa and contextual PPV analysis are STATISTICALLY DISTINCT.
# For advanced single-stage debugging only:
#   Sys.setenv(JBS_ONE_CLICK="0", JBS_STAGE="analysis")
#   source(this_file, encoding="UTF-8")
# One-click defaults: B=999 real-data local bootstrap; 1999 multiplier draws;
# pilot MC=50; full MC=1000 per cell; pilot/full PV contextual MC;
# local coverage MC=200 independent samples with B=199 each (still preliminary).
# The one-click pipeline does NOT imply nominal statistical coverage, originality,
# clinical utility, or guaranteed journal acceptance.
# ==============================================================================

v133_default_env <- function(key, value) {
  # Ignore stale variables from earlier V12/V13.2 interactive sessions.
  # User overrides, if required, must use version-specific JBS_V133_* names.
  override_key <- paste0("JBS_V133_", sub("^JBS_", "", key))
  chosen <- Sys.getenv(override_key, unset=as.character(value))
  do.call(Sys.setenv, stats::setNames(list(chosen), key))
}

v133_mc_complete <- function(mode) {
  root <- file.path(CFG$out_dir, paste0("MONTE_CARLO_", toupper(mode)))
  f <- file.path(root, "MC_COMPLETENESS_AUDIT.txt")
  if (!file.exists(f)) stop("Missing MC audit: ", f)
  lines <- readLines(f, warn=FALSE)
  fetch <- function(prefix) {
    x <- grep(paste0("^", prefix), lines, value=TRUE)
    if (length(x)!=1L) stop("Malformed MC audit (", prefix, "): ", f)
    z <- suppressWarnings(as.integer(trimws(sub("^[^:]+:", "", x))))
    if (!is.finite(z)) stop("Non-numeric MC audit: ", x)
    z
  }
  complete <- fetch("Successful scenario cells:")
  expected <- fetch("Expected scenario cells:")
  if (expected<1L || complete!=expected)
    stop("MC ",mode," incomplete: ",complete,"/",expected," cells. Inspect MC_CELL_ERRORS.csv; saved cells remain reusable")
  invisible(list(completed=complete, expected=expected))
}

v133_output_check <- function(stage) {
  r <- CFG$out_dir
  req <- switch(stage,
    selftest=character(0),
    protocol="V13_METHODS_MAP_REFERENCE_VS_KAPPA.csv",
    analysis="06_Integrated_Weighted_Kappa_and_Paired_Comparison.csv",
    sensitivity="V13_Nuisance_Model_Sensitivity.csv",
    generate_data=file.path("SYNTHETIC_DATA", "SYNTHETIC_null_linear_n858.csv"),
    mc_pilot=file.path("MONTE_CARLO_PILOT", "MC_COMPLETENESS_AUDIT.txt"),
    mc_full=file.path("MONTE_CARLO_FULL", "MC_COMPLETENESS_AUDIT.txt"),
    local_calibration="Table16_Local_Band_MC_Coverage.csv",
    decision="V13_MANUSCRIPT_READINESS_CHECKLIST.csv",
    character(0))
  missing <- req[!file.exists(file.path(r,req))]
  if (length(missing)) stop("Phase ",stage," missing expected output: ",paste(missing,collapse=", "))
  if (stage %in% c("mc_pilot","mc_full")) {
    v133_mc_complete(if(stage=="mc_pilot") "pilot" else "full")
    pv_file<-file.path(r, paste0("REFERENCE_PV_MC_",toupper(if(stage=="mc_pilot") "pilot" else "full")),
      "TablePV_Context_MC_TypeI_Power_Failures.csv")
    if (!file.exists(pv_file)) stop("Contextual predictive-value MC table missing: ",pv_file)
  }
  invisible(TRUE)
}

v133_one_click <- function() {
  # Preserve an explicit user-provided configuration. Never invent patient data.
  oldstage <- Sys.getenv("JBS_STAGE",unset="")
  on.exit({
    if(nzchar(oldstage)) Sys.setenv(JBS_STAGE=oldstage) else Sys.unsetenv("JBS_STAGE")
  }, add=TRUE)
  Sys.setenv(JBS_STAGE="all_stages")
  excel <- Sys.getenv("JBS_EXCEL_PATH",unset="")
  if (!nzchar(excel) || !file.exists(excel)) {
    if (!interactive()) stop("Provide absolute JBS_EXCEL_PATH when executing noninteractively")
    message("ONE CLICK: Select Cervical_Cancer_858_Paired_Kappa_RStudio.xlsx ONCE")
    excel <- file.choose()
  }
  excel <- normalizePath(excel, winslash="/",mustWork=TRUE)
  Sys.setenv(JBS_EXCEL_PATH=excel)
  # Use a NEW results directory by default even when a V12 folder is still
  # defined in RStudio's environment; never silently overwrite prior results.
  out <- Sys.getenv("JBS_V133_OUTPUT_DIR",unset="")
  if (!nzchar(out)) out<-file.path(dirname(excel),"JBS_V13_3_RESULTS")
  out <- v13_safe_output(out)
  Sys.setenv(JBS_OUTPUT_DIR=out)
  CFG$out_dir <<- out
  CFG$excel <<- excel

  # These are full-research defaults, not performance guarantees.
  # Change before source() only if a different prespecified design is justified.
  v133_default_env("JBS_BOOTSTRAP_B", "999")
  v133_default_env("JBS_MULTIPLIER_B", "1999")
  v133_default_env("JBS_MC_PILOT_REPS", "50")
  v133_default_env("JBS_MC_FULL_REPS", "1000")
  v133_default_env("JBS_PV_PILOT_REPS", "50")
  v133_default_env("JBS_PV_FULL_REPS", "1000")
  v133_default_env("JBS_LOCAL_REPS", "200")
  v133_default_env("JBS_LOCAL_B", "199")
  v133_default_env("JBS_DEPENDENCE_LEVELS", "0.25,0.55,0.85")

  stages <- c("selftest", "protocol", "analysis", "sensitivity", "generate_data",
              "mc_pilot", "mc_full", "local_calibration", "decision")
  run_file <- file.path(out,"ONE_CLICK_STAGE_STATUS.csv")
  started_at <- Sys.time()
  status <- data.frame(stage=stages,status=rep("PENDING",length(stages)),
          start_time=rep("",length(stages)),end_time=rep("",length(stages)),
          elapsed_minutes=rep(NA_real_,length(stages)),
          warning_count=rep(NA_integer_,length(stages)),
          detail=rep("",length(stages)),stringsAsFactors=FALSE)
  save_status <- function() utils::write.csv(status,run_file,row.names=FALSE,na="",fileEncoding="UTF-8")
  save_status()

  env_keys<-c("JBS_EXCEL_PATH", "JBS_OUTPUT_DIR", "JBS_MC_PILOT_REPS", "JBS_MC_FULL_REPS",
     "JBS_PV_PILOT_REPS", "JBS_PV_FULL_REPS", "JBS_BOOTSTRAP_B", "JBS_MULTIPLIER_B",
     "JBS_LOCAL_REPS", "JBS_LOCAL_B", "JBS_DEPENDENCE_LEVELS")
  manifest<-data.frame(setting=env_keys, value=Sys.getenv(env_keys),stringsAsFactors=FALSE)
  utils::write.csv(manifest,file.path(out,"ONE_CLICK_PRESPECIFIED_CONFIGURATION.csv"),
                   row.names=FALSE,fileEncoding="UTF-8")
  utils::write.csv(data.frame(file=excel,md5=unname(tools::md5sum(excel))),
                   file.path(out,"ONE_CLICK_INPUT_MD5.csv"),row.names=FALSE)
  writeLines(c("V13.3 ONE-CLICK RESEARCH RUN",
    "The FULL option may take many hours or days depending on the computer.",
    "Stage status is saved after every stage. Monte Carlo cells are checkpointed.",
    "If rerunning after interruption, keep SAME replication counts and folder; completed MC cells resume.",
    "The local calibration uses 200 samples x 199 bootstrap draws by default; this is exploratory.",
    "Cross-fitting, model-based Delta, and PV-context estimands remain distinct.",
    "Completing the run cannot prove nonregular uniform size or journal acceptance."),
    file.path(out,"READ_FIRST_ONE_CLICK.txt"))
  logline("V13.3 ALL-STAGES single run; Excel=",excel)
  logline("Output=",out," ; primary full MC reps=",Sys.getenv("JBS_MC_FULL_REPS"),
          " ; dependence=",Sys.getenv("JBS_DEPENDENCE_LEVELS"))

  for (i in seq_along(stages)) {
    stage <- stages[i]
    Sys.setenv(JBS_STAGE="all_stages") # protects full-MC defaults during nested stages
    status$status[i] <- "RUNNING"
    started<-Sys.time()
    status$start_time[i]<-format(started,"%Y-%m-%d %H:%M:%S %Z")
    save_status()
    logline(sprintf("=== ONE CLICK PHASE %d/%d: %s ===",i,length(stages),stage))
    warns<-character(0)
    result<-tryCatch(
      withCallingHandlers({
        v13_run(stage)
        v133_output_check(stage)
        TRUE
      },warning=function(w) {
        warns <<- c(warns,conditionMessage(w))
        # Preserve the warning on the console; DO NOT silence research diagnostics.
      }),error=function(e)e)
    finished<-Sys.time()
    status$end_time[i]<-format(finished,"%Y-%m-%d %H:%M:%S %Z")
    status$elapsed_minutes[i]<-round(as.numeric(difftime(finished,started,units="mins")),3)
    status$warning_count[i]<-length(warns)
    if (inherits(result,"error")) {
      status$status[i]<-"FAILED"
      status$detail[i]<-conditionMessage(result)
      save_status()
      writeLines(c("FAILED PHASE: ",stage,"Error:",conditionMessage(result),
        "Already completed MC scenario cells are saved on disk.",
        "Rerun with identical settings to resume completed MC cells."),
        file.path(out,"ONE_CLICK_LAST_ERROR.txt"))
      logline("ONE CLICK STOPPED in phase=",stage," ; ",conditionMessage(result))
      stop("One-click run stopped in ",stage,": ",conditionMessage(result),call.=FALSE)
    }
    status$status[i]<-"COMPLETE"
    status$detail[i]<-if(length(warns))paste(unique(warns),collapse=" | ") else "No stage warnings"
    save_status()
    logline("PHASE COMPLETE: ",stage," ; elapsed minutes=",status$elapsed_minutes[i])
  }
  total<-round(as.numeric(difftime(Sys.time(),started_at,units="hours")),3)
  assessment<-file.path(out,"V13_RESEARCH_DECISION.txt")
  verdict<-if(file.exists(assessment))readLines(assessment,warn=FALSE)else "No research decision generated"
  writeLines(c("ALL COMPUTATIONAL STAGES COMPLETED",paste("Elapsed hours:",total),
     paste("Excel MD5:",unname(tools::md5sum(excel))),
     "Scientific caveats: review Table13, Table14, sample size, local coverage, and unproved nonregular inference.",
     "This is not a journal acceptance certificate.","",verdict),
     file.path(out,"ONE_CLICK_FINAL_SUMMARY.txt"))
  logline("ONE CLICK ALL PHASES COMPLETE; hours=",total,"; results=",out)
  invisible(status)
}

# Default source() action: run everything once, without extra commands.
# To run only one stage for debugging, set JBS_ONE_CLICK='0' and JBS_STAGE.
if (Sys.getenv("JBS_ONE_CLICK",unset="1") != "0") {
  tryCatch(v133_one_click(),error=function(e) {
    message("V13.3 ONE-CLICK STOPPED: ",conditionMessage(e))
    stop(e)
  })
} else {
  tryCatch(v13_run(),error=function(e) {
    message("V13.3 SINGLE STAGE STOPPED: ",conditionMessage(e))
    stop(e)
  })
}
