################################################################################
# JBS_V18_SECOND_APPLICATION_AND_SENSITIVITY.R
#
# Covariate-dependent loss-weighted kappa for paired binary diagnostic tests.
# Analyses added in manuscript version V18 (R only):
#   PART A  Sections 6.6-6.7 / Table 9, Tables S5-S7  role of age and further
#           sensitivity analyses, cervical data
#   PART B  Section 7 / Tables 10-12, S8, Figure 5a  coronary application
#   PART C  Section 7.6 / Table 13, Table S9, Figure 5b  calibrated simulation
#   PART D  Section 7.7 / Appendix M, Table S10  categorical covariate (colorectal)
#
# Author : Ismail Radi Ismail Fkhadah, Department of Biostatistics, School of
#          Medicine, University of Granada, Spain (ismailfk2021@correo.ugr.es)
# License: MIT
#
# HOW TO RUN (RStudio). Set the working directory to the repository root, then
#   Sys.setenv(JBS_V18_PARTS = "cervical,heart,crc,sim")   # any subset
#   Sys.setenv(JBS_V18_OUT   = "new_results_V18")          # output folder
#   source("R/JBS_V18_SECOND_APPLICATION_AND_SENSITIVITY.R", encoding = "UTF-8")
# PART C (4 x 1,000 replicates) takes about 20-60 minutes; PARTs A, B, D take
# about one minute. Optional: JBS_V18_SIM_REPS (default 1000) and
# JBS_V18_SIM_NS (default "301,602,1204,2408").
#
# The estimators, folds, policies, box, endpoint interval, smoothing values and
# multiplier bootstrap are those of R/JBS_V13_3_ONE_CLICK_ALL_STAGES.R (the
# functions below reproduce fixed_policy_inference(), endpoint_iut_v11(),
# prepare_envelope(), smooth_dominance_test(), upper_dominance_test_v11(),
# multiplier_integrated_band_v11() and regression_delta_v11()). Nuisance models
# are bias-reduced logistic regressions on splines::ns(X, df = 3), fitted with
# brglm2 (type = "AS_mean") when the package is installed, as in the main
# program; otherwise an equivalent built-in Firth fit is used (for the logistic
# model, mean bias reduction is Firth's penalized likelihood; the two agree to
# the convergence tolerance, about 1e-9 in the predicted probabilities).
################################################################################

options(stringsAsFactors = FALSE, scipen = 6)
V18 <- list(
  seed = 20261017L, K = 5L, df = 3L, alpha = 0.05,
  taus = c(0.0001, 0.0005, 0.002), box = c(0.40, 0.90), endpoints = c(0.10, 0.90),
  policy_center = 35, policy_c35 = 0.65, policy_slope = 0.025,
  multiplier_B = 1999L, grid = seq(0.05, 0.95, by = 0.025),
  parts = strsplit(Sys.getenv("JBS_V18_PARTS", "cervical,heart,crc,sim"), ",")[[1]],
  out = Sys.getenv("JBS_V18_OUT", "new_results_V18"),
  sim_reps = as.integer(Sys.getenv("JBS_V18_SIM_REPS", "1000")),
  sim_ns = as.integer(strsplit(Sys.getenv("JBS_V18_SIM_NS", "301,602,1204,2408"), ",")[[1]]),
  use_brglm2 = requireNamespace("brglm2", quietly = TRUE) &&
               !identical(Sys.getenv("JBS_V18_FORCE_INTERNAL_FIRTH"), "1")
)
dir.create(V18$out, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
say("nuisance fitter:", if (V18$use_brglm2) "brglm2::brglmFit (AS_mean)" else "built-in Firth (equivalent)")

# ------------------------------------------------------------------ helpers
expit <- function(z) stats::plogis(z)
clip01 <- function(p, eps = 1e-7) pmax(eps, pmin(1 - eps, as.numeric(p)))
age_policy <- function(x, slope = V18$policy_slope, center = V18$policy_center, c35 = V18$policy_c35)
  stats::plogis(stats::qlogis(c35) + slope * (x - center))

# Bias-reduced logistic regression of y on the design matrix B (B includes the intercept).
firth_fit <- function(B, y) {
  if (V18$use_brglm2) {
    fit <- suppressWarnings(stats::glm(y ~ 0 + B, family = stats::binomial(),
             method = brglm2::brglmFit, type = "AS_mean",
             control = stats::glm.control(maxit = 200L)))
    b <- stats::coef(fit)
    if (!isTRUE(fit$converged) || any(!is.finite(b))) stop("bias-reduced fit did not converge")
    return(unname(b))
  }
  b <- rep(0, ncol(B))
  pen <- function(b) { eta <- drop(B %*% b); mu <- expit(eta)
    sum(y * eta - log1p(exp(eta))) + 0.5 * as.numeric(determinant(crossprod(B, B * (mu * (1 - mu))))$modulus) }
  cur <- pen(b)
  for (it in 1:500) {
    mu <- expit(drop(B %*% b)); w <- mu * (1 - mu)
    inv <- solve(crossprod(B, B * w))
    h <- rowSums((B %*% inv) * B) * w
    step <- drop(inv %*% crossprod(B, y - mu + h * (0.5 - mu)))
    t <- 1
    repeat { nb <- b + t * step; nw <- pen(nb); if (nw >= cur - 1e-12 || t < 1e-8) break; t <- t / 2 }
    b <- nb; cur <- nw
    if (max(abs(t * step)) < 1e-12) break
  }
  b
}

nuisance_design <- function(model, xtr, xte, df = V18$df) {
  if (model == "spline") { bs <- splines::ns(xtr, df = df)
    return(list(tr = cbind(1, bs), te = cbind(1, stats::predict(bs, xte)))) }
  if (model == "intercept") return(list(tr = matrix(1, length(xtr), 1), te = matrix(1, length(xte), 1)))
  if (model == "factor") { lev <- sort(unique(xtr))[-1]
    ind <- function(x) sapply(lev, function(l) as.numeric(x == l))
    return(list(tr = cbind(1, ind(xtr)), te = cbind(1, ind(xte)))) }
  stop("unknown nuisance model")
}

# Cross-fitted p, q1, q2 with the fold rule of crossfit_probabilities().
crossfit <- function(X, D, T1, T2, K = V18$K, seed = V18$seed, model = "spline", df = V18$df, fold = NULL) {
  n <- length(X)
  if (is.null(fold)) { set.seed(seed); fold <- sample(rep(seq_len(K), length.out = n), size = n, replace = FALSE) }
  pred <- matrix(NA_real_, n, 3); ev <- matrix(NA_integer_, K, 3)
  for (k in seq_len(K)) {
    id <- which(fold == k); des <- nuisance_design(model, X[-id], X[id], df)
    for (j in 1:3) { y <- list(D, T1, T2)[[j]]
      pred[id, j] <- clip01(expit(drop(des$te %*% firth_fit(des$tr, y[-id]))))
      ev[k, j] <- sum(y[-id]) }
  }
  list(p = pred[, 1], q1 = pred[, 2], q2 = pred[, 3], fold = fold, train_events = ev)
}

# Fixed policy: A_j, C_j, kappa_j, Delta, paired influence functions (fixed_policy_inference()).
fixed_policy <- function(D, T1, T2, nu, cval, s = rep(1, length(D)), alpha = V18$alpha) {
  n <- length(D); cval <- rep_len(cval, n)
  a1 <- s * (cval * D * (1 - T1) + (1 - cval) * (1 - D) * T1)
  a2 <- s * (cval * D * (1 - T2) + (1 - cval) * (1 - D) * T2)
  u1 <- s * (D - nu$p) * (T1 - nu$q1); u2 <- s * (D - nu$p) * (T2 - nu$q2)
  A1 <- mean(a1); A2 <- mean(a2); C1 <- mean(u1); C2 <- mean(u2); R1 <- A1 + C1; R2 <- A2 + C2
  if (min(R1, R2) <= 1e-8) stop("non-positive reference loss")
  phi1 <- (A1 * (u1 - C1) - C1 * (a1 - A1)) / R1^2; phi2 <- (A2 * (u2 - C2) - C2 * (a2 - A2)) / R2^2
  phid <- phi1 - phi2; se <- function(f) stats::sd(f) / sqrt(n); z <- stats::qnorm(1 - alpha / 2)
  d <- C1 / R1 - C2 / R2; L <- a1 - a2
  list(A1 = A1, A2 = A2, C1 = C1, C2 = C2, R1 = R1, R2 = R2, kappa1 = C1 / R1, kappa2 = C2 / R2,
       delta = d, se_kappa1 = se(phi1), se_kappa2 = se(phi2), se_delta = se(phid),
       ci_low = d - z * se(phid), ci_high = d + z * se(phid), p = 2 * stats::pnorm(-abs(d / se(phid))),
       loss_diff = mean(L), se_loss = se(L - mean(L)), loss_low = mean(L) - z * se(L - mean(L)),
       loss_high = mean(L) + z * se(L - mean(L)), p_loss = 2 * stats::pnorm(-abs(mean(L) / se(L - mean(L)))),
       phi_delta = phid)
}
row_of <- function(r, label) data.frame(label = label, kappa1 = r$kappa1, kappa2 = r$kappa2, delta = r$delta,
  se = r$se_delta, ci_low = r$ci_low, ci_high = r$ci_high, p = r$p, C1 = r$C1, C2 = r$C2,
  loss_diff = r$loss_diff, loss_low = r$loss_low, loss_high = r$loss_high, p_loss = r$p_loss)

endpoint_iut <- function(D, T1, T2, nu, cL = V18$endpoints[1], cU = V18$endpoints[2]) {
  e <- lapply(c(cL, cU), function(cc) fixed_policy(D, T1, T2, nu, cc))
  z <- sapply(e, function(x) x$delta / x$se_delta)
  data.frame(cL = cL, cU = cU, delta_L = e[[1]]$delta, se_L = e[[1]]$se_delta, delta_U = e[[2]]$delta,
             se_U = e[[2]]$se_delta, p_IUT_T1 = max(stats::pnorm(z, lower.tail = FALSE)), p_IUT_T2 = max(stats::pnorm(z)))
}

envelope <- function(D, T1, T2, nu, L, U, s = rep(1, length(D))) {
  n <- length(D); L <- rep_len(L, n); U <- rep_len(U, n); p <- nu$p; q1 <- nu$q1; q2 <- nu$q2
  eD <- D - p; e1 <- T1 - q1; e2 <- T2 - q2; uC1 <- s * eD * e1; uC2 <- s * eD * e2
  C1 <- mean(uC1); C2 <- mean(uC2)
  v1 <- s * ((1 - p) * q1 - q1 * eD + (1 - p) * e1); v2 <- s * ((1 - p) * q2 - q2 * eD + (1 - p) * e2)
  V1 <- mean(v1); V2 <- mean(v2); h <- s * (C1 * (p - q2) - C2 * (p - q1)); b0 <- C1 * V2 - C2 * V1
  d1 <- s * (p - q1); d2 <- s * (p - q2)
  list(s = s, L = L, U = U, p = p, q1 = q1, q2 = q2, C1 = C1, C2 = C2, V1 = V1, V2 = V2,
       phiC1 = uC1 - C1, phiC2 = uC2 - C2, phiV1 = v1 - V1, phiV2 = v2 - V2, h = h, b0 = b0,
       lower = b0 + mean(pmin(L * h, U * h)), upper = b0 + mean(pmax(L * h, U * h)),
       minR1 = V1 + mean(pmin(L * d1, U * d1)), minR2 = V2 + mean(pmin(L * d2, U * d2)),
       xi = (C1 - C2) * eD + C2 * e1 - C1 * e2)
}
# Six fixed-smoothing tests (3 tau x 2 directions), each one-sided at alpha/6.
robust_tests <- function(env, taus = V18$taus, alpha = V18$alpha) {
  a_each <- alpha / (2 * length(taus)); zc <- stats::qnorm(1 - a_each); n <- length(env$h)
  out <- list()
  for (tau in taus) for (dir in c("T1_over_T2", "T2_over_T1")) {
    pr <- stats::plogis((env$U - env$L) * env$h / tau)
    if (dir == "T1_over_T2") {
      a <- -env$L * env$h / tau; b <- -env$U * env$h / tau; m <- pmax(a, b)
      g <- -tau * (m + log(exp(a - m) + exp(b - m))); w <- env$L * pr + env$U * (1 - pr)
    } else {
      a <- env$L * env$h / tau; b <- env$U * env$h / tau; m <- pmax(a, b)
      g <- tau * (m + log(exp(a - m) + exp(b - m))); w <- env$L * (1 - pr) + env$U * pr
    }
    al1 <- mean(w * env$s * (env$p - env$q2)); al2 <- mean(w * env$s * (env$p - env$q1))
    theta <- env$b0 + mean(g + w * env$s * env$xi)
    phi <- (env$V2 + al1) * env$phiC1 - (env$V1 + al2) * env$phiC2 + env$C1 * env$phiV2 -
           env$C2 * env$phiV1 + g - mean(g) + w * env$s * env$xi
    phi <- phi - mean(phi); se <- stats::sd(phi) / sqrt(n)
    ok <- min(env$minR1, env$minR2) > 1e-8
    bound <- if (dir == "T1_over_T2") theta - zc * se else theta + zc * se
    rej <- ok && (if (dir == "T1_over_T2") bound > 0 else bound < 0)
    out[[length(out) + 1]] <- data.frame(direction = dir, tau = tau, theta = theta, se = se, bound = bound, certificate = rej)
  }
  do.call(rbind, out)
}

multiplier_band <- function(D, T1, T2, nu, grid = V18$grid, B = V18$multiplier_B, seed = V18$seed + 111L, alpha = V18$alpha) {
  n <- length(D); obj <- lapply(grid, function(cv) fixed_policy(D, T1, T2, nu, cv))
  mat <- do.call(cbind, lapply(obj, `[[`, "phi_delta")); se_cols <- apply(mat, 2, stats::sd)
  set.seed(seed); crit <- numeric(B); block <- 250L
  for (start in seq.int(1L, B, by = block)) {
    m <- min(block, B - start + 1L); noise <- matrix(stats::rnorm(n * m), nrow = n, ncol = m)
    crit[start:(start + m - 1L)] <- apply(abs(sweep(crossprod(noise, mat), 2, sqrt(n) * se_cols, "/")), 1, max)
  }
  cv <- unname(stats::quantile(crit, 1 - alpha)); dl <- sapply(obj, `[[`, "delta"); se <- se_cols / sqrt(n)
  z <- stats::qnorm(1 - alpha / 2)
  data.frame(c = grid, delta = dl, se = se, pointwise_low = dl - z * se, pointwise_high = dl + z * se,
             simultaneous_low = dl - cv * se, simultaneous_high = dl + cv * se, critical = cv)
}

# Regression-standardized benchmark with stacked subject-level sandwich (regression_delta_v11()).
regression_benchmark <- function(X, D, T1, T2, policies, s = rep(1, length(X)), alpha = V18$alpha) {
  n <- length(X); Z <- (X - mean(X)) / stats::sd(X)
  XD <- cbind(1, Z); XT <- cbind(1, D, Z, D * Z); designs <- list(XD, XT, XT); ys <- list(D, T1, T2)
  fits <- lapply(1:3, function(j) {
    f <- suppressWarnings(stats::glm.fit(designs[[j]], ys[[j]], family = stats::binomial(),
                                         control = stats::glm.control(maxit = 150L, epsilon = 1e-9)))
    if (!isTRUE(f$converged) || f$rank != ncol(designs[[j]]) || any(!is.finite(f$coefficients)))
      stop("regression benchmark not estimable")
    f })
  beta <- unlist(lapply(fits, `[[`, "coefficients"), use.names = FALSE)
  pr <- lapply(1:3, function(j) expit(drop(designs[[j]] %*% fits[[j]]$coefficients)))
  S <- do.call(cbind, lapply(1:3, function(j) designs[[j]] * (ys[[j]] - pr[[j]])))
  cuts <- cumsum(c(0, sapply(designs, ncol))); H <- matrix(0, length(beta), length(beta))
  for (j in 1:3) { idx <- (cuts[j] + 1):cuts[j + 1]
    H[idx, idx] <- crossprod(designs[[j]], designs[[j]] * (pr[[j]] * (1 - pr[[j]]))) / n }
  Hinv <- solve(H); X0 <- cbind(1, 0, Z, 0); X1 <- cbind(1, 1, Z, Z)
  comp <- function(par, cv) {
    p <- expit(drop(XD %*% par[1:2])); t10 <- expit(drop(X0 %*% par[3:6])); t11 <- expit(drop(X1 %*% par[3:6]))
    t20 <- expit(drop(X0 %*% par[7:10])); t21 <- expit(drop(X1 %*% par[7:10]))
    cbind(s * (cv * p * (1 - t11) + (1 - cv) * (1 - p) * t10), s * p * (1 - p) * (t11 - t10),
          s * (cv * p * (1 - t21) + (1 - cv) * (1 - p) * t20), s * p * (1 - p) * (t21 - t20)) }
  out <- list()
  for (nm in names(policies)) {
    cv <- rep_len(policies[[nm]], n); pp <- comp(beta, cv); th <- colMeans(pp)
    R1 <- th[1] + th[2]; R2 <- th[3] + th[4]; delta <- th[2] / R1 - th[4] / R2
    G <- matrix(NA_real_, 4, length(beta))
    for (m in seq_along(beta)) { st <- 1e-5 * max(1, abs(beta[m])); bp <- bm <- beta
      bp[m] <- bp[m] + st; bm[m] <- bm[m] - st
      G[, m] <- (colMeans(comp(bp, cv)) - colMeans(comp(bm, cv))) / (2 * st) }
    phi_t <- sweep(pp, 2, th, "-") + S %*% t(Hinv) %*% t(G)
    phd <- drop(phi_t %*% c(-th[2] / R1^2, th[1] / R1^2, th[4] / R2^2, -th[3] / R2^2))
    se <- stats::sd(phd) / sqrt(n); z <- stats::qnorm(1 - alpha / 2)
    out[[nm]] <- data.frame(label = nm, kappa1 = th[2] / R1, kappa2 = th[4] / R2, delta = unname(delta), se = se,
                            ci_low = delta - z * se, ci_high = delta + z * se, p = 2 * stats::pnorm(-abs(delta / se)))
  }
  do.call(rbind, out)
}
marginal_wk <- function(D, T, cc) { p <- mean(D); q <- mean(T)
  1 - mean(cc * D * (1 - T) + (1 - cc) * (1 - D) * T) / (cc * p * (1 - q) + (1 - cc) * (1 - p) * q) }
cp_ci <- function(k, m) c(if (k > 0) stats::qbeta(.025, k, m - k + 1) else 0, if (k < m) stats::qbeta(.975, k + 1, m - k) else 1)
wcsv <- function(x, f) utils::write.csv(x, file.path(V18$out, f), row.names = FALSE)
box_list <- list(c(0.40, 0.90), c(0.30, 0.90), c(0.50, 0.90), c(0.40, 0.80), c(0.25, 0.75), c(0.50, 0.95), c(0.10, 0.90))
box_row <- function(D, T1, T2, nu, L, U, label) {
  e <- envelope(D, T1, T2, nu, L, U); rr <- robust_tests(e)
  data.frame(box = label, lower_B = e$lower, upper_B = e$upper, minR1 = e$minR1, minR2 = e$minR2,
             smallest_UCB_T2 = min(rr$bound[rr$direction == "T2_over_T1"]),
             largest_LCB_T1 = max(rr$bound[rr$direction == "T1_over_T2"]),
             certificate_T1 = any(rr$certificate[rr$direction == "T1_over_T2"]),
             certificate_T2 = any(rr$certificate[rr$direction == "T2_over_T1"]))
}

################################################################################
# PART A -- cervical data: fold seeds, ten folds, policy parameters, boxes, intervals
################################################################################
if ("cervical" %in% V18$parts) {
  say("PART A: cervical sensitivity analyses")
  d <- utils::read.csv("data/cervical_858_age_biopsy_hinselmann_schiller.csv")
  X <- d$X; D <- d$D; T1 <- d$T1; T2 <- d$T2
  nu <- crossfit(X, D, T1, T2)
  rows <- list(row_of(fixed_policy(D, T1, T2, nu, age_policy(X)), "primary: slope 0.025, centre 35, c(35)=0.65, K=5"))
  rows[[2]] <- row_of(fixed_policy(D, T1, T2, crossfit(X, D, T1, T2, K = 10L), age_policy(X)), "K = 10 folds")
  for (sl in c(0, 0.0125, 0.05)) rows[[length(rows) + 1]] <- row_of(fixed_policy(D, T1, T2, nu, age_policy(X, slope = sl)), paste("slope", sl))
  for (ce in c(25, 45)) rows[[length(rows) + 1]] <- row_of(fixed_policy(D, T1, T2, nu, age_policy(X, center = ce)), paste("centre", ce))
  for (c0 in c(0.55, 0.75)) rows[[length(rows) + 1]] <- row_of(fixed_policy(D, T1, T2, nu, age_policy(X, c35 = c0)), paste("c(35) =", c0))
  wcsv(do.call(rbind, rows), "A1_cervical_policy_sensitivity.csv")
  seeds <- V18$seed + 1000L * (1:50)
  dd <- sapply(seeds, function(sd) fixed_policy(D, T1, T2, crossfit(X, D, T1, T2, seed = sd), age_policy(X))$delta)
  wcsv(data.frame(seed = seeds, delta = dd), "A2_cervical_fold_seed_sensitivity.csv")
  bx <- do.call(rbind, lapply(box_list, function(b) box_row(D, T1, T2, nu, b[1], b[2], sprintf("[%.2f,%.2f]", b[1], b[2]))))
  ca <- age_policy(X); bx <- rbind(bx, box_row(D, T1, T2, nu, pmax(ca - .15, 0), pmin(ca + .15, 1), "c(x) +/- 0.15"))
  wcsv(bx, "A3_cervical_box_sensitivity.csv")
  wcsv(do.call(rbind, lapply(c(.10, .20, .25, .30, .40), function(a) endpoint_iut(D, T1, T2, nu, a, .90))), "A4_cervical_endpoint_intervals.csv")
  # Section 6.6: information carried by age (Table 9), age-free nuisances (Table S5),
  # age distribution and age-restricted analyses (Table S6). Age-free nuisances are
  # cross-fitted intercept-only bias-reduced fits on the same folds.
  nu0 <- crossfit(X, D, T1, T2, model = "intercept", fold = nu$fold)
  z <- (X - mean(X)) / stats::sd(X)
  info <- do.call(rbind, lapply(list(list("biopsy", D, nu$p, nu0$p), list("Hinselmann", T1, nu$q1, nu0$q1),
                                     list("Schiller", T2, nu$q2, nu0$q2)), function(o) {
    y <- o[[2]]; f0 <- stats::glm(y ~ 1, family = stats::binomial()); f1 <- stats::glm(y ~ z, family = stats::binomial())
    f3 <- stats::glm(y ~ z + I(z^2) + I(z^3), family = stats::binomial())
    ll <- function(f) -mean(y * log(f) + (1 - y) * log(1 - f))
    data.frame(outcome = o[[1]], p_linear = stats::pchisq(2 * as.numeric(stats::logLik(f1) - stats::logLik(f0)), 1, lower.tail = FALSE),
               p_cubic = stats::pchisq(2 * as.numeric(stats::logLik(f3) - stats::logLik(f0)), 3, lower.tail = FALSE),
               OR_per_10y = exp(unname(stats::coef(f1)[2]) / stats::sd(X) * 10),
               brier_spline = mean((y - o[[3]])^2), brier_intercept = mean((y - o[[4]])^2),
               logloss_spline = ll(o[[3]]), logloss_intercept = ll(o[[4]]),
               observed = mean(y), mean_prediction = mean(o[[3]])) }))
  wcsv(info, "A5_cervical_age_information.csv")
  pols <- list("equal, c = 0.50" = 0.5, "FN-weighted, c = 0.75" = 0.75, "age-dependent c(x)" = age_policy(X),
               "constant c = mean of c(x)" = mean(age_policy(X)))
  wcsv(do.call(rbind, lapply(names(pols), function(k) { a <- fixed_policy(D, T1, T2, nu, pols[[k]]); b <- fixed_policy(D, T1, T2, nu0, pols[[k]])
    data.frame(policy = k, delta_age_adjusted = a$delta, se_age_adjusted = a$se_delta, delta_age_free = b$delta,
               se_age_free = b$se_delta, difference = a$delta - b$delta) })), "A6_cervical_age_free_nuisances.csv")
  rest <- lapply(list("all ages" = rep(TRUE, length(X)), "age <= 45" = X <= 45, "age <= 40" = X <= 40), function(m) {
    nm <- crossfit(X[m], D[m], T1[m], T2[m], model = "intercept", fold = nu$fold[m])
    r <- fixed_policy(D[m], T1[m], T2[m], nm, age_policy(X[m]))
    data.frame(n = sum(m), biopsy_positive = sum(D[m]), discordant = sum(T1[m] != T2[m]),
               Lambda = r$loss_diff, Lambda_low = r$loss_low, Lambda_high = r$loss_high,
               delta_age_free = r$delta, delta_low = r$ci_low, delta_high = r$ci_high) })
  grp <- lapply(list("<= 30" = X <= 30, "31-45" = X > 30 & X <= 45, "> 45" = X > 45), function(m)
    data.frame(n = sum(m), biopsy_positive = sum(D[m]), discordant = sum(T1[m] != T2[m]), Lambda = NA, Lambda_low = NA,
               Lambda_high = NA, delta_age_free = NA, delta_low = NA, delta_high = NA))
  tab <- do.call(rbind, c(rest, grp)); tab <- cbind(population = c(names(rest), paste("age group", names(grp))), tab)
  wcsv(tab, "A7_cervical_age_restricted.csv")
  say("PART A done; fold-seed range", sprintf("%.4f to %.4f", min(dd), max(dd)))
}

################################################################################
# PART B -- coronary application (UCI Heart Disease, Cleveland; n = 301)
################################################################################
if ("heart" %in% V18$parts) {
  say("PART B: coronary application")
  d <- utils::read.csv("data/heart_cleveland_301_extended.csv")
  X <- d$X; D <- d$D; T1 <- d$T1; T2 <- d$T2; n <- nrow(d)
  POL <- list(equal = 0.5, fn_075 = 0.75, age_dependent = age_policy(X))
  # descriptive statistics and unadjusted accuracy
  desc <- do.call(rbind, lapply(list(all = rep(TRUE, n), no_CAD = D == 0, CAD = D == 1), function(m)
    data.frame(n = sum(m), age_mean = mean(X[m]), age_sd = stats::sd(X[m]), age_median = stats::median(X[m]),
               age_q1 = stats::quantile(X[m], .25), age_q3 = stats::quantile(X[m], .75), male = sum(d$sex_male[m]),
               T1_pos = sum(T1[m]), T2_pos = sum(T2[m]))))
  desc$group <- c("all", "no_CAD", "CAD"); wcsv(desc, "B1_heart_descriptives.csv")
  acc <- do.call(rbind, lapply(list(T1 = T1, T2 = T2), function(T) {
    tp <- sum(T == 1 & D == 1); fn <- sum(T == 0 & D == 1); tn <- sum(T == 0 & D == 0); fp <- sum(T == 1 & D == 0)
    data.frame(Se = tp / (tp + fn), Se_low = cp_ci(tp, tp + fn)[1], Se_high = cp_ci(tp, tp + fn)[2],
               Sp = tn / (tn + fp), Sp_low = cp_ci(tn, tn + fp)[1], Sp_high = cp_ci(tn, tn + fp)[2],
               PPV = tp / (tp + fp), PPV_low = cp_ci(tp, tp + fp)[1], PPV_high = cp_ci(tp, tp + fp)[2],
               NPV = tn / (tn + fn), NPV_low = cp_ci(tn, tn + fn)[1], NPV_high = cp_ci(tn, tn + fn)[2], FN = fn, FP = fp) }))
  acc$test <- c("exercise ECG", "thallium")
  mc <- function(m) { b10 <- sum(T1 == 1 & T2 == 0 & m); b01 <- sum(T1 == 0 & T2 == 1 & m)
    c(b10, b01, stats::binom.test(b10, b10 + b01)$p.value) }
  acc$mcnemar_T1only <- c(mc(D == 1)[1], mc(D == 0)[1]); acc$mcnemar_T2only <- c(mc(D == 1)[2], mc(D == 0)[2])
  acc$mcnemar_p <- c(mc(D == 1)[3], mc(D == 0)[3]); acc$mcnemar_for <- c("Se (diseased)", "Sp (non-diseased)")
  wcsv(acc, "B2_heart_accuracy.csv")
  # information carried by age
  nu <- crossfit(X, D, T1, T2); nu0 <- crossfit(X, D, T1, T2, model = "intercept", fold = nu$fold)
  z <- (X - mean(X)) / stats::sd(X)
  info <- do.call(rbind, lapply(list(D = list(D, nu$p, nu0$p), T1 = list(T1, nu$q1, nu0$q1), T2 = list(T2, nu$q2, nu0$q2)), function(o) {
    y <- o[[1]]; f0 <- stats::glm(y ~ 1, family = stats::binomial()); f1 <- stats::glm(y ~ z, family = stats::binomial())
    f3 <- stats::glm(y ~ z + I(z^2) + I(z^3), family = stats::binomial())
    ll <- function(f) -mean(y * log(f) + (1 - y) * log(1 - f))
    data.frame(p_linear = stats::pchisq(2 * (stats::logLik(f1) - stats::logLik(f0)), 1, lower.tail = FALSE),
               p_cubic = stats::pchisq(2 * (stats::logLik(f3) - stats::logLik(f0)), 3, lower.tail = FALSE),
               OR_per_10y = exp(stats::coef(f1)[2] / stats::sd(X) * 10),
               brier_spline = mean((y - o[[2]])^2), brier_intercept = mean((y - o[[3]])^2),
               logloss_spline = ll(o[[2]]), logloss_intercept = ll(o[[3]]), observed = mean(y), mean_prediction = mean(o[[2]])) }))
  info$outcome <- c("CAD", "exercise ECG", "thallium"); wcsv(info, "B3_heart_age_information.csv")
  wcsv(data.frame(X = X, D = D, T1 = T1, T2 = T2, fold = nu$fold, p_hat = nu$p, q1_hat = nu$q1, q2_hat = nu$q2), "B0_heart_out_of_fold_predictions.csv")
  # policy-specific estimates: age-conditional, age-free, regression benchmark
  est <- rbind(do.call(rbind, lapply(names(POL), function(k) row_of(fixed_policy(D, T1, T2, nu, POL[[k]]), paste("age-conditional", k)))),
               do.call(rbind, lapply(names(POL), function(k) row_of(fixed_policy(D, T1, T2, nu0, POL[[k]]), paste("age-free", k)))))
  wcsv(est, "B4_heart_policy_estimates.csv")
  rb <- regression_benchmark(X, D, T1, T2, POL); wcsv(rb, "B5_heart_regression_benchmark.csv")
  cov <- data.frame(C1 = est$C1[1], C2 = est$C2[1], marginal_cov1 = mean((D - mean(D)) * (T1 - mean(T1))),
                    marginal_cov2 = mean((D - mean(D)) * (T2 - mean(T2))))
  cov$between1 <- cov$marginal_cov1 - cov$C1; cov$between2 <- cov$marginal_cov2 - cov$C2
  cov$share1 <- cov$between1 / cov$marginal_cov1; cov$share2 <- cov$between2 / cov$marginal_cov2
  cov$mwk1_c050 <- marginal_wk(D, T1, .5); cov$mwk2_c050 <- marginal_wk(D, T2, .5)
  cov$mwk1_c075 <- marginal_wk(D, T1, .75); cov$mwk2_c075 <- marginal_wk(D, T2, .75)
  cov$policy_min <- min(POL$age_dependent); cov$policy_max <- max(POL$age_dependent); cov$policy_mean <- mean(POL$age_dependent)
  wcsv(cov, "B6_heart_covariance_decomposition.csv")
  # policy-robust analyses
  wcsv(endpoint_iut(D, T1, T2, nu), "B7_heart_endpoint_IUT.csv")
  env <- envelope(D, T1, T2, nu, V18$box[1], V18$box[2]); rt <- robust_tests(env)
  rt$lower_B_plugin <- env$lower; rt$upper_B_plugin <- env$upper; rt$minR1 <- env$minR1; rt$minR2 <- env$minR2
  wcsv(rt, "B8_heart_dominance_tests.csv")
  band <- multiplier_band(D, T1, T2, nu)
  band$delta_age_free <- sapply(band$c, function(cc) fixed_policy(D, T1, T2, nu0, cc)$delta)
  wcsv(band, "B9_heart_constant_c_curve.csv")
  cg <- seq(0.001, 0.999, by = 0.001); dl <- sapply(cg, function(cc) fixed_policy(D, T1, T2, nu, cc)$delta)
  wcsv(data.frame(delta_min = min(dl), delta_max = max(dl), sign_changes = sum(diff(sign(dl)) != 0)), "B10_heart_constant_c_range.csv")
  # sensitivity analyses (Table S8)
  sens <- list()
  for (k in 1:5) sens[[length(sens) + 1]] <- row_of(fixed_policy(D, T1, T2, crossfit(X, D, T1, T2, df = k, fold = nu$fold), POL$age_dependent), paste("natural spline,", k, "df"))
  sens[[length(sens) + 1]] <- row_of(fixed_policy(D, T1, T2, crossfit(X, D, T1, T2, K = 10L), POL$age_dependent), "K = 10 folds")
  for (sl in c(0, 0.0125, 0.05)) sens[[length(sens) + 1]] <- row_of(fixed_policy(D, T1, T2, nu, age_policy(X, slope = sl)), paste("slope", sl))
  for (ce in c(45, 55)) sens[[length(sens) + 1]] <- row_of(fixed_policy(D, T1, T2, nu, age_policy(X, center = ce)), paste("centre", ce))
  for (c0 in c(0.55, 0.75)) sens[[length(sens) + 1]] <- row_of(fixed_policy(D, T1, T2, nu, age_policy(X, c35 = c0)), paste("c(35) =", c0))
  T1b <- d$T1_flat_down
  sens[[length(sens) + 1]] <- row_of(fixed_policy(D, T1b, T2, crossfit(X, D, T1b, T2), POL$age_dependent), "exercise ECG >= 1 mm with flat/downsloping ST")
  wcsv(do.call(rbind, sens), "B11_heart_sensitivity.csv")
  seeds <- V18$seed + 1000L * (1:50)
  dd <- sapply(seeds, function(sd) fixed_policy(D, T1, T2, crossfit(X, D, T1, T2, seed = sd), POL$age_dependent)$delta)
  wcsv(data.frame(seed = seeds, delta = dd), "B12_heart_fold_seed_sensitivity.csv")
  bx <- do.call(rbind, lapply(box_list, function(b) box_row(D, T1, T2, nu, b[1], b[2], sprintf("[%.2f,%.2f]", b[1], b[2]))))
  ca <- POL$age_dependent; bx <- rbind(bx, box_row(D, T1, T2, nu, pmax(ca - .15, 0), pmin(ca + .15, 1), "c(x) +/- 0.15"))
  wcsv(bx, "B13_heart_box_sensitivity.csv")
  say("PART B done; age-dependent delta", sprintf("%.6f (SE %.6f)", est$delta[3], est$se[3]))
}

################################################################################
# PART C -- simulation calibrated to the coronary data
################################################################################
if ("sim" %in% V18$parts) {
  say("PART C: calibrated simulation,", V18$sim_reps, "replicates for n =", paste(V18$sim_ns, collapse = ", "))
  d0 <- utils::read.csv("data/heart_cleveland_301_age_cad_exercise_ecg_thallium.csv")
  X0 <- d0$X; D0 <- d0$D; T10 <- d0$T1; T20 <- d0$T2; mx <- mean(X0); sx <- stats::sd(X0)
  ns0 <- splines::ns(X0, df = 3); bD <- firth_fit(cbind(1, ns0), D0)
  pfun <- function(x) expit(drop(cbind(1, stats::predict(ns0, x)) %*% bD))
  XT <- cbind(1, D0, (X0 - mx) / sx, D0 * (X0 - mx) / sx)
  b1 <- firth_fit(XT, T10); b2 <- firth_fit(XT, T20)
  rate <- function(b, x, dd) { zz <- (x - mx) / sx; expit(b[1] + b[2] * dd + b[3] * zz + b[4] * dd * zz) }
  r1o <- rate(b1, X0, D0); r2o <- rate(b2, X0, D0)
  RHO <- stats::uniroot(function(rho) sum(rho * pmin(r1o, r2o) + (1 - rho) * r1o * r2o) - sum(T10 * T20), c(0, 1), tol = 1e-12)$root
  truth <- function(x) {
    p <- pfun(x); s1 <- rate(b1, x, 1); s2 <- rate(b2, x, 1); f1 <- rate(b1, x, 0); f2 <- rate(b2, x, 0)
    q1 <- p * s1 + (1 - p) * f1; q2 <- p * s2 + (1 - p) * f2
    C1 <- mean(p * (1 - p) * (s1 - f1)); C2 <- mean(p * (1 - p) * (s2 - f2))
    Cm1 <- mean(p * s1) - mean(p) * mean(q1); Cm2 <- mean(p * s2) - mean(p) * mean(q2)
    out <- do.call(rbind, lapply(list(equal = rep(.5, length(x)), fn_075 = rep(.75, length(x)), age_dependent = age_policy(x)), function(cc) {
      A1 <- mean(cc * p * (1 - s1) + (1 - cc) * (1 - p) * f1); A2 <- mean(cc * p * (1 - s2) + (1 - cc) * (1 - p) * f2)
      data.frame(k1 = C1 / (A1 + C1), k2 = C2 / (A2 + C2), delta = C1 / (A1 + C1) - C2 / (A2 + C2),
                 k1_marg = Cm1 / (A1 + Cm1), k2_marg = Cm2 / (A2 + Cm2), delta_marg = Cm1 / (A1 + Cm1) - Cm2 / (A2 + Cm2)) }))
    out$policy <- c("equal", "fn_075", "age_dependent")
    V1 <- mean((1 - p) * q1); V2 <- mean((1 - p) * q2); h <- C1 * (p - q2) - C2 * (p - q1); b0 <- C1 * V2 - C2 * V1
    attr(out, "box") <- c(B_L = b0 + .1 * mean(h), B_U = b0 + .9 * mean(h), lower = b0 + mean(pmin(.4 * h, .9 * h)), upper = b0 + mean(pmax(.4 * h, .9 * h)))
    out }
  TR <- truth(X0)
  wcsv(cbind(TR, rho = RHO), "C0_calibrated_truth.csv")
  wcsv(data.frame(coefficient = c(paste0("bD", 1:4), paste0("b1_", 1:4), paste0("b2_", 1:4), names(attr(TR, "box"))),
                  value = c(bD, b1, b2, attr(TR, "box"))), "C0_calibrated_truth_model.csv")
  one_rep <- function(n, r) {
    set.seed(20261019L + 7919L * n + r)
    x <- sample(X0, n, replace = TRUE); p <- pfun(x); dd <- as.numeric(stats::runif(n) < p)
    r1 <- rate(b1, x, dd); r2 <- rate(b2, x, dd); u1 <- stats::runif(n)
    share <- stats::runif(n) < RHO; u2 <- ifelse(share, u1, stats::runif(n))
    t1 <- as.numeric(u1 < r1); t2 <- as.numeric(u2 < r2)
    pol <- list(equal = 0.5, fn_075 = 0.75, age_dependent = age_policy(x))
    cols <- c(as.vector(outer(c("cf_", "cf_se_", "af_", "af_se_", "rb_", "rb_se_"), names(pol), paste0)),
              "cf_k1", "cf_se_k1", "cf_k2", "cf_se_k2", "af_k1", "af_se_k1", "af_k2", "af_se_k2",
              "iut_T2", "iut_T1", "box_T2", "box_T1", "ok_cf", "ok_rb")
    row <- c(list(n = n, rep = r), stats::setNames(as.list(rep(NA_real_, length(cols))), cols))
    res <- tryCatch({
      nu <- crossfit(x, dd, t1, t2, seed = 20261017L + r); nu0 <- crossfit(x, dd, t1, t2, model = "intercept", fold = nu$fold)
      for (k in names(pol)) { a <- fixed_policy(dd, t1, t2, nu, pol[[k]]); b <- fixed_policy(dd, t1, t2, nu0, pol[[k]])
        row[[paste0("cf_", k)]] <- a$delta; row[[paste0("cf_se_", k)]] <- a$se_delta
        row[[paste0("af_", k)]] <- b$delta; row[[paste0("af_se_", k)]] <- b$se_delta
        if (k == "age_dependent") { row$cf_k1 <- a$kappa1; row$cf_se_k1 <- a$se_kappa1; row$cf_k2 <- a$kappa2; row$cf_se_k2 <- a$se_kappa2
          row$af_k1 <- b$kappa1; row$af_se_k1 <- b$se_kappa1; row$af_k2 <- b$kappa2; row$af_se_k2 <- b$se_kappa2 } }
      iu <- endpoint_iut(dd, t1, t2, nu); row$iut_T2 <- as.integer(iu$p_IUT_T2 < V18$alpha / 2); row$iut_T1 <- as.integer(iu$p_IUT_T1 < V18$alpha / 2)
      rt <- robust_tests(envelope(dd, t1, t2, nu, V18$box[1], V18$box[2]))
      row$box_T2 <- as.integer(any(rt$certificate[rt$direction == "T2_over_T1"])); row$box_T1 <- as.integer(any(rt$certificate[rt$direction == "T1_over_T2"]))
      row$ok_cf <- 1L; TRUE }, error = function(e) { row$ok_cf <<- 0L; FALSE })
    rbk <- tryCatch(regression_benchmark(x, dd, t1, t2, pol), error = function(e) NULL)
    row$ok_rb <- as.integer(!is.null(rbk))
    if (!is.null(rbk)) for (k in names(pol)) { row[[paste0("rb_", k)]] <- rbk$delta[rbk$label == k]; row[[paste0("rb_se_", k)]] <- rbk$se[rbk$label == k] }
    as.data.frame(row)
  }
  reps <- list()
  for (n in V18$sim_ns) {
    t0 <- Sys.time()
    reps[[as.character(n)]] <- do.call(rbind, lapply(0:(V18$sim_reps - 1L), function(r) one_rep(n, r)))
    say("  n =", n, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")
    wcsv(do.call(rbind, reps), "C1_calibrated_sim_replicates.csv")
  }
  R <- do.call(rbind, reps); z <- stats::qnorm(.975); summ <- list()
  for (n in V18$sim_ns) { g <- R[R$n == n, ]
    for (k in c("equal", "fn_075", "age_dependent")) { tr <- TR$delta[TR$policy == k]; trm <- TR$delta_marg[TR$policy == k]
      for (e in c("cf", "af", "rb")) { est <- g[[paste0(e, "_", k)]]; se <- g[[paste0(e, "_se_", k)]]; ok <- is.finite(est)
        est <- est[ok]; se <- se[ok]
        summ[[length(summ) + 1]] <- data.frame(n = n, target = paste0("delta:", k), estimator = e, n_ok = sum(ok), truth = tr,
          bias = mean(est) - tr, rmse = sqrt(mean((est - tr)^2)), sd = stats::sd(est), mean_se = mean(se), se_sd = mean(se) / stats::sd(est),
          coverage = mean(abs(est - tr) <= z * se), coverage_marginal = mean(abs(est - trm) <= z * se), rejection = mean(abs(est / se) > z)) } }
    for (kk in c("k1", "k2")) { tr <- TR[[kk]][TR$policy == "age_dependent"]; trm <- TR[[paste0(kk, "_marg")]][TR$policy == "age_dependent"]
      for (e in c("cf", "af")) { est <- g[[paste0(e, "_", kk)]]; se <- g[[paste0(e, "_se_", kk)]]; ok <- is.finite(est); est <- est[ok]; se <- se[ok]
        summ[[length(summ) + 1]] <- data.frame(n = n, target = paste0(kk, ":age_dependent"), estimator = e, n_ok = sum(ok), truth = tr,
          bias = mean(est) - tr, rmse = sqrt(mean((est - tr)^2)), sd = stats::sd(est), mean_se = mean(se), se_sd = mean(se) / stats::sd(est),
          coverage = mean(abs(est - tr) <= z * se), coverage_marginal = mean(abs(est - trm) <= z * se), rejection = NA_real_) } }
    summ[[length(summ) + 1]] <- data.frame(n = n, target = "certification", estimator = "cf", n_ok = sum(g$ok_cf), truth = NA_real_,
      bias = NA_real_, rmse = NA_real_, sd = NA_real_, mean_se = NA_real_, se_sd = NA_real_,
      coverage = mean(g$iut_T2, na.rm = TRUE), coverage_marginal = mean(g$box_T2, na.rm = TRUE),
      rejection = mean(pmax(g$iut_T1, g$box_T1), na.rm = TRUE)) }
  S <- do.call(rbind, summ)
  wcsv(S, "C2_calibrated_sim_summary.csv")
  say("PART C done (certification rows: coverage = IUT, coverage_marginal = box, rejection = wrong direction)")
}

################################################################################
# PART D -- categorical covariate: colorectal-cancer counts (Roldan-Nofuentes and Regad, 2025)
################################################################################
if ("crc" %in% V18$parts) {
  say("PART D: categorical-covariate example")
  cnt <- data.frame(X = rep(c(1, 0), each = 8), D = rep(rep(c(1, 0), each = 4), 2),
                    T1 = rep(c(1, 1, 0, 0), 4), T2 = rep(c(1, 0, 1, 0), 4),
                    count = c(51, 13, 1, 9, 3, 1, 1, 34, 17, 5, 0, 4, 1, 0, 1, 27))
  d <- cnt[rep(seq_len(nrow(cnt)), cnt$count), c("X", "D", "T1", "T2")]; rownames(d) <- NULL
  wcsv(d, "D0_crc_diabetes_168_fit_fobt_biopsy.csv")
  X <- d$X; D <- d$D; T1 <- d$T1; T2 <- d$T2
  nu <- crossfit(X, D, T1, T2, model = "factor"); nu0 <- crossfit(X, D, T1, T2, model = "intercept", fold = nu$fold)
  POL <- list(equal = 0.5, fn_075 = 0.75, diabetes_dependent = ifelse(X == 1, 0.75, 0.65))
  est <- rbind(do.call(rbind, lapply(names(POL), function(k) row_of(fixed_policy(D, T1, T2, nu, POL[[k]]), paste("diabetes-conditional", k)))),
               do.call(rbind, lapply(names(POL), function(k) row_of(fixed_policy(D, T1, T2, nu0, POL[[k]]), paste("covariate-free", k)))))
  wcsv(est, "D1_crc_estimates.csv")
  wcsv(endpoint_iut(D, T1, T2, nu), "D2_crc_endpoint_IUT.csv")
  env <- envelope(D, T1, T2, nu, V18$box[1], V18$box[2]); rt <- robust_tests(env)
  rt$lower_B_plugin <- env$lower; rt$upper_B_plugin <- env$upper; rt$minR1 <- env$minR1; rt$minR2 <- env$minR2
  wcsv(rt, "D3_crc_dominance_tests.csv")
  say("PART D done; diabetes-dependent delta", sprintf("%.4f", est$delta[3]))
}
say("All requested parts finished. Output folder:", normalizePath(V18$out))
