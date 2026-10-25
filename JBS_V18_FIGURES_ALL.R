################################################################################
# JBS_V18_FIGURES_ALL.R
#
# Produces Figures 1-5 of the manuscript
#   "Covariate-Dependent Loss-Weighted Kappa for Comparing Paired Binary
#    Diagnostic Tests: Orthogonal Estimation and Policy-Robust Dominance"
# from the archived results of the R programs, using base R graphics only
# (no additional packages).
#
# Author : Ismail Radi Ismail Fkhadah, Department of Biostatistics, School of
#          Medicine, University of Granada, Spain (ismailfk2021@correo.ugr.es)
# License: MIT
#
# HOW TO RUN (RStudio). Set the working directory to the repository root, then
#   source("R/JBS_V18_FIGURES_ALL.R", encoding = "UTF-8")
# Optional environment variables:
#   JBS_RESULTS  folder with the archived results (default "results")
#   JBS_FIG_OUT  output folder (default "figures")
#
# For each figure three files are written:
#   FigureK_<name>.tiff  600 dpi, LZW-compressed, RGB  (journal submission file)
#   FigureK_<name>.pdf   vector graphics with embedded fonts (cairo_pdf)
#   FigureK_<name>.png   300 dpi (the image embedded in the manuscript file)
# Sizes (inches): Figure 1 6.6 x 3.0; Figure 2 6.6 x 4.4; Figure 3 5.4 x 3.3;
# Figure 4 6.2 x 3.4; Figure 5 6.8 x 3.1. All figures are within the
# 17.4-cm (6.85-inch) full-page width.
#
# Inputs (all written by the R programs of this repository):
#   Figure 1  results/JBS_V13_3_RESULTS/MONTE_CARLO_FULL/CHECKPOINT_MC_SUMMARY.csv
#   Figure 2  results/JBS_V14_BOUNDARY_RESULTS/V14_RAW_*.csv and
#             results/JBS_V14_CONFIRM_1000_N1500_RHO025/V14_RAW_*.csv
#   Figure 3  results/JBS_V15_NUISANCE_RESULTS/V15_NUISANCE_RATES_RAW.csv
#   Figure 4  results/JBS_V13_3_RESULTS/JBS_Tables/V11_Table11_Integrated_Contrast_Multiplier_Grid_Band.csv
#   Figure 5  results/V18_R_RESULTS/B9_heart_constant_c_curve.csv and
#             results/V18_R_RESULTS/C2_calibrated_sim_summary.csv
#             (written by R/JBS_V18_SECOND_APPLICATION_AND_SENSITIVITY.R)
################################################################################

RES <- Sys.getenv("JBS_RESULTS", "results")
OUT <- Sys.getenv("JBS_FIG_OUT", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(...) utils::read.csv(file.path(RES, ...), stringsAsFactors = FALSE)

# Colours and graphical style (identical in all figures)
COL <- c(blue = "#2a78d6", orange = "#eb6834", green = "#1baf7a", violet = "#6250d6",
         ink = "#52514e", grid = "#e4e3df", band = "#2a78d624", mcse = "#e4e3dfB3")
style <- function() par(family = "sans", ps = 9, cex = 1, mgp = c(1.9, 0.5, 0), tcl = -0.3, bty = "n",
                        las = 1, col.axis = COL[["ink"]], col.lab = "black", fg = COL[["ink"]], lend = 1)
# Numeric axis labels with a true minus sign (plotmath), identical in TIFF, PNG and PDF output
MINUS <- if (has_cairo_early <- isTRUE(capabilities("cairo"))) "\u2212" else "-"
num_lab <- function(v, d) sub("^-", MINUS, sprintf(paste0("%.", d, "f"), v))
frame_axes <- function(xat = NULL, xlab = NULL, yat = NULL, ylab = NULL, xlog = FALSE, ylog = FALSE, ygrid = TRUE, xgrid = TRUE) {
  if (xgrid && !is.null(xat)) abline(v = xat, col = COL[["grid"]], lwd = 0.6)
  if (ygrid && !is.null(yat)) abline(h = yat, col = COL[["grid"]], lwd = 0.6)
  box(bty = "l", col = COL[["ink"]], lwd = 0.8)
}
panel_title <- function(txt) mtext(txt, side = 3, line = 0.35, adj = 0, cex = 1)
# Cairo devices embed the fonts in the PDF and take Greek letters and the minus sign
# from the same sans-serif font as the text (cairoSymbolFont, R >= 4.0).
has_cairo <- isTRUE(capabilities("cairo"))
sym <- if (has_cairo && exists("cairoSymbolFont", envir = asNamespace("grDevices")))
  grDevices::cairoSymbolFont("sans", usePUA = FALSE) else NULL
save_fig <- function(name, w, h, plotfun) {
  tf <- file.path(OUT, paste0(name, ".tiff"))
  if (has_cairo) {
    args <- list(units = "in", bg = "white", type = "cairo"); if (!is.null(sym)) args$symbolfamily <- sym
    do.call(grDevices::tiff, c(list(tf, width = w, height = h, res = 600, compression = "lzw"), args))
    style(); plotfun(); grDevices::dev.off()
    do.call(grDevices::png, c(list(file.path(OUT, paste0(name, ".png")), width = w, height = h, res = 300), args))
    style(); plotfun(); grDevices::dev.off()
    pargs <- list(file.path(OUT, paste0(name, ".pdf")), width = w, height = h, family = "sans"); if (!is.null(sym)) pargs$symbolfamily <- sym
    do.call(grDevices::cairo_pdf, pargs)
    style(); plotfun(); grDevices::dev.off()
  } else {
    grDevices::tiff(tf, width = w, height = h, units = "in", res = 600, compression = "lzw", bg = "white")
    style(); plotfun(); grDevices::dev.off()
    grDevices::png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "in", res = 300, bg = "white")
    style(); plotfun(); grDevices::dev.off()
    grDevices::pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h, family = "Helvetica", useDingbats = FALSE)
    style(); plotfun(); grDevices::dev.off()
  }
  cat("written:", name, "(.tiff, .pdf, .png)\n")
}
lab_box <- function(x, y, lab, cex = 0.85) {   # text on a translucent white box
  w <- strwidth(lab, cex = cex); h <- strheight(lab, cex = cex)
  rect(x - 0.004, y - 0.75 * h, x + w + 0.004, y + 0.75 * h, col = "#FFFFFFD9", border = NA)
  text(x, y, lab, adj = c(0, 0.5), cex = cex)
}
cp_ci <- function(k, m) c(if (k > 0) stats::qbeta(.025, k, m - k + 1) else 0, if (k < m) stats::qbeta(.975, k + 1, m - k) else 1)

########################################################################## Figure 1
mc <- rd("JBS_V13_3_RESULTS", "MONTE_CARLO_FULL", "CHECKPOINT_MC_SUMMARY.csv")
mc <- mc[mc$method == "crossfit", ]
ns <- c(250, 500, 858, 1500)
fig1 <- function() {
  layout(matrix(1:2, 1)); par(mar = c(3.2, 3.6, 1.6, 0.6), oma = c(0, 0, 0, 0))
  sc <- c("null", "t2_local", "t1_strong", "t2_strong", "tradeoff")
  lab <- c("Null", "Weak T2", "Strong T1", "Strong T2", "Trade-off")
  cl <- COL[c("blue", "orange", "green", "violet", "ink")]; pc <- c(16, 15, 17, 18, 25); lt <- c(1, 2, 4, 3, 1)
  plot(NA, xlim = c(0.7, 4.3), ylim = c(0.903, 0.975), xaxt = "n", yaxt = "n", xlab = expression(italic(n)), ylab = "95% CI coverage")
  mcse <- 1.96 * sqrt(.95 * .05 / 1000)
  rect(0, .95 - mcse, 5, .95 + mcse, col = COL[["mcse"]], border = NA)
  frame_axes(1:4, yat = seq(0.92, 0.97, 0.01))
  abline(h = 0.95, col = COL[["ink"]], lwd = 0.8)
  axis(1, at = 1:4, labels = ns); axis(2, at = seq(0.92, 0.97, 0.01), labels = sprintf("%.2f", seq(0.92, 0.97, 0.01)))
  for (k in seq_along(sc)) {
    g <- mc[mc$scenario == sc[k], ]
    m <- sapply(ns, function(v) mean(g$coverage[g$n == v])); lo <- sapply(ns, function(v) min(g$coverage[g$n == v]))
    hi <- sapply(ns, function(v) max(g$coverage[g$n == v])); x <- 1:4 + (k - 3) * 0.07
    segments(x, lo, x, hi, col = cl[k], lwd = 0.8)
    lines(x, m, col = cl[k], lty = lt[k], lwd = 1.2); points(x, m, pch = pc[k], col = cl[k], bg = cl[k], cex = 0.85)
  }
  legend("bottom", legend = lab, col = cl, pch = pc, pt.bg = cl, lty = lt, lwd = 1.2, ncol = 3, bty = "n", cex = 0.78, seg.len = 2, x.intersp = 0.6, inset = c(0, 0.005),
         text.width = 1.12 * max(strwidth(lab, cex = 0.78)))
  panel_title("(a) Coverage, mean and range over 6 cells")
  plot(NA, xlim = c(0.7, 4.3), ylim = c(0, 0.82), xaxt = "n", yaxt = "n", xlab = expression(italic(n)), ylab = "Two-sided rejection rate")
  frame_axes(1:4, yat = seq(0, 0.8, 0.2)); abline(h = 0.05, col = COL[["ink"]], lwd = 0.8)
  axis(1, at = 1:4, labels = ns); axis(2, at = seq(0, 0.8, 0.2), labels = sprintf("%.1f", seq(0, 0.8, 0.2)))
  g <- mc[mc$scenario == "null", ]
  m <- sapply(ns, function(v) mean(g$rejection_rate[g$n == v]))
  segments(1:4, sapply(ns, function(v) min(g$rejection_rate[g$n == v])), 1:4, sapply(ns, function(v) max(g$rejection_rate[g$n == v])), col = COL[["blue"]], lwd = 0.8)
  lines(1:4, m, col = COL[["blue"]], lwd = 1.4); points(1:4, m, pch = 16, col = COL[["blue"]], cex = 0.85)
  rh <- c(0.25, 0.55, 0.85); pr <- c(15, 17, 18); lr <- c(2, 4, 3)
  for (i in 1:3) { g <- mc[mc$scenario == "t2_local" & abs(mc$dependence - rh[i]) < 1e-8, ]
    m <- sapply(ns, function(v) mean(g$rejection_rate[g$n == v]))
    lines(1:4, m, col = COL[["orange"]], lty = lr[i], lwd = 1.2); points(1:4, m, pch = pr[i], col = COL[["orange"]], cex = 0.85) }
  legend("topleft", legend = c("Null: type I error", expression("Weak T2 power, " * rho == 0.25), expression("Weak T2 power, " * rho == 0.55), expression("Weak T2 power, " * rho == 0.85)),
         col = COL[c("blue", "orange", "orange", "orange")], pch = c(16, pr), lty = c(1, lr), lwd = 1.2, bty = "n", cex = 0.8, seg.len = 2.2)
  panel_title("(b) Type I error and power")
}
save_fig("Figure1_simulation_coverage_rejection", 6.6, 3.0, fig1)

########################################################################## Figure 2
files <- list.files(file.path(RES, "JBS_V14_BOUNDARY_RESULTS"), pattern = "^V14_RAW_.*\\.csv$", full.names = TRUE)
cells <- do.call(rbind, lapply(files, function(f) {
  r <- utils::read.csv(f); m <- regmatches(basename(f), regexec("RAW_(lower_zero|upper_zero)_n([0-9]+)_rho0([0-9]+)_(oracle|crossfit)", basename(f)))[[1]]
  anyc <- tapply(r$T1_certificate + r$T2_certificate, r$rep, max) > 0
  data.frame(direction = m[2], n = as.integer(m[3]), rho = as.numeric(paste0("0.", m[4])), mode = m[5], k = sum(anyc), reps = length(anyc))
}))
cells <- cells[order(cells$direction, cells$n, cells$rho), ]
cf <- cells[cells$mode == "crossfit", ]
cf$lab <- paste0(ifelse(cf$direction == "lower_zero", "L", "U"), " n=", cf$n, " ρ=", cf$rho)
conf <- sapply(c("lower_zero", "upper_zero"), function(dname) {
  r <- utils::read.csv(file.path(RES, "JBS_V14_CONFIRM_1000_N1500_RHO025", paste0("V14_RAW_", dname, "_n1500_rho025_crossfit.csv")))
  sum(tapply(r$T1_certificate + r$T2_certificate, r$rep, max) > 0) })
fig2 <- function() {
  par(mar = c(5.2, 7.2, 0.8, 0.8))
  ny <- nrow(cf)
  plot(NA, xlim = c(-0.2, 8.2), ylim = c(ny + 0.7, 0.3), yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  frame_axes(0:8, yat = NULL, ygrid = FALSE)
  axis(1, at = 0:8); mtext("Familywise false-dominance rate (%) with exact 95% CI", side = 1, line = 1.8)
  labs <- sapply(seq_len(ny), function(i) as.expression(bquote(.(ifelse(cf$direction[i] == "lower_zero", "L", "U")) ~ "n=" * .(cf$n[i]) ~ rho * "=" * .(cf$rho[i]))))
  axis(2, at = 1:ny, labels = labs, cex.axis = 0.75, las = 1)
  for (md in c("oracle", "crossfit")) {
    d <- cells[cells$mode == md, ]; d <- d[order(d$direction, d$n, d$rho), ]
    off <- if (md == "oracle") -0.18 else 0.18; cl <- if (md == "oracle") COL[["violet"]] else COL[["orange"]]
    ci <- t(mapply(cp_ci, d$k, d$reps)) * 100
    segments(ci[, 1], 1:ny + off, ci[, 2], 1:ny + off, col = cl, lwd = 0.8)
    points(d$k / d$reps * 100, 1:ny + off, pch = if (md == "oracle") 16 else 15, col = cl, cex = 0.8)
  }
  for (i in 1:2) { dname <- c("lower_zero", "upper_zero")[i]
    yy <- which(cf$direction == dname & cf$n == 1500 & abs(cf$rho - 0.25) < 1e-8) + 0.45
    ci <- cp_ci(conf[i], 1000) * 100
    segments(ci[1], yy, ci[2], yy, col = COL[["green"]], lwd = 1.2); points(conf[i] / 10, yy, pch = 18, col = COL[["green"]], cex = 1.1) }
  abline(v = 5, col = COL[["ink"]], lwd = 0.9); text(5.1, 0.1, "nominal 5%", adj = c(0, 0.5), cex = 0.83, col = COL[["ink"]], xpd = TRUE)
  ytop <- -(2.9 * par("csi")) / par("pin")[2]          # 2.9 text lines below the plot region
  par(usr = c(0, 1, 0, 1), xpd = NA)
  labs2 <- c("Oracle nuisances (200 reps)", "Cross-fitted nuisances (200 reps)", "Cross-fitted confirmation (1,000 reps)")
  legend(x = grconvertX(0.5, "ndc", "user"), y = ytop, xjust = 0.5, yjust = 1, horiz = TRUE, legend = labs2,
         text.width = strwidth(labs2, cex = 0.76) + 0.03,
         col = COL[c("violet", "orange", "green")], pch = c(16, 15, 18), lty = 1, lwd = c(0.8, 0.8, 1.2), bty = "n", cex = 0.76, seg.len = 1.1, x.intersp = 0.45)
  par(xpd = FALSE)
}
save_fig("Figure2_boundary_familywise_error", 6.6, 4.4, fig2)

########################################################################## Figure 3
nr <- rd("JBS_V15_NUISANCE_RESULTS", "V15_NUISANCE_RATES_RAW.csv")
nr <- nr[nr$success %in% c(TRUE, "TRUE"), ]
fig3 <- function() {
  par(mar = c(3.3, 4.0, 0.8, 0.8))
  plot(NA, xlim = c(230, 3250), ylim = c(0.009, 0.13), log = "xy", xaxt = "n", yaxt = "n",
       xlab = expression("Sample size " * italic(n) * " (training size 4" * italic(n) * "/5, rounded down)"), ylab = "")
  mtext(expression("Mean " * italic(Q)[italic(n)] == sqrt(italic(n)) ~ italic(e)[italic(p)] ~ max[italic(j)] ~ italic(e)[italic(q)[italic(j)]]), side = 2, line = 2.6, las = 0)
  xt <- c(250, 500, 858, 1500, 3000); yt <- c(0.01, 0.02, 0.05, 0.10)
  frame_axes(xt, yat = yt); axis(1, at = xt, labels = xt); axis(2, at = yt, labels = sprintf("%.2f", yt))
  mods <- c("linear", "spline_df3", "spline_df6"); cl <- COL[c("blue", "orange", "green")]; pc <- c(16, 15, 17)
  for (i in 1:3) { g <- nr[nr$model == mods[i], ]; nn <- sort(unique(g$n))
    m <- sapply(nn, function(v) mean(g$sqrt_n_max_prod[g$n == v]))
    lines(nn, m, col = cl[i], lwd = 1.6); points(nn, m, pch = pc[i], col = cl[i], cex = 0.9) }
  lines(c(250, 3000), 0.064 * (c(250, 3000) / 250)^-0.5, col = COL[["ink"]], lty = 3, lwd = 0.9)
  legend("topright", legend = c("Logistic, linear", "Natural spline, 3 df (primary)", "Natural spline, 6 df", paste0("Reference slope ", MINUS, "1/2")),
         col = c(cl, COL[["ink"]]), pch = c(pc, NA), lty = c(1, 1, 1, 3), lwd = c(1.6, 1.6, 1.6, 0.9), bty = "n", cex = 0.83, seg.len = 2.2)
}
save_fig("Figure3_nuisance_rates", 5.4, 3.3, fig3)

########################################################################## Figure 4
b4 <- rd("JBS_V13_3_RESULTS", "JBS_Tables", "V11_Table11_Integrated_Contrast_Multiplier_Grid_Band.csv")
fig4 <- function() {
  par(mar = c(3.3, 3.9, 0.8, 0.8))
  yl <- range(c(b4$simultaneous_grid_low, b4$simultaneous_grid_high)) + c(-0.01, 0.01)
  plot(NA, xlim = c(0.03, 0.97), ylim = yl, xlab = expression("Constant false-negative loss share " * italic(c)), ylab = "", yaxt = "n")
  mtext(bquote(Delta(italic(c)) * .(paste0(", Hinselmann ", MINUS, " Schiller"))), side = 2, line = 2.7, las = 0)
  yt <- pretty(yl); frame_axes(pretty(c(0, 1)), yat = yt); axis(2, at = yt, labels = num_lab(yt, 1))
  polygon(c(b4$c, rev(b4$c)), c(b4$simultaneous_grid_low, rev(b4$simultaneous_grid_high)), col = COL[["band"]], border = NA)
  lines(b4$c, b4$pointwise_low, col = COL[["blue"]], lty = 2, lwd = 1); lines(b4$c, b4$pointwise_high, col = COL[["blue"]], lty = 2, lwd = 1)
  lines(b4$c, b4$delta, col = COL[["blue"]], lwd = 2); abline(h = 0, col = COL[["ink"]], lwd = 0.9)
  for (c0 in c(0.50, 0.75)) { dv <- b4$delta[abs(b4$c - c0) < 1e-9]
    points(c0, dv, pch = 21, bg = COL[["orange"]], col = "white", cex = 1.3, lwd = 1.2)
    lab_box(c0 + 0.018, dv + 0.072, bquote(italic(c) == .(sprintf("%.2f", c0))))
    lab_box(c0 + 0.018, dv + 0.037, bquote(Delta == .(sub("^-", MINUS, sprintf("%.3f", dv))))) }
  legend("bottomleft", legend = c("95% simultaneous finite-grid band", "95% pointwise interval", "Cross-fitted estimate"),
         col = c(COL[["band"]], COL[["blue"]], COL[["blue"]]), lty = c(1, 2, 1), lwd = c(9, 1, 2), bty = "n", cex = 0.85, seg.len = 2.2)
}
save_fig("Figure4_application_constant_policies", 6.2, 3.4, fig4)

########################################################################## Figure 5
b5 <- rd("V18_R_RESULTS", "B9_heart_constant_c_curve.csv")
s5 <- rd("V18_R_RESULTS", "C2_calibrated_sim_summary.csv")
fig5 <- function() {
  layout(matrix(1:2, 1)); par(mar = c(3.3, 3.9, 1.6, 0.6))
  plot(NA, xlim = c(0.03, 0.97), ylim = c(-0.50, 0.06), xlab = expression("Constant false-negative loss share " * italic(c)), ylab = "", yaxt = "n")
  mtext(bquote(Delta(italic(c)) * .(paste0(", exercise ECG ", MINUS, " thallium"))), side = 2, line = 2.7, las = 0)
  yt <- seq(-0.5, 0, 0.1); frame_axes(pretty(c(0, 1)), yat = yt); axis(2, at = yt, labels = num_lab(yt, 1))
  polygon(c(b5$c, rev(b5$c)), c(b5$simultaneous_low, rev(b5$simultaneous_high)), col = COL[["band"]], border = NA)
  lines(b5$c, b5$pointwise_low, col = COL[["blue"]], lty = 2, lwd = 1); lines(b5$c, b5$pointwise_high, col = COL[["blue"]], lty = 2, lwd = 1)
  lines(b5$c, b5$delta, col = COL[["blue"]], lwd = 2); lines(b5$c, b5$delta_age_free, col = COL[["orange"]], lty = 4, lwd = 1.6)
  abline(h = 0, col = COL[["ink"]], lwd = 0.9)
  legend("bottomleft", legend = c("95% simultaneous band", "95% pointwise interval", "Proposed (age-conditional)", "Age-free"),
         col = c(COL[["band"]], COL[["blue"]], COL[["blue"]], COL[["orange"]]), lty = c(1, 2, 1, 4), lwd = c(8, 1, 2, 1.6),
         ncol = 1, bty = "n", cex = 0.74, seg.len = 2, x.intersp = 0.7, y.intersp = 0.95)
  panel_title(expression("(a) Data, " * italic(n) == 301))
  nsim <- sort(unique(s5$n)); xi <- seq_along(nsim)
  plot(NA, xlim = c(0.7, length(nsim) + 0.3), ylim = c(45, 100), xaxt = "n", yaxt = "n", xlab = expression(italic(n)), ylab = "Coverage of 95% CI (%)")
  mcse <- 1.96 * sqrt(.95 * .05 / 1000) * 100
  rect(0, 95 - mcse, 10, 95 + mcse, col = COL[["mcse"]], border = NA)
  frame_axes(xi, yat = seq(50, 100, 10)); abline(h = 95, col = COL[["ink"]], lwd = 0.8)
  axis(1, at = xi, labels = nsim); axis(2, at = seq(50, 100, 10))
  tg <- c("k1:age_dependent", "k2:age_dependent", "delta:age_dependent"); pc <- c(16, 15, 17)
  for (i in 1:3) for (e in c("cf", "af")) {
    g <- s5[s5$target == tg[i] & s5$estimator == e, ]; g <- g[order(g$n), ]
    cl <- if (e == "cf") COL[["blue"]] else COL[["orange"]]
    lines(xi, 100 * g$coverage, col = cl, lty = if (e == "cf") 1 else 2, lwd = 1.3); points(xi, 100 * g$coverage, pch = pc[i], col = cl, cex = 0.9) }
  legend("bottomleft", legend = c(expression(kappa[1] * " (exercise ECG), proposed"), expression(kappa[1] * " (exercise ECG), age-free"),
                                  expression(kappa[2] * " (thallium), proposed"), expression(kappa[2] * " (thallium), age-free"),
                                  expression(Delta * ", proposed"), expression(Delta * ", age-free")),
         col = rep(COL[c("blue", "orange")], 3), pch = rep(pc, each = 2), lty = rep(c(1, 2), 3), lwd = 1.3, bty = "n", cex = 0.72, seg.len = 2.2)
  panel_title("(b) Calibrated simulation")
}
save_fig("Figure5_coronary_application_and_calibrated_simulation", 6.8, 3.1, fig5)
cat("All figures written to", normalizePath(OUT), "\n")
