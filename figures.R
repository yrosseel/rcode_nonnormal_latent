## figures.R ----------------------------------------------------------------
## All figures for paper 1, written as PDFs into ../figures/ (or into
## ./figures/ when that directory does not exist).
##   fig_contours.pdf  : bivariate shifted log-normal equidensity contours
##   fig_scoredist.pdf : observed vs predicted raw-score distributions (M0/M5)
##   fig_deviances.pdf : deviance differences M0-M5, both item blocks
##   fig_epil_dens.pdf : epilepsy random-effect densities (normal vs skew)
## The last three need replica_results.rds / epil_results.rds.
## Run: Rscript figures.R
## --------------------------------------------------------------------------
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "bcm.R"))
figdir <- file.path(here, "..", "figures")
if (!dir.exists(figdir)) {               # stand-alone copy of the code
  figdir <- file.path(here, "figures")
  dir.create(figdir, showWarnings = FALSE)
}
dir.create(figdir, showWarnings = FALSE)

## ===== 1. bivariate contours ==============================================
## standardized: a = 0, b = 1 per coordinate; density eq. (mv-density)
biv_dens <- function(th1, th2, l1, l2, rho) {
  zi <- function(th, l) if (abs(l) < 1e-8) th else log(pmax(l*th + 1, 1e-300))/l
  ok <- if (l1 > 0) th1 > -1/l1 else th1 < -1/l1
  ok2 <- if (l2 > 0) th2 > -1/l2 else th2 < -1/l2
  z1 <- zi(th1, l1); z2 <- zi(th2, l2)
  Q  <- (z1^2 - 2*rho*z1*z2 + z2^2) / (1 - rho^2)
  d  <- exp(-0.5*Q - l1*z1 - l2*z2) / (2*pi*sqrt(1 - rho^2))
  d[!(ok & ok2)] <- 0
  d
}
Et_ <- function(l) if (abs(l) < 1e-8) 0 else (exp(l^2/2)-1)/l
mode_biv <- function(l1, l2, rho) {
  zm <- -c(rho*l2 + l1, rho*l1 + l2)  # -R lambda
  c(bcm_t(zm[1], l1), bcm_t(zm[2], l2))
}
## conditional mean curve E(Th1 | th2) (eq. cond-mean), standardized a=0,b=1
cond_mean <- function(th2, l1, l2, rho) {
  u1 <- -1/l1; u2 <- -1/l2
  h  <- exp((1 - rho^2) * l1^2)
  u1 + (1/l1) * sqrt(h) * (l2 * (th2 - u2))^(rho * l1 / l2)
}

panels <- list(list(l1 = 0.4, l2 =  0.2, rho =  0.0),
               list(l1 = 0.4, l2 = -0.2, rho =  0.0),
               list(l1 = 0.5, l2 =  0.5, rho =  0.8),
               list(l1 = 0.5, l2 = -0.5, rho =  0.8),
               list(l1 = 0.5, l2 =  0.5, rho = -0.8),
               list(l1 = 0.5, l2 = -0.5, rho = -0.8))

pdf(file.path(figdir, "fig_contours.pdf"), width = 9, height = 6.4)
op <- par(mfrow = c(2, 3), mar = c(3.2, 3.2, 2.2, 0.8),
          mgp = c(2.0, 0.7, 0), cex.main = 1.0)
for (pn in panels) {
  gx <- seq(-2.6, 3.4, length.out = 301)
  gy <- seq(-2.6, 3.4, length.out = 301)
  D  <- outer(gx, gy, biv_dens, l1 = pn$l1, l2 = pn$l2, rho = pn$rho)
  lev <- max(D) * c(.90, .74, .58, .42, .26, .10)
  contour(gx, gy, D, levels = lev, drawlabels = FALSE,
          xlab = expression(theta[1]), ylab = expression(theta[2]),
          main = bquote(lambda[1] == .(pn$l1) ~ "," ~ lambda[2] == .(pn$l2)
                        ~ "," ~ rho == .(pn$rho)))
  abline(h = 0, v = 0, col = "grey80")
  mo <- mode_biv(pn$l1, pn$l2, pn$rho)
  points(mo[1], mo[2], pch = 16, cex = 1.1)                       # mode
  points(Et_(pn$l1), Et_(pn$l2), pch = 3, cex = 1.2, lwd = 2)     # mean
  if (abs(pn$rho) > 0) {                                          # regression curves
    ## support: theta > u = -1/lambda if lambda > 0, theta < u if lambda < 0
    supp <- function(l) if (l > 0) seq(-1/l + 1e-3, 6, length.out = 400)
                        else       seq(-6, -1/l - 1e-3, length.out = 400)
    t2 <- supp(pn$l2)   # E(Th1 | th2): th2 runs over support of Theta2
    cm <- tryCatch(cond_mean(t2, pn$l1, pn$l2, pn$rho), error = function(e) NULL)
    if (!is.null(cm)) lines(cm, t2, lty = 1, col = "blue")
    t1 <- supp(pn$l1)   # E(Th2 | th1)
    cm2 <- tryCatch(cond_mean(t1, pn$l2, pn$l1, pn$rho), error = function(e) NULL)
    if (!is.null(cm2)) lines(t1, cm2, lty = 2, col = "red")
  }
}
par(op); dev.off()
cat("fig_contours.pdf written\n")

## ===== 2-3. replica-based figures =========================================
rrds <- file.path(here, "replica_results.rds")
if (file.exists(rrds)) {
  res <- readRDS(rrds)

  ## ---- score distributions (estimation block) ----
  r <- res$blocks$est
  ks <- length(r$obs_raw_total) - 1L
  pdf(file.path(figdir, "fig_scoredist.pdf"), width = 9, height = 3.8)
  op <- par(mfrow = c(1, 2), mar = c(3.4, 3.6, 2.2, 0.8), mgp = c(2.2, 0.7, 0))
  for (m in c("M0", "M5")) {
    plot(0:ks, r$obs_raw_total, type = "s", lwd = 1, col = "black",
         xlab = "raw score", ylab = "frequency",
         main = sprintf("%s  (Pearson X2 = %.1f)", m, r$pearson_x2[m]))
    lines(0:ks, r$pred_raw[[m]], col = "red", lwd = 2)
    legend("topleft", bty = "n", lwd = c(1, 2), col = c("black", "red"),
           legend = c("observed", "predicted"))
  }
  par(op); dev.off()
  cat("fig_scoredist.pdf written\n")

  ## ---- deviance plot (both blocks) ----
  pdf(file.path(figdir, "fig_deviances.pdf"), width = 7.6, height = 4.2)
  op <- par(mfrow = c(1, 2), mar = c(2.6, 3.6, 2.2, 0.8), mgp = c(2.2, 0.7, 0))
  ## partial order edges (1-df nested pairs)
  edges <- rbind(c("M0","M1"), c("M0","M2"), c("M2","M3"), c("M2","M4"),
                 c("M1","M4"), c("M3","M5"), c("M4","M5"))
  xpos  <- c(M0 = 1.5, M1 = 2.2, M2 = 1.2, M3 = 1.0, M4 = 2.0, M5 = 1.6)
  for (bl in c("est", "xval")) {
    dv <- res$blocks[[bl]]$deviance_diff
    plot(NA, xlim = c(0.6, 2.6), ylim = c(-30, max(dv) * 1.12),
         xaxt = "n", xlab = "", ylab = "deviance difference from M0",
         main = if (bl == "est") "items 1-60" else "items 61-120")
    for (e in seq_len(nrow(edges)))
      lines(xpos[edges[e, ]], dv[edges[e, ]], col = "red")
    points(xpos[names(dv)], dv, pch = 15, col = "red")
    text(xpos[names(dv)], dv, names(dv), pos = 3, cex = 0.85)
  }
  par(op); dev.off()
  cat("fig_deviances.pdf written\n")
} else cat("replica_results.rds not found -- skipping score/deviance figures\n")

## ===== 4. epilepsy density figure =========================================
erds <- file.path(here, "epil_results.rds")
if (file.exists(erds)) {
  er <- readRDS(erds)
  pdf(file.path(figdir, "fig_epil_dens.pdf"), width = 5.4, height = 4.2)
  op <- par(mar = c(3.4, 3.6, 1.0, 0.8), mgp = c(2.2, 0.7, 0))
  with(er$dens, {
    plot(theta, normal, type = "l", lwd = 2,
         xlab = expression(theta), ylab = "density", xlim = c(-2, 2))
    lines(theta, skew, lwd = 2, lty = 2, col = "red")
  })
  legend("topright", bty = "n", lwd = 2, lty = c(1, 2), col = c("black", "red"),
         legend = c(expression(lambda == 0 ~ "(normal)"),
                    expression(hat(lambda) ~ "(BCM)")))
  par(op); dev.off()
  cat("fig_epil_dens.pdf written\n")
} else cat("epil_results.rds not found -- skipping epil figure\n")

cat("FIGURES DONE\n")
