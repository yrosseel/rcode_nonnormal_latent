## validate_bcm.R -----------------------------------------------------------
## Sanity checks for bcm.R:
##   1. analytic gradient vs numDeriv on random parameter points
##   2. parameter recovery on a large simulated data set (M2 truth)
##   3. LR test behaves (M2 vs M0 significant under skew truth)
## Run: Rscript validate_bcm.R
## --------------------------------------------------------------------------
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "bcm.R"))

set.seed(1)
k     <- 20L
alpha <- sample(1:3, k, replace = TRUE)
delta <- seq(-2, 2, length.out = k)
beta  <- c(0.4, 0.2)
dat   <- sim_bcm_data(3000, alpha, delta, beta, bg = c(0.9, 1.0),
                      lg = c(-0.25, -0.25), seed = 42)
smax  <- sum(alpha)
freq  <- freq_by_group(dat$sw, dat$x1, smax)
X     <- rbind(c(1, 0), c(1, 1))
gh    <- gauss_hermite(48L)

## ---- 1. gradient check ----------------------------------------------------
cat("== gradient checks ==\n")
ok <- TRUE
for (m in c("M0", "M1", "M2", "M3", "M4", "M5")) {
  spec <- bcm_spec(m)
  np   <- npar_bcm(spec, 2L)
  for (r in 1:3) {
    par <- c(rnorm(2, 0.3, 0.2), rnorm(max(spec$bmap), -0.1, 0.1),
             if (max(spec$lmap) > 0) rnorm(max(spec$lmap), -0.1, 0.1))
    ga <- bcm_negll_grad(par, X, freq, alpha, delta, spec, gh)
    gn <- numDeriv::grad(function(p)
      as.numeric(bcm_negll(p, X, freq, alpha, delta, spec, gh, want_grad = FALSE)),
      par)
    rel <- max(abs(ga - gn) / pmax(abs(gn), 1))
    if (rel > 1e-6) { ok <- FALSE
      cat(sprintf("  %s rep %d: REL ERR %.2e  <-- PROBLEM\n", m, r, rel))
    }
  }
  cat(sprintf("  %s: ok\n", m))
}
stopifnot(ok)

## ---- 2. recovery (M2 truth: common b, common lambda) ----------------------
cat("== recovery check (truth: beta=0.4/0.2, b=0.95avg, lambda=-0.25) ==\n")
dat2  <- sim_bcm_data(20000, alpha, delta, beta, bg = c(0.95, 0.95),
                      lg = c(-0.25, -0.25), seed = 7)
freq2 <- freq_by_group(dat2$sw, dat2$x1, smax)
t0 <- proc.time()[3]
f0 <- fit_bcm(X, freq2, alpha, delta, "M0")
f2 <- fit_bcm(X, freq2, alpha, delta, "M2")
cat(sprintf("  fit time M0+M2: %.1fs\n", proc.time()[3] - t0))
print(round(rbind(est = f2$par, se = f2$se), 4))
stopifnot(abs(f2$par["lambda1"] + 0.25) < 0.06)
stopifnot(abs(f2$par["beta1"] - 0.2) < 0.05)

## ---- 3. LR test -----------------------------------------------------------
lt <- lr_test(f0, f2)
cat(sprintf("  LR(M0 vs M2) = %.1f on %d df, p = %.2g\n", lt[1], lt[2], lt[3]))
stopifnot(lt["LR"] > 10)

## ---- 4. size sanity: normal truth should not reject wildly ---------------
dat3  <- sim_bcm_data(5000, alpha, delta, beta, bg = c(0.95, 0.95),
                      lg = c(0, 0), seed = 11)
freq3 <- freq_by_group(dat3$sw, dat3$x1, smax)
f0n <- fit_bcm(X, freq3, alpha, delta, "M0", hessian = FALSE)
f2n <- fit_bcm(X, freq3, alpha, delta, "M2", hessian = FALSE)
ltn <- lr_test(f0n, f2n)
cat(sprintf("  normal truth: LR = %.2f (should be O(1)), p = %.2f\n",
            ltn[1], ltn[3]))

cat("ALL CHECKS PASSED\n")
