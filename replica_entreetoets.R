## replica_entreetoets.R ----------------------------------------------------
## Synthetic replica of the Entreetoets latent regression design (paper §6.1).
## 120 items (blocks 1-60 = estimation, 61-120 = cross-validation),
## n = 16,000, gender covariate, left-skew latent distribution.
## Fits M0-M5 on both item blocks; saves everything needed for the paper
## (deviance table, parameter tables, latent moments, predicted raw-score
## distributions, Pearson X2, parameter correlation matrix).
## Run: Rscript replica_entreetoets.R      (writes replica_results.rds)
## --------------------------------------------------------------------------
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "bcm.R"))

set.seed(2009)

## ---- design ---------------------------------------------------------------
n      <- 16000L
kblock <- 60L
alpha  <- sample(1:4, 2 * kblock, replace = TRUE, prob = c(.2, .35, .3, .15))
delta  <- c(seq(-2.5, 2.5, length.out = kblock)[sample(kblock)],
            seq(-2.5, 2.5, length.out = kblock)[sample(kblock)])
## truth: median regression a = 1.0 + 0.10 * girl; group scales/skews
beta_t <- c(1.00, 0.10)
b_t    <- c(1.05, 0.95)     # boys, girls
lam_t  <- c(-0.30, -0.20)   # boys, girls  (skew to the left)

dat <- sim_bcm_data(n, alpha, delta, beta_t, b_t, lam_t, pgirl = 0.493,
                    seed = 20090115)
X   <- rbind(c(1, 0), c(1, 1))
wg  <- c(mean(dat$x1 == 0), mean(dat$x1 == 1))

blocks <- list(est = 1:kblock, xval = kblock + 1:kblock)
models <- c("M0", "M1", "M2", "M3", "M4", "M5")

results <- list(truth = list(beta = beta_t, b = b_t, lambda = lam_t,
                             alpha = alpha, delta = delta, wg = wg),
                blocks = list())

for (bl in names(blocks)) {
  items <- blocks[[bl]]
  al    <- alpha[items]; de <- delta[items]
  sw    <- as.numeric(dat$Y[, items] %*% al)
  sr    <- rowSums(dat$Y[, items])
  freq  <- freq_by_group(sw, dat$x1, sum(al))
  fits  <- list()
  for (m in models) {
    cat(sprintf("[%s] fitting %s ...\n", bl, m))
    t0 <- proc.time()[3]
    fits[[m]] <- fit_bcm(X, freq, al, de, model = m, nq = 180L,
                         hessian = (m %in% c("M2", "M5")))
    cat(sprintf("      logLik = %.2f  (%.1fs, converged: %s)\n",
                fits[[m]]$loglik, proc.time()[3] - t0, fits[[m]]$converged))
  }
  ## deviance differences from M0
  dev <- vapply(models, function(m) 2 * (fits[[m]]$loglik - fits$M0$loglik),
                numeric(1))
  ## latent moments per model
  mom <- t(vapply(models, function(m) latent_moments(fits[[m]], wg), numeric(2)))
  ## observed & predicted raw-score distributions (M0 and M5)
  obs_raw <- lapply(0:1, function(g)
    tabulate(sr[dat$x1 == g] + 1L, nbins = length(items) + 1L))
  ng      <- c(sum(dat$x1 == 0), sum(dat$x1 == 1))
  pred    <- lapply(c("M0", "M5"), function(m) {
    P <- predicted_raw(fits[[m]])
    sweep(P, 2, ng, "*")     # expected frequencies per group
  })
  names(pred) <- c("M0", "M5")
  ## total (both groups) observed and expected; Pearson X2 with tail pooling
  obs_tot  <- obs_raw[[1]] + obs_raw[[2]]
  x2 <- vapply(pred, function(E) {
    e_tot <- rowSums(E)
    keep  <- e_tot >= 5
    ## pool the low-expected tail cells (guard: only if any)
    o <- obs_tot[keep]; e <- e_tot[keep]
    if (any(!keep)) { o <- c(sum(obs_tot[!keep]), o); e <- c(sum(e_tot[!keep]), e) }
    sum((o - e)^2 / e)
  }, numeric(1))
  results$blocks[[bl]] <- list(
    fits = lapply(fits, function(f) f[c("par", "se", "loglik", "npar",
                                        "model", "converged", "pcor")]),
    deviance_diff = dev, latent_moments = mom,
    obs_raw_total = obs_tot, pred_raw = lapply(pred, rowSums),
    pearson_x2 = x2, ng = ng)
}

saveRDS(results, file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "replica_results.rds"))

## ---- console summary ------------------------------------------------------
for (bl in names(results$blocks)) {
  r <- results$blocks[[bl]]
  cat(sprintf("\n===== block %s =====\n", bl))
  cat("deviance differences from M0:\n")
  print(round(r$deviance_diff, 1))
  cat("latent moments:\n"); print(round(r$latent_moments, 3))
  cat(sprintf("Pearson X2: M0 = %.1f, M5 = %.1f\n",
              r$pearson_x2["M0"], r$pearson_x2["M5"]))
  cat("M5 estimates:\n")
  print(round(rbind(est = r$fits$M5$par, se = r$fits$M5$se), 4))
  cat("M2 max |param corr|:",
      round(max(abs(r$fits$M2$pcor[upper.tri(r$fits$M2$pcor)])), 3), "\n")
}
cat("\nDONE\n")
