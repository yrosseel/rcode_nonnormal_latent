## sim_S1S2.R ---------------------------------------------------------------
## Simulation studies S1 (size/power of the LR test, recovery of lambda)
## and S2 (robustness of structural estimates under misspecification).
## Design (paper §7): k = 40 items, restricted 2PLM with known item
## parameters; binary covariate, beta = (0.5, 0.25); Var(Theta) = 1;
## common b and lambda across groups (M2-type truth).
## Cells: lambda in {0, ±0.10, ±0.25} x n in {500, 2000, 8000}.
## Per replicate: fit M0 and M2 (with SEs), record everything.
## Env: REPS (default 500), NQ (default 96), CORES (default 18).
## Run: Rscript sim_S1S2.R                (writes sim_S1S2.rds)
## --------------------------------------------------------------------------
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "bcm.R"))

REPS  <- as.integer(Sys.getenv("REPS",  "500"))
NQ    <- as.integer(Sys.getenv("NQ",    "96"))
CORES <- as.integer(Sys.getenv("CORES", "18"))

## fixed design
set.seed(4040)
k     <- 40L
alpha <- sample(1:4, k, replace = TRUE, prob = c(.2, .35, .3, .15))
delta <- seq(-2, 2, length.out = k)[sample(k)]
beta_t <- c(0.5, 0.25)
X     <- rbind(c(1, 0), c(1, 1))

cells <- expand.grid(lambda = c(0, -0.10, 0.10, -0.25, 0.25),
                     n = c(500L, 2000L, 8000L))

one_rep <- function(cell, rep) {
  lam <- cells$lambda[cell]; n <- cells$n[cell]
  b   <- 1 / sqrt(Vart_bcm(lam))          # Var(Theta) = 1
  dat <- sim_bcm_data(n, alpha, delta, beta_t, bg = c(b, b), lg = c(lam, lam),
                      pgirl = 0.5, seed = cell * 100000L + rep)
  freq <- freq_by_group(dat$sw, dat$x1, sum(alpha))
  t0  <- proc.time()[3]
  f0  <- tryCatch(fit_bcm(X, freq, alpha, delta, "M0", nq = NQ),
                  error = function(e) NULL)
  f2  <- tryCatch(fit_bcm(X, freq, alpha, delta, "M2", nq = NQ),
                  error = function(e) NULL)
  if (is.null(f0) || is.null(f2))
    return(data.frame(cell = cell, rep = rep, lambda_t = lam, n = n,
                      ok = FALSE))
  lt  <- lr_test(f0, f2)
  data.frame(cell = cell, rep = rep, lambda_t = lam, n = n, ok = TRUE,
             conv0 = f0$converged, conv2 = f2$converged,
             LR = lt["LR"], reject = (lt["p"] < 0.05),
             lam_hat = f2$par["lambda1"], lam_se = f2$se[4],
             b1_0 = f0$par["beta1"], b1se_0 = f0$se[2],
             b1_2 = f2$par["beta1"], b1se_2 = f2$se[2],
             mean_0 = latent_moments(f0, c(.5, .5))["mean"],
             var_0  = latent_moments(f0, c(.5, .5))["variance"],
             mean_2 = latent_moments(f2, c(.5, .5))["mean"],
             var_2  = latent_moments(f2, c(.5, .5))["variance"],
             maxcor2 = if (!is.null(f2$pcor))
                         max(abs(f2$pcor[upper.tri(f2$pcor)])) else NA_real_,
             secs = proc.time()[3] - t0)
}

grid <- expand.grid(rep = seq_len(REPS), cell = seq_len(nrow(cells)))
cat(sprintf("running %d replicates on %d cores (REPS=%d, NQ=%d)\n",
            nrow(grid), CORES, REPS, NQ))
t0  <- proc.time()[3]
res <- parallel::mclapply(seq_len(nrow(grid)), function(i)
         tryCatch(one_rep(grid$cell[i], grid$rep[i]),
                  error = function(e) NULL),
         mc.cores = CORES, mc.preschedule = TRUE)
res <- do.call(rbind, Filter(Negate(is.null), res))
cat(sprintf("done in %.1f min; %d/%d rows ok\n",
            (proc.time()[3] - t0) / 60, sum(res$ok), nrow(grid)))
saveRDS(list(design = list(alpha = alpha, delta = delta, beta = beta_t),
             cells = cells, res = res),
        file.path(here, "sim_S1S2.rds"))

## ---- aggregate summary ----------------------------------------------------
agg <- do.call(rbind, lapply(split(res[res$ok, ], list(res$lambda_t[res$ok],
                                                       res$n[res$ok])),
  function(d) {
    if (!nrow(d)) return(NULL)
    lam <- d$lambda_t[1]
    data.frame(lambda = lam, n = d$n[1], nrep = nrow(d),
               conv = mean(d$conv0 & d$conv2),
               reject = mean(d$reject),
               bias = mean(d$lam_hat) - lam,
               rmse = sqrt(mean((d$lam_hat - lam)^2)),
               cover = mean(abs(d$lam_hat - lam) <= 1.96 * d$lam_se),
               b1bias_M0 = mean(d$b1_0) - 0.25,
               b1bias_M2 = mean(d$b1_2) - 0.25,
               se_ratio_M0 = mean(d$b1se_0) / sd(d$b1_0),
               se_ratio_M2 = mean(d$b1se_2) / sd(d$b1_2),
               maxcor = mean(d$maxcor2, na.rm = TRUE))
  }))
agg <- agg[order(agg$lambda, agg$n), ]
print(agg, digits = 3, row.names = FALSE)
write.csv(agg, file.path(here, "sim_S1S2_summary.csv"), row.names = FALSE)
cat("SIM S1/S2 DONE\n")
