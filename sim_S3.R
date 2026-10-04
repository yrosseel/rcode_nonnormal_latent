## sim_S3.R -----------------------------------------------------------------
## Simulation study S3: head-to-head with the three-parameter skew-normal
## latent distribution (Azzalini / eq. (sn-basic) of the paper, with
## location and scale).
## Same generating design as S1 (BCM truth), cells:
##   lambda_true in {0, -0.25} x n in {2000, 8000};  REPS (default 500).
## Per replicate, fit (i) BCM M2 and (ii) the skew-normal latent regression
## (mu, ln sigma, alpha) by ML; record the LR tests of normality of both
## models, the estimate and standard error of the skewness parameter, the
## Hessian conditioning, and max |parameter correlation|.
## The first NFINE replicates of each cell are refitted with 2*NQ nodes as
## a quadrature check.  NQ defaults to 192 here, not the 96 of S1/S2: with
## 96 nodes the skew-normal LR statistic is off by up to 1% and alpha-hat by
## 0.05-0.1 in this design; 192 nodes agree with 384 and with a fine
## trapezoid grid.
## Run: Rscript sim_S3.R                  (writes sim_S3.rds)
##
## 2026-10-03: the skew-normal arm previously implemented the four-parameter
## density of eq. (sn-four) as printed, which does not integrate to one; the
## results of that version are in archive/2026-10-03_pre-S3-fix/.
## --------------------------------------------------------------------------
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "bcm.R"))

REPS  <- as.integer(Sys.getenv("REPS",  "500"))
NQ    <- as.integer(Sys.getenv("NQ",    "192"))
CORES <- as.integer(Sys.getenv("CORES", "18"))
NFINE <- as.integer(Sys.getenv("NFINE", "25"))

set.seed(4040)
k     <- 40L
alpha <- sample(1:4, k, replace = TRUE, prob = c(.2, .35, .3, .15))
delta <- seq(-2, 2, length.out = k)[sample(k)]
beta_t <- c(0.5, 0.25)
X     <- rbind(c(1, 0), c(1, 1))
gh    <- gauss_hermite(NQ)
gh2   <- gauss_hermite(2L * NQ)

## ---- skew-normal latent regression (no covariate effect on shape) --------
## latent density: g(th) = (2/sig) phi(z) Phi(a z),  z = (th - mu)/sig,
## with the latent regression mu_g = x_g' beta.  Parameters:
## par = c(beta0, beta1, ln sigma, a);  a = 0 is the normal model M0.
sn_negll <- function(par, freq, gh) {
  beta <- par[1:2]; sig <- exp(par[3]); a <- par[4]
  mu   <- as.numeric(X %*% beta)
  selw <- 2 * pnorm(a * gh$z)          # selection on the standardized node
  nll  <- 0
  for (g in 1:2) {
    th   <- mu[g] + sig * gh$z
    sa   <- score_all_w(th, alpha, delta)
    Pw   <- pmax(as.numeric(sa$P %*% (gh$w * selw)), 1e-300)
    nll  <- nll - sum(freq[[g]] * log(Pw))
  }
  nll
}

## Quadrature error in the mean of the standardized skew-normal,
## E(Z) = sqrt(2/pi) a / sqrt(1 + a^2).  (The total mass sum(w * selw) is 1
## for every a by symmetry of the nodes, so it cannot serve as a check.)
sn_qerr <- function(a, gh)
  sum(gh$w * 2 * pnorm(a * gh$z) * gh$z) - sqrt(2 / pi) * a / sqrt(1 + a^2)
stopifnot(abs(sn_qerr(4, gh)) < 1e-5)

fit_sn <- function(freq, start, gh, hessian = TRUE) {
  op <- tryCatch(optim(start, sn_negll, freq = freq, gh = gh, method = "BFGS",
                       control = list(maxit = 500, reltol = 1e-10)),
                 error = function(e) NULL)
  if (is.null(op)) return(list(status = "error", loglik = -Inf))
  out <- list(par = op$par, loglik = -op$value,
              status = if (op$convergence == 0) "converged" else "maxit",
              mineig = NA_real_, hess_cond = NA_real_, maxcor = NA_real_,
              se = rep(NA_real_, length(start)))
  if (!hessian) return(out)
  H <- tryCatch(numDeriv::hessian(sn_negll, op$par, freq = freq, gh = gh),
                error = function(e) NULL)
  if (!is.null(H) && all(is.finite(H))) {
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    out$mineig <- min(ev)
    if (min(ev) > 0) {
      out$hess_cond <- max(ev) / min(ev)
      V <- solve(H)
      pc <- cov2cor(V)
      out$se     <- sqrt(diag(V))
      out$maxcor <- max(abs(pc[upper.tri(pc)]))
    }
  }
  out
}

## alpha = 0 is a stationary point of the skew-normal likelihood, and the
## likelihood is typically bimodal in alpha near normality: start on both
## sides and keep the better fit (as in Molenaar, Dolan & Verhelst, 2010).
fit_sn2 <- function(freq, f0par, gh, hessian = TRUE) {
  fits <- lapply(c(-1, 1), function(a0)
            fit_sn(freq, c(f0par[1:3], a0), gh, hessian = hessian))
  ll   <- sapply(fits, function(f) f$loglik)
  fs   <- fits[[which.max(ll)]]
  fs$lldiff <- abs(diff(ll))           # gap between the two local solutions
  fs
}

cells <- expand.grid(lambda = c(0, -0.25), n = c(2000L, 8000L))

one_rep <- function(cell, rep) {
  lam <- cells$lambda[cell]; n <- cells$n[cell]
  b   <- 1 / sqrt(Vart_bcm(lam))
  dat <- sim_bcm_data(n, alpha, delta, beta_t, bg = c(b, b), lg = c(lam, lam),
                      pgirl = 0.5, seed = 7700000L + cell * 100000L + rep)
  freq <- freq_by_group(dat$sw, dat$x1, sum(alpha))
  ## normal fit: null model of both LR tests, and source of starting values
  f0 <- tryCatch(fit_bcm(X, freq, alpha, delta, "M0", nq = NQ, hessian = FALSE),
                 error = function(e) NULL)
  if (is.null(f0)) return(NULL)
  ## BCM arm
  fb <- tryCatch(fit_bcm(X, freq, alpha, delta, "M2", nq = NQ),
                 error = function(e) NULL)
  bcm_ok   <- !is.null(fb)
  bcm_conv <- bcm_ok && fb$converged
  bcm_maxcor <- if (bcm_ok && !is.null(fb$pcor))
                  max(abs(fb$pcor[upper.tri(fb$pcor)])) else NA_real_
  ## SN arm
  fs <- fit_sn2(freq, f0$par, gh)
  if (fs$status == "error") return(NULL)
  out <- data.frame(cell = cell, rep = rep, lambda_t = lam, n = n,
             bcm_conv = bcm_conv, bcm_maxcor = bcm_maxcor,
             bcm_lambda = if (bcm_ok) fb$par[4] else NA_real_,
             bcm_se_l   = if (bcm_ok && !is.null(fb$se)) fb$se[4] else NA_real_,
             LR_bcm     = if (bcm_ok) 2 * (fb$loglik - f0$loglik) else NA_real_,
             sn_status = fs$status,
             sn_hess_pd = is.finite(fs$mineig) && fs$mineig > 0,
             sn_mineig = fs$mineig, sn_cond = fs$hess_cond,
             sn_maxcor = fs$maxcor,
             sn_alpha = fs$par[4], sn_se_a = fs$se[4],
             sn_lldiff = fs$lldiff, sn_qerr = sn_qerr(fs$par[4], gh),
             LR_sn = 2 * (fs$loglik - f0$loglik),
             sn_alpha_fine = NA_real_, LR_sn_fine = NA_real_)
  ## quadrature check: refit normal and skew-normal with 2*NQ nodes
  if (rep <= NFINE) {
    f0f <- tryCatch(fit_bcm(X, freq, alpha, delta, "M0", nq = 2L * NQ,
                            hessian = FALSE), error = function(e) NULL)
    if (!is.null(f0f)) {
      fsf <- fit_sn2(freq, f0f$par, gh2, hessian = FALSE)
      if (fsf$status != "error") {
        out$sn_alpha_fine <- fsf$par[4]
        out$LR_sn_fine    <- 2 * (fsf$loglik - f0f$loglik)
      }
    }
  }
  out
}

grid <- expand.grid(rep = seq_len(REPS), cell = seq_len(nrow(cells)))
cat(sprintf("running %d replicates on %d cores\n", nrow(grid), CORES))
t0  <- proc.time()[3]
res <- parallel::mclapply(seq_len(nrow(grid)), function(i)
         tryCatch(one_rep(grid$cell[i], grid$rep[i]), error = function(e) NULL),
         mc.cores = CORES)
nfail <- sum(!vapply(res, is.data.frame, logical(1)))
res <- do.call(rbind, Filter(is.data.frame, res))
rownames(res) <- NULL
cat(sprintf("done in %.1f min; %d replicates failed outright\n",
            (proc.time()[3] - t0) / 60, nfail))
saveRDS(res, file.path(here, "sim_S3.rds"))

## Per-cell summary.  Replicates with a non-positive-definite skew-normal
## Hessian are counted (sn_nonpd), not dropped silently; sn_maxcor and
## sn_se_a are over the positive-definite replicates only.
crit <- qchisq(.95, 1)
qq   <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
agg <- do.call(rbind, lapply(split(res, list(res$lambda_t, res$n)), function(d) {
  if (!nrow(d)) return(NULL)
  data.frame(lambda = d$lambda_t[1], n = d$n[1], nrep = nrow(d),
             bcm_conv = mean(d$bcm_conv),
             bcm_maxcor = mean(d$bcm_maxcor, na.rm = TRUE),
             bcm_lam_sd = sd(d$bcm_lambda, na.rm = TRUE),
             bcm_lam_se = mean(d$bcm_se_l, na.rm = TRUE),
             bcm_rej = mean(d$LR_bcm > crit, na.rm = TRUE),
             sn_nonconv = sum(d$sn_status != "converged"),
             sn_nonpd = sum(!d$sn_hess_pd),
             sn_rej = mean(d$LR_sn > crit),
             LR_cor = cor(d$LR_sn, d$LR_bcm, use = "complete.obs"),
             LR_agree = mean((d$LR_sn > crit) == (d$LR_bcm > crit), na.rm = TRUE),
             sn_a_q25 = qq(d$sn_alpha, .25), sn_a_med = qq(d$sn_alpha, .5),
             sn_a_q75 = qq(d$sn_alpha, .75), sn_a_sd = sd(d$sn_alpha),
             sn_a_se = mean(d$sn_se_a, na.rm = TRUE),
             sn_a_se_med = qq(d$sn_se_a, .5),
             sn_mineig_med = qq(d$sn_mineig, .5),
             sn_cond_med = qq(d$sn_cond, .5), sn_cond_max = qq(d$sn_cond, 1),
             sn_maxcor = mean(d$sn_maxcor, na.rm = TRUE),
             sn_maxcor_med = qq(d$sn_maxcor, .5),
             sn_cor95 = mean(d$sn_maxcor > 0.95, na.rm = TRUE),
             sn_lldiff_med = qq(d$sn_lldiff, .5),
             sn_qerr_max = max(abs(d$sn_qerr)),
             fine_n = sum(!is.na(d$LR_sn_fine)),
             fine_dalpha = max(abs(d$sn_alpha_fine - d$sn_alpha), na.rm = TRUE),
             fine_dLR = max(abs(d$LR_sn_fine - d$LR_sn), na.rm = TRUE))
}))
agg <- agg[order(agg$lambda, agg$n), ]
print(t(format(agg, digits = 3)), quote = FALSE)
write.csv(agg, file.path(here, "sim_S3_summary.csv"), row.names = FALSE)
cat("SIM S3 DONE\n")
