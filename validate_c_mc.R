## validate_c_mc.R -----------------------------------------------------------
## Validation of compute_c.R beyond the k=5 Rasch checks:
##   (A) exact-side checks on a WEIGHTED design (score functions + CML info
##       vs numDeriv; pattern enumeration)
##   (B) Monte Carlo end-to-end: simulate under the null, CML-estimate the
##       item parameters, plug in, compute the LR for lambda; the mean LR
##       should equal c (and the control arm with TRUE item parameters
##       should give mean LR = 1).
## Arms: (1) k=5 Rasch single group; (2) k=6 weighted two-group.
## Env: REPS (default 20000), CORES (default 18).
## --------------------------------------------------------------------------
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "compute_c.R"))   # runs its own validations + designs

REPS  <- as.integer(Sys.getenv("REPS",  "20000"))
CORES <- as.integer(Sys.getenv("CORES", "18"))
gridMC <- make_grid_unif(0.1, 6)

## ===== (A) exact-side checks on a weighted design ==========================
cat("\n== (A) weighted-design exact checks (k=6, alpha=(1,2,1,3,1,2)) ==\n")
kw <- 6L; alw <- c(1L,2L,1L,3L,1L,2L); dlw <- seq(-1.5, 1.5, length.out = kw)
grpw <- list(list(x = c(1, 0), wt = 0.5), list(x = c(1, 1), wt = 0.5))
betw <- c(0.3, 0.25); bw <- 1

## (A1) marginal score functions vs numDeriv (group 1)
logPs_w <- function(par) {  # par = (beta0, beta1, lnb, lambda, delta)
  a <- par[1] + par[2] * 0
  th <- a + exp(par[3]) * bcm_t(grid$z, par[4])
  SA <- score_all_w(th, alw, par[5:10])
  log(pmax(as.numeric(SA$P %*% grid$w), 1e-300))
}
J <- numDeriv::jacobian(logPs_w, c(betw, 0, 0, dlw))
a <- betw[1]; th <- a + grid$z
SA <- score_all_w(th, alw, dlw); Pw <- as.numeric(SA$P %*% grid$w)
smaxw <- sum(alw)
R <- outer(0:smaxw, SA$ES, "-"); PR <- SA$P * R
Gd <- sapply(seq_len(kw), function(i) {
  SAi <- score_all_w(th, alw[-i], dlw[-i])
  pq <- plogis(alw[i] * (th - dlw[i])); d <- -alw[i] * pq * (1 - pq)
  Plo <- rbind(matrix(0, alw[i], length(grid$z)), SAi$P)
  Phi <- rbind(SAi$P, matrix(0, alw[i], length(grid$z)))
  as.numeric(((Plo - Phi) * rep(d, each = smaxw + 1)) %*% grid$w) / Pw
})
G_an <- cbind(as.numeric(PR %*% grid$w) / Pw,           # beta0 (x=1)
              0,                                        # beta1 (x=0 in group 1)
              as.numeric(PR %*% (grid$w * grid$z)) / Pw,
              as.numeric(PR %*% (grid$w * grid$z^2 / 2)) / Pw,
              Gd)
cat(sprintf("  (A1) weighted score functions vs numDeriv: max abs diff %.2e\n",
            max(abs(G_an - J))))
stopifnot(max(abs(G_an - J)) < 1e-6)

## (A2) CML info vs numDeriv via pattern enumeration (mixture-of-groups P0)
ciw <- cml_info(grpw, betw, bw, alw, dlw, grid)
patw <- as.matrix(expand.grid(rep(list(0:1), kw)))
cond_negll_w <- function(dl) {
  val <- 0
  for (g in grpw) {
    ag <- sum(g$x * betw)
    for (r in seq_len(nrow(patw))) {
      y <- patw[r, ]; s <- sum(y * alw)
      py_t <- sapply(grid$z, function(z)
        prod(plogis(alw * (ag + z - dl0))^y *
             (1 - plogis(alw * (ag + z - dl0)))^(1 - y)))
      ## true prob of pattern under group g (dl0 = truth, fixed below)
      py <- sum(py_t * grid$w)
      lg <- log_esf(alw * (0 - dl), alw, smaxw)     # log gamma at theta_ref=0
      val <- val - g$wt * py * (sum(y * alw * (0 - dl)) - lg[s + 1L])
    }
  }
  val
}
dl0 <- dlw
Hw <- numDeriv::hessian(cond_negll_w, dlw)
cat(sprintf("  (A2) weighted CML info vs numDeriv: max abs diff %.2e (max entry %.3f)\n",
            max(abs(Hw - ciw$IC)), max(abs(ciw$IC))))
stopifnot(max(abs(Hw - ciw$IC)) < 1e-5)
cat("  (A) PASSED\n")

## exact c for the two MC arms
cc_r <- compute_c(list(list(x = 1, wt = 1)), 0, 1,
                  rep(1L, 5), seq(-2, 2, length.out = 5), grid)
cc_w <- compute_c(grpw, betw, bw, alw, dlw, grid)
cat(sprintf("\nexact c: Rasch k=5 arm %.4f | weighted k=6 arm %.4f\n",
            cc_r$c, cc_w$c))

## ===== (B) Monte Carlo end-to-end ==========================================
## marginal negll for grouped score frequencies, params (beta..., lnb[,lam])
negll_marg <- function(par, freqs, Xg, alpha, delta, gh, est_lam) {
  p <- ncol(Xg); nll <- 0
  lam <- if (est_lam) par[p + 2] else 0
  for (g in seq_len(nrow(Xg))) {
    a <- sum(Xg[g, ] * par[1:p])
    th <- a + exp(par[p + 1]) * bcm_t(gh$z, lam)
    SA <- score_all_w(th, alpha, delta)
    Pw <- pmax(as.numeric(SA$P %*% gh$w), 1e-300)
    nll <- nll - sum(freqs[[g]] * log(Pw))
  }
  nll
}
## CML negll from sufficient statistics (item totals + score counts)
negll_cml <- function(dfree, Ti, ns, alpha) {
  dl <- c(dfree, -sum(dfree))
  lg <- log_esf(alpha * (0 - dl), alpha, sum(alpha))
  ## sufficient statistic for delta_i is alpha_i * T_i (weighted score model)
  sum(dl * alpha * Ti) + sum(ns * lg[as.numeric(names(ns)) + 1L])
}

one_rep <- function(rep, arm) {
  if (arm == "rasch") {
    k <- 5L; al <- rep(1L, k); dl <- seq(-2, 2, length.out = k)
    Xg <- matrix(1, 1, 1); bet <- 0; n <- 2000L
    x1 <- rep(0L, n); a_i <- rep(0, n)
  } else {
    k <- 6L; al <- alw; dl <- dlw
    Xg <- rbind(c(1, 0), c(1, 1)); bet <- betw; n <- 3000L
    x1 <- rbinom(n, 1L, 0.5); a_i <- betw[1] + betw[2] * x1
  }
  th <- a_i + rnorm(n)
  P  <- plogis(outer(th, dl, "-") * rep(al, each = n))
  Y  <- matrix(rbinom(n * k, 1L, P), n, k)
  sw <- as.numeric(Y %*% al)
  ## CML estimate of delta (sum-zero)
  Ti <- colSums(Y) ; ns <- table(sw)
  op <- optim(dl[-k], negll_cml, Ti = Ti, ns = ns, alpha = al,
              method = "BFGS", control = list(reltol = 1e-12, maxit = 500))
  dhat <- c(op$par, -sum(op$par))
  smax <- sum(al)
  freqs <- lapply(seq_len(nrow(Xg)), function(g)
    tabulate(sw[x1 == (g - 1L)] + 1L, nbins = smax + 1L))
  fitLR <- function(dd) {
    p <- ncol(Xg)
    st0 <- c(bet, 0)
    f0 <- optim(st0, negll_marg, freqs = freqs, Xg = Xg, alpha = al,
                delta = dd, gh = gridMC, est_lam = FALSE,
                method = "BFGS", control = list(reltol = 1e-12, maxit = 500))
    f1 <- optim(c(f0$par, 0.05), negll_marg, freqs = freqs, Xg = Xg,
                alpha = al, delta = dd, gh = gridMC, est_lam = TRUE,
                method = "BFGS", control = list(reltol = 1e-12, maxit = 500))
    2 * (f0$value - f1$value)
  }
  c(LR_hat = fitLR(dhat), LR_true = fitLR(dl))
}

for (arm in c("rasch", "weighted")) {
  t0 <- proc.time()[3]
  res <- parallel::mclapply(seq_len(REPS), function(r) {
    set.seed(31000000 + r + ifelse(arm == "rasch", 0, 50000000))
    tryCatch(one_rep(r, arm), error = function(e) c(NA, NA))
  }, mc.cores = CORES)
  M <- do.call(rbind, res); M <- M[complete.cases(M) & M[,1] >= 0 & M[,2] >= 0, ]
  se <- apply(M, 2, sd) / sqrt(nrow(M))
  cexp <- if (arm == "rasch") cc_r$c else cc_w$c
  cat(sprintf("MC arm %-8s (%d reps ok, %.1f min): mean LR plug-in %.4f (se %.4f) [exact c %.4f] | true-delta %.4f (se %.4f) [expect 1]\n",
      arm, nrow(M), (proc.time()[3] - t0)/60,
      mean(M[,1]), se[1], cexp, mean(M[,2]), se[2]))
}
cat("MC DONE\n")
