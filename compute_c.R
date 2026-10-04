## compute_c.R ---------------------------------------------------------------
## Exact computation of the plug-in LR inflation factor
##
##   LR  =>  c * chi^2_1,   c = 1 + Itilde_{ld} V_d Itilde_{dl} / Itilde_l
##
## for the latent-normality test with CML-estimated item difficulties
## treated as known (paper §6.1 remark).  All quantities are expected
## values at the null (lambda = 0), computed by numerical integration --
## no simulation.
##
## Blocks (per person):
##   I_full : expected information of the MARGINAL model over
##            eta = (nu, lambda, delta), nu the freely re-estimated
##            second-stage parameters (a, ln b) or (beta0, beta1, ln b),
##            via I = sum_s P0(s) g(s) g(s)' with analytic score
##            functions (delta-scores via leave-one-out score dists).
##   I_C    : CML information for delta, via theta-free conditional
##            moments E[Y_i | S=s], Cov(Y_i, Y_j | S=s).
##   Both are singular in the all-ones delta direction; a sum-to-zero
##   contrast basis C removes it (c is invariant to the basis choice).
##
## Itilde_l  = I_ll - I_ln Inn^-1 I_nl              (efficient lambda info)
## Itilde_ld = (I_ld - I_ln Inn^-1 I_nd) C          (efficient cross block)
## V_d       = (C' I_C C)^-1
##
## Run: Rscript compute_c.R          (validations first, then the designs)
## --------------------------------------------------------------------------
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "bcm.R"))

## uniform (trapezoid) grid on the basic scale -- the reliable default
make_grid_unif <- function(step = 0.05, lim = 7) {
  z <- seq(-lim, lim, by = step)
  w <- dnorm(z) * step
  list(z = z, w = w / sum(w))
}

## ===== marginal score functions at lambda = 0 ===============================
## groups: list of list(x = covariate row for nu-mapping a_g = x %*% beta,
##                      wt = population weight)
## nu = (beta..., lnb): a_g = x_g' beta, common b.
## Returns I_full over (nu, lambda, delta) and P0 per group.
marginal_info <- function(groups, beta, b, alpha, delta, grid) {
  k <- length(alpha); p <- length(beta)
  npar <- p + 1 + 1 + k                      # beta, lnb, lambda, delta
  I <- matrix(0, npar, npar)
  z <- grid$z; w <- grid$w
  for (g in groups) {
    a  <- sum(g$x * beta)
    th <- a + b * z
    SA  <- score_all_w(th, alpha, delta)     # P: (smax+1) x Q, ES: Q
    P   <- SA$P; ES <- SA$ES
    Pw  <- pmax(as.numeric(P %*% w), 1e-300)
    smax <- length(Pw) - 1L
    R   <- outer(0:smax, ES, "-")            # s - E(S|theta_q)
    PR  <- P * R
    ## nu and lambda scores: int [s-ES] dtheta/deta  (posterior-weighted)
    g_a   <- as.numeric(PR %*% w) / Pw                    # dtheta/da = 1
    g_lnb <- as.numeric(PR %*% (w * b * z)) / Pw          # dtheta/dlnb = b z
    g_lam <- as.numeric(PR %*% (w * b * z^2 / 2)) / Pw    # dtheta/dlam|0 = b z^2/2
    ## delta scores via leave-one-out score distributions
    Gd <- matrix(0, smax + 1L, k)
    for (i in seq_len(k)) {
      SAi <- score_all_w(th, alpha[-i], delta[-i])        # (smax-alpha_i+1) x Q
      Pi  <- SAi$P
      pi_q <- plogis(alpha[i] * (th - delta[i]))
      dpdd <- -alpha[i] * pi_q * (1 - pi_q)               # dp_i/ddelta_i, length Q
      ## dP(s)/ddelta_i = dpdd * [P^{(-i)}(s - alpha_i) - P^{(-i)}(s)]
      Pi_pad_lo <- rbind(matrix(0, alpha[i], length(z)), Pi)          # s - alpha_i
      Pi_pad_hi <- rbind(Pi, matrix(0, alpha[i], length(z)))          # s
      D <- Pi_pad_lo - Pi_pad_hi                                      # (smax+1) x Q
      Gd[, i] <- as.numeric((D * rep(dpdd, each = smax + 1L)) %*% w) / Pw
    }
    ## assemble G: columns = (beta_1..p via chain rule, lnb, lambda, delta)
    G <- cbind(outer(g_a, g$x), g_lnb, g_lam, Gd)
    I <- I + g$wt * crossprod(G * sqrt(Pw) , G * sqrt(Pw))
  }
  list(I = I, p = p)
}

## ===== CML information ======================================================
## Conditional moments E[Y_i|s], E[Y_i Y_j|s] are theta-free; we compute
## them via weighted elementary symmetric functions in LOG space, which is
## numerically stable at every score value (the naive probability-scale
## ratio underflows in the tails for long tests).
##   gamma_s = ESF of order s of the generating polynomial
##   prod_i (1 + u_i x^{alpha_i}),  u_i = exp(alpha_i (theta_ref - delta_i));
##   E[Y_i|s]      = u_i gamma^{(-i)}_{s-alpha_i} / gamma_s
##   E[Y_i Y_j|s]  = u_i u_j gamma^{(-ij)}_{s-alpha_i-alpha_j} / gamma_s
log_esf <- function(lu, alpha, smax) {
  ## returns log gamma_s, s = 0..smax (=-Inf where unattainable)
  lg <- c(0, rep(-Inf, smax))
  for (i in seq_along(lu)) {
    a <- alpha[i]
    shifted <- c(rep(-Inf, a), lg[seq_len(smax + 1L - a)]) + lu[i]
    m <- pmax(lg, shifted)
    lg <- ifelse(is.finite(m),
                 m + log(exp(lg - m) + exp(shifted - m)), -Inf)
  }
  lg
}
cml_info <- function(groups, beta, b, alpha, delta, grid, theta_ref = 0) {
  k <- length(alpha); smax <- sum(alpha)
  lu <- alpha * (theta_ref - delta)             # log u_i
  lg  <- log_esf(lu, alpha, smax)
  lgi <- lapply(seq_len(k), function(i)
           log_esf(lu[-i], alpha[-i], smax))    # padded to smax for indexing
  lEY <- matrix(-Inf, smax + 1L, k)             # log E[Y_i | s]
  for (i in seq_len(k)) {
    idx <- (alpha[i] + 1L):(smax + 1L)
    v <- lu[i] + lgi[[i]][idx - alpha[i]] - lg[idx]
    lEY[idx, i] <- ifelse(is.nan(v), -Inf, v)   # unattainable scores
  }
  EY <- exp(lEY)
  EYY <- array(0, c(smax + 1L, k, k))
  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    lgij <- log_esf(lu[-c(i, j)], alpha[-c(i, j)], smax)
    sh <- alpha[i] + alpha[j]
    idx <- (sh + 1L):(smax + 1L)
    lv <- lu[i] + lu[j] + lgij[idx - sh] - lg[idx]
    v  <- exp(ifelse(is.nan(lv), -Inf, lv))
    EYY[idx, i, j] <- EYY[idx, j, i] <- v
  }
  ## true marginal distribution of S (mixture over groups)
  P0 <- numeric(smax + 1L)
  for (g in groups) {
    a  <- sum(g$x * beta)
    SA <- score_all_w(a + b * grid$z, alpha, delta)
    P0 <- P0 + g$wt * as.numeric(SA$P %*% grid$w)
  }
  ## I_C = sum_s P0(s) * [alpha_i alpha_j Cov(Y_i, Y_j | s)]
  IC <- matrix(0, k, k)
  for (s in 0:smax) {
    if (P0[s + 1L] < 1e-14) next
    ey <- EY[s + 1L, ]
    CV <- EYY[s + 1L, , ] - tcrossprod(ey)
    diag(CV) <- ey * (1 - ey)
    IC <- IC + P0[s + 1L] * (outer(alpha, alpha) * CV)
  }
  list(IC = IC, EY = EY, P0 = P0)
}

## ===== the inflation factor =================================================
compute_c <- function(groups, beta, b, alpha, delta, grid) {
  k <- length(alpha)
  mi <- marginal_info(groups, beta, b, alpha, delta, grid)
  I  <- mi$I; p <- mi$p
  i_nu <- 1:(p + 1); i_l <- p + 2; i_d <- (p + 3):(p + 2 + k)
  Inn <- I[i_nu, i_nu]; Iln <- I[i_l, i_nu, drop = FALSE]
  Ild <- I[i_l, i_d, drop = FALSE];  Ind <- I[i_nu, i_d, drop = FALSE]
  Itil_l  <- I[i_l, i_l] - Iln %*% solve(Inn, t(Iln))
  Itil_ld <- Ild - Iln %*% solve(Inn, Ind)
  ## contrast basis orthogonal complement of all-ones (sum-to-zero coding)
  C <- rbind(diag(k - 1), -1)
  ci <- cml_info(groups, beta, b, alpha, delta, grid)
  Vd <- solve(t(C) %*% ci$IC %*% C)
  lc <- Itil_ld %*% C
  c_val <- 1 + as.numeric(lc %*% Vd %*% t(lc)) / as.numeric(Itil_l)
  list(c = c_val, Itil_l = as.numeric(Itil_l), I = I, IC = ci$IC)
}

grid <- make_grid_unif(0.05, 7)

## ===========================================================================
## VALIDATIONS (small design: k = 5 Rasch)
## ===========================================================================
cat("== validations (k = 5 Rasch, delta equally spaced [-2,2], a=0, b=1) ==\n")
k5 <- 5L; alpha5 <- rep(1L, k5); delta5 <- seq(-2, 2, length.out = k5)
grp1 <- list(list(x = 1, wt = 1))

## (1) marginal score functions vs numDeriv on log P(s)
mi5 <- marginal_info(grp1, beta = 0, b = 1, alpha5, delta5, grid)
logPs <- function(par) {   # par = (a, lnb, lambda, delta)
  th <- par[1] + exp(par[2]) * bcm_t(grid$z, par[3])
  SA <- score_all_w(th, alpha5, par[4:8])
  log(pmax(as.numeric(SA$P %*% grid$w), 1e-300))
}
J <- numDeriv::jacobian(logPs, c(0, 0, 0, delta5))
## rebuild our G for the same point
a <- 0; b <- 1; th <- a + b * grid$z
SA <- score_all_w(th, alpha5, delta5); Pw <- as.numeric(SA$P %*% grid$w)
R <- outer(0:k5, SA$ES, "-"); PR <- SA$P * R
G_an <- cbind(as.numeric(PR %*% grid$w) / Pw,
              as.numeric(PR %*% (grid$w * grid$z)) / Pw,
              as.numeric(PR %*% (grid$w * grid$z^2 / 2)) / Pw,
              sapply(seq_len(k5), function(i) {
                SAi <- score_all_w(th, alpha5[-i], delta5[-i])
                pq <- plogis(th - delta5[i]); d <- -pq * (1 - pq)
                Plo <- rbind(0, SAi$P); Phi <- rbind(SAi$P, 0)
                as.numeric(((Plo - Phi) * rep(d, each = k5 + 1)) %*% grid$w) / Pw
              }))
cat(sprintf("  (1) score functions vs numDeriv: max abs diff %.2e\n",
            max(abs(G_an - J))))
stopifnot(max(abs(G_an - J)) < 1e-6)

## (2) CML info vs numDeriv hessian of expected conditional negll
ci5 <- cml_info(grp1, 0, 1, alpha5, delta5, grid)
patterns <- as.matrix(expand.grid(rep(list(0:1), k5)))
cond_negll <- function(dl) {
  ## E over y ~ truth of -log P(y|s; dl), s-weights from true P0
  u  <- exp(-dl); val <- 0
  for (r in seq_len(nrow(patterns))) {
    y <- patterns[r, ]; s <- sum(y)
    ## true prob of y: integrate prod p^y (1-p)^(1-y) over theta
    py_t <- sapply(grid$z, function(z) prod(plogis(z - delta5)^y *
                     (1 - plogis(z - delta5))^(1 - y)))
    py <- sum(py_t * grid$w)
    gam_s <- sum(apply(patterns[rowSums(patterns) == s, , drop = FALSE], 1,
                       function(yy) prod(u^yy)))
    val <- val - py * (sum(-y * dl) - log(gam_s))
  }
  val
}
H5 <- numDeriv::hessian(cond_negll, delta5)
cat(sprintf("  (2) CML info vs numDeriv: max abs diff %.2e (max entry %.3f)\n",
            max(abs(H5 - ci5$IC)), max(abs(ci5$IC))))
stopifnot(max(abs(H5 - ci5$IC)) < 1e-5)

## (3) theta-invariance of conditional moments + score identity
ci5b <- cml_info(grp1, 0, 1, alpha5, delta5, grid, theta_ref = 0.7)
cat(sprintf("  (3) E[Y|s] theta-invariance: max diff %.2e | sum_i a_i E[Y_i|s] = s: %.2e\n",
            max(abs(ci5$EY - ci5b$EY)),
            max(abs(ci5$EY %*% alpha5 - 0:k5))))
stopifnot(max(abs(ci5$EY - ci5b$EY)) < 1e-10)

## (4) basis invariance of c
cc5 <- compute_c(grp1, 0, 1, alpha5, delta5, grid)
C2 <- qr.Q(qr(cbind(rep(1, k5), matrix(rnorm(k5 * (k5 - 1)), k5))))[, 2:k5]
mi <- cc5$I; p <- 1
i_nu <- 1:2; i_l <- 3; i_d <- 4:8
Iln <- mi[i_l, i_nu, drop = FALSE]; Inn <- mi[i_nu, i_nu]
Itl <- mi[i_l, i_l] - Iln %*% solve(Inn, t(Iln))
Itld <- (mi[i_l, i_d, drop = FALSE] -
         Iln %*% solve(Inn, mi[i_nu, i_d, drop = FALSE])) %*% C2
c_alt <- 1 + as.numeric(Itld %*% solve(t(C2) %*% cc5$IC %*% C2) %*% t(Itld)) /
             as.numeric(Itl)
cat(sprintf("  (4) basis invariance: c = %.6f vs %.6f\n", cc5$c, c_alt))
stopifnot(abs(cc5$c - c_alt) < 1e-8)

## (5) grid refinement
cc5f <- compute_c(grp1, 0, 1, alpha5, delta5, make_grid_unif(0.025, 7.5))
cat(sprintf("  (5) grid refinement: c = %.6f vs %.6f\n", cc5$c, cc5f$c))
stopifnot(abs(cc5$c - cc5f$c) < 1e-5)
cat("  ALL VALIDATIONS PASSED\n\n")

## ===========================================================================
## RESULTS
## ===========================================================================
cat("== Rasch designs, delta equally spaced [-2,2], a = 0, b = 1 ==\n")
for (k in c(5, 10, 20, 40, 60)) {
  al <- rep(1L, k); dl <- seq(-2, 2, length.out = k)
  cc <- compute_c(grp1, 0, 1, al, dl, grid)
  cat(sprintf("  k = %2d: c = %.4f  (size of nominal 5%% test: %.4f)  Itilde_l = %.4f\n",
      k, cc$c, pchisq(qchisq(.95, 1) / cc$c, 1, lower.tail = FALSE), cc$Itil_l))
}

cat("\n== S1 simulation design (k = 40 weighted, two groups) ==\n")
set.seed(4040)
kS <- 40L
alS <- sample(1:4, kS, replace = TRUE, prob = c(.2, .35, .3, .15))
dlS <- seq(-2, 2, length.out = kS)[sample(kS)]
grp2 <- list(list(x = c(1, 0), wt = 0.5), list(x = c(1, 1), wt = 0.5))
ccS <- compute_c(grp2, beta = c(0.5, 0.25), b = 1, alS, dlS, grid)
cat(sprintf("  c = %.4f  (size of nominal 5%% test: %.4f)  Itilde_l = %.4f\n",
    ccS$c, pchisq(qchisq(.95, 1) / ccS$c, 1, lower.tail = FALSE), ccS$Itil_l))

cat("\n== Application design (k = 60 weighted, replica items 1-60, gender) ==\n")
rep60 <- readRDS(file.path(here, "replica_results.rds"))$truth
alA <- rep60$alpha[1:60]; dlA <- rep60$delta[1:60]
grpA <- list(list(x = c(1, 0), wt = 0.507), list(x = c(1, 1), wt = 0.493))
ccA <- compute_c(grpA, beta = c(1.0, 0.1), b = 1, alA, dlA, grid)
cat(sprintf("  c = %.4f  (size of nominal 5%% test: %.4f)  Itilde_l = %.4f\n",
    ccA$c, pchisq(qchisq(.95, 1) / ccA$c, 1, lower.tail = FALSE), ccA$Itil_l))

saveRDS(list(rasch = "see log", S1 = ccS[c("c", "Itil_l")],
             app = ccA[c("c", "Itil_l")]),
        file.path(here, "compute_c_results.rds"))
cat("\nDONE\n")
