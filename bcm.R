## bcm.R -------------------------------------------------------------------
## Core machinery for paper 1: latent regression with a Box-Cox-Manly
## (shifted log-normal) latent distribution and an OPLM/2PLM measurement
## model with known item parameters.
##
## Model (three layers), per group g defined by the covariate pattern:
##   z ~ N(0,1)                                (basic variable)
##   theta = a_g + b_g * t(z; lambda_g),       t(z;l) = (exp(l z)-1)/l
##   a_g = x_g' beta                           (latent regression)
##   P(S = s | theta): weighted-score distribution of the OPLM
##                     (integer discriminations alpha_i, difficulties
##                      delta_i, fixed and known)
##
## Marginal likelihood by Gauss-Hermite quadrature on the z-scale
## (parameter-free integration limits, eq. (lik-z) of the paper).
## Analytic gradient via the exponential-family score identity
## d log P(s|theta)/d theta = s - E(S|theta).
##
## Model sequence M0-M5 (Table `models' of the paper) is encoded via
## index maps: bmap[g] gives the index of the ln-b parameter of group g,
## lmap[g] the index of the lambda parameter of group g (0 = fixed at 0).
## --------------------------------------------------------------------------

## ===== Gauss-Hermite for N(0,1) ============================================
gauss_hermite <- function(nq) {
  i  <- 1:(nq - 1)
  bb <- sqrt(i / 2)
  J  <- matrix(0, nq, nq)
  J[cbind(i, i + 1L)] <- bb
  J[cbind(i + 1L, i)] <- bb
  e  <- eigen(J, symmetric = TRUE)
  o  <- order(e$values)
  ## nodes/weights for integral wrt standard normal density
  list(z = sqrt(2) * e$values[o],
       w = (e$vectors[1, ]^2)[o])          # weights sum to 1
}

## ===== Box-Cox-Manly transform and lambda-derivative =======================
bcm_t <- function(z, lam) {
  if (abs(lam) < 1e-6) z + 0.5 * lam * z^2 + lam^2 * z^3 / 6
  else (exp(lam * z) - 1) / lam
}
bcm_dt <- function(z, lam) {                 # d t(z;lam) / d lam
  if (abs(lam) < 1e-4) z^2 / 2 + lam * z^3 / 3
  else { E <- exp(lam * z); (lam * z * E - E + 1) / lam^2 }
}

## ===== OPLM weighted-score distribution ====================================
## P(S=s | theta) for s = 0..sum(alpha), via convolution of per-item
## generating polynomials (1-p_i) + p_i x^{alpha_i}.  Returns the full
## vector; also E(S|theta) when asked.
score_dist_w <- function(theta, alpha, delta) {
  p    <- plogis(alpha * (theta - delta))
  poly <- 1
  for (i in seq_along(alpha)) {
    f <- numeric(alpha[i] + 1L)
    f[1L]            <- 1 - p[i]
    f[alpha[i] + 1L] <- p[i]
    poly <- as.numeric(convolve(poly, rev(f), type = "open"))
  }
  pmax(poly, 0)
}
escore_w <- function(theta, alpha, delta)    # E(S | theta)
  sum(alpha * plogis(alpha * (theta - delta)))

## Vectorized version: P(S=s | theta_q) for ALL quadrature nodes at once.
## Returns list(P = (smax+1) x Q matrix, ES = length-Q vector E(S|theta_q)).
score_all_w <- function(theta, alpha, delta) {
  Q    <- length(theta)
  smax <- sum(alpha)
  P    <- matrix(0, smax + 1L, Q); P[1L, ] <- 1
  ES   <- numeric(Q)
  for (i in seq_along(alpha)) {
    p  <- plogis(alpha[i] * (theta - delta[i]))       # length Q
    a  <- alpha[i]
    Pshift <- rbind(matrix(0, a, Q), P[seq_len(smax + 1L - a), , drop = FALSE])
    P  <- sweep(P, 2L, 1 - p, "*") + sweep(Pshift, 2L, p, "*")
    ES <- ES + a * p
  }
  list(P = P, ES = ES)
}

## Unweighted (raw, number-correct) score distribution, for fit plots.
score_dist_raw <- function(theta, alpha, delta) {
  p    <- plogis(alpha * (theta - delta))
  poly <- 1
  for (i in seq_along(p)) {
    poly <- as.numeric(convolve(poly, rev(c(1 - p[i], p[i])), type = "open"))
  }
  pmax(poly, 0)
}

## ===== model specification =================================================
## groups: data.frame/list with, per group g:
##   X[g, ]  : covariate row (first column = 1)
##   freq[[g]]: named or plain vector of weighted-score frequencies,
##              length sum(alpha)+1 (scores 0..sum(alpha))
## spec: list(bmap = integer vector (per group), lmap = integer vector
##   (per group, 0 => lambda fixed at 0)).  Parameter vector:
##   par = c(beta[1..p], lnb[1..max(bmap)], lambda[1..max(lmap)])
bcm_spec <- function(model = c("M0","M1","M2","M3","M4","M5"), ngroups = 2L) {
  model <- match.arg(model)
  stopifnot(ngroups == 2L)
  switch(model,
    M0 = list(bmap = c(1L,1L), lmap = c(0L,0L)),
    M1 = list(bmap = c(1L,2L), lmap = c(0L,0L)),
    M2 = list(bmap = c(1L,1L), lmap = c(1L,1L)),
    M3 = list(bmap = c(1L,1L), lmap = c(1L,2L)),
    M4 = list(bmap = c(1L,2L), lmap = c(1L,1L)),
    M5 = list(bmap = c(1L,2L), lmap = c(1L,2L)))
}

npar_bcm <- function(spec, p) p + max(spec$bmap) + max(c(spec$lmap, 0L))

unpack_bcm <- function(par, spec, p) {
  nb <- max(spec$bmap); nl <- max(c(spec$lmap, 0L))
  list(beta = par[seq_len(p)],
       lnb  = par[p + seq_len(nb)],
       lam  = if (nl > 0) par[p + nb + seq_len(nl)] else numeric(0))
}

## ===== log-likelihood and analytic gradient ================================
## Returns -logLik; attribute "gradient" = -gradient (for optim).
bcm_negll <- function(par, X, freq, alpha, delta, spec, gh,
                      want_grad = TRUE) {
  p  <- ncol(X); G <- nrow(X)
  pp <- unpack_bcm(par, spec, p)
  a  <- as.numeric(X %*% pp$beta)
  nll <- 0
  grad <- numeric(length(par))
  for (g in seq_len(G)) {
    bg   <- exp(pp$lnb[spec$bmap[g]])
    lg   <- if (spec$lmap[g] > 0L) pp$lam[spec$lmap[g]] else 0
    ns   <- freq[[g]]
    smax <- length(ns) - 1L
    ## quadrature accumulation (vectorized over nodes)
    tq   <- bcm_t(gh$z, lg)
    dtq  <- bcm_dt(gh$z, lg)
    thq  <- a[g] + bg * tq
    sa   <- score_all_w(thq, alpha, delta)
    P    <- sa$P; ES <- sa$ES
    Pw   <- pmax(as.numeric(P %*% gh$w), 1e-300)  # marginal P(S=s), all s
    nll  <- nll - sum(ns * log(Pw))
    if (want_grad) {
      ## dP(s)/d eta = sum_q w_q P(s|th_q) (s - ES_q) dtheta/deta
      resid <- outer(0:smax, ES, "-")          # (s - E(S|theta_q))
      cf    <- ns / Pw                          # zero for unobserved scores
      base  <- as.numeric(t(cf) %*% (P * resid))  # length Q
      d_a   <- sum(base * gh$w)
      d_lnb <- sum(base * gh$w * (bg * tq))
      d_lam <- sum(base * gh$w * (bg * dtq))
      grad[seq_len(p)] <- grad[seq_len(p)] - d_a * X[g, ]
      grad[p + spec$bmap[g]] <- grad[p + spec$bmap[g]] - d_lnb
      if (spec$lmap[g] > 0L)
        grad[p + max(spec$bmap) + spec$lmap[g]] <-
          grad[p + max(spec$bmap) + spec$lmap[g]] - d_lam
    }
  }
  if (want_grad) attr(nll, "gradient") <- grad
  nll
}

bcm_negll_grad <- function(par, ...) attr(bcm_negll(par, ..., want_grad = TRUE), "gradient")

## ===== moments helpers =====================================================
Et_bcm   <- function(lam) if (abs(lam) < 1e-8) 0 else (exp(lam^2/2) - 1)/lam
Vart_bcm <- function(lam) if (abs(lam) < 1e-8) 1 else (exp(2*lam^2) - exp(lam^2))/lam^2

## solve (a, b) from (mu, sig2) given lambda  [step 2, eq. (step2)]
ab_from_musig <- function(mu, sig2, lam) {
  b <- sqrt(sig2 / Vart_bcm(lam))
  a <- mu - b * Et_bcm(lam)
  c(a = a, b = b)
}

## ===== the three-step fitting routine ======================================
## X: G x p covariate matrix (groups); freq: list of score-frequency vectors.
fit_bcm <- function(X, freq, alpha, delta, model = "M0", nq = 48L,
                    lam_grid = seq(-0.5, 0.5, by = 0.05),
                    start = NULL, hessian = TRUE) {
  spec <- bcm_spec(model)
  p    <- ncol(X); gh <- gauss_hermite(nq)
  np   <- npar_bcm(spec, p)

  obj <- function(par) as.numeric(bcm_negll(par, X, freq, alpha, delta, spec, gh,
                                            want_grad = FALSE))
  grf <- function(par) bcm_negll_grad(par, X, freq, alpha, delta, spec, gh)

  if (is.null(start)) {
    if (max(spec$lmap) == 0L) {
      ## normal model: moment-based starting values from observed scores
      start <- c(rep(0, p), rep(0, max(spec$bmap)))
      ## crude: mean weighted score -> theta location via escore inverse
      start[1] <- 0.5
    } else {
      ## step 1: fit the corresponding normal model (same bmap, lmap = 0)
      m_norm <- if (max(spec$bmap) == 1L) "M0" else "M1"
      fit0   <- fit_bcm(X, freq, alpha, delta, model = m_norm, nq = nq,
                        hessian = FALSE)
      pp0    <- unpack_bcm(fit0$par, bcm_spec(m_norm), p)
      ## implied group means/sds under normal fit
      a0   <- as.numeric(X %*% pp0$beta)
      b0   <- exp(pp0$lnb[bcm_spec(m_norm)$bmap])
      ## step 2: grid search over common lambda
      best <- NULL
      for (lam in lam_grid) {
        if (abs(lam) < 1e-9) next
        ## per lambda-parameter: use the same lambda for all its groups
        lamv  <- rep(lam, max(spec$lmap))
        ## solve (a, b) per group to preserve group mean/variance
        ab    <- t(vapply(seq_len(nrow(X)), function(g)
                    ab_from_musig(a0[g], b0[g]^2, lam), numeric(2)))
        ## map back to parameters: beta from a's (X assumed full rank G = p)
        beta  <- solve(t(X) %*% X, t(X) %*% ab[, 1])
        lnb   <- log(tapply(ab[, 2], spec$bmap, mean))
        cand  <- c(as.numeric(beta), as.numeric(lnb), lamv)
        val   <- obj(cand)
        if (is.null(best) || val < best$val) best <- list(par = cand, val = val)
      }
      start <- best$par
    }
  }

  op <- optim(start, obj, grf, method = "BFGS",
              control = list(maxit = 500, reltol = 1e-12))
  out <- list(par = op$par, loglik = -op$value, npar = np, model = model,
              spec = spec, converged = (op$convergence == 0),
              X = X, alpha = alpha, delta = delta, nq = nq)
  if (hessian) {
    H <- tryCatch(numDeriv::hessian(obj, op$par), error = function(e) NULL)
    if (!is.null(H)) {
      V <- tryCatch(solve(H), error = function(e) NULL)
      out$vcov <- V
      out$se   <- if (!is.null(V)) sqrt(pmax(diag(V), 0)) else rep(NA_real_, np)
      out$pcor <- if (!is.null(V)) cov2cor(V) else NULL
    }
  }
  ## labelled parameters
  nb <- max(spec$bmap); nl <- max(c(spec$lmap, 0L))
  names(out$par) <- c(paste0("beta", seq_len(p) - 1L),
                      paste0("lnb", seq_len(nb)),
                      if (nl > 0) paste0("lambda", seq_len(nl)))
  out
}

## LR test between two nested fits
lr_test <- function(fit0, fit1) {
  stat <- 2 * (fit1$loglik - fit0$loglik)
  df   <- fit1$npar - fit0$npar
  c(LR = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE))
}

## ===== implied latent moments of a fitted model ============================
## mixture over groups with weights wg
latent_moments <- function(fit, wg) {
  p  <- ncol(fit$X)
  pp <- unpack_bcm(fit$par, fit$spec, p)
  a  <- as.numeric(fit$X %*% pp$beta)
  G  <- nrow(fit$X)
  mg <- vg <- numeric(G)
  for (g in seq_len(G)) {
    b  <- exp(pp$lnb[fit$spec$bmap[g]])
    l  <- if (fit$spec$lmap[g] > 0L) pp$lam[fit$spec$lmap[g]] else 0
    mg[g] <- a[g] + b * Et_bcm(l)
    vg[g] <- b^2 * Vart_bcm(l)
  }
  mean <- sum(wg * mg)
  var  <- sum(wg * (vg + mg^2)) - mean^2
  c(mean = mean, variance = var)
}

## ===== predicted raw-score distribution of a fitted model ==================
## Returns matrix (scores 0..k) x groups of predicted probabilities.
predicted_raw <- function(fit, nq = 96L) {
  gh <- gauss_hermite(nq)
  p  <- ncol(fit$X)
  pp <- unpack_bcm(fit$par, fit$spec, p)
  a  <- as.numeric(fit$X %*% pp$beta)
  k  <- length(fit$alpha)
  out <- matrix(0, k + 1L, nrow(fit$X))
  for (g in seq_len(nrow(fit$X))) {
    b <- exp(pp$lnb[fit$spec$bmap[g]])
    l <- if (fit$spec$lmap[g] > 0L) pp$lam[fit$spec$lmap[g]] else 0
    for (q in seq_along(gh$z)) {
      th <- a[g] + b * bcm_t(gh$z[q], l)
      out[, g] <- out[, g] + gh$w[q] * score_dist_raw(th, fit$alpha, fit$delta)
    }
  }
  out
}

## ===== data simulation =====================================================
sim_bcm_data <- function(n, alpha, delta, beta, bg, lg, pgirl = 0.493,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1  <- rbinom(n, 1L, pgirl)
  z   <- rnorm(n)
  a   <- beta[1] + beta[2] * x1
  th  <- numeric(n)
  for (g in 0:1) {
    idx <- which(x1 == g)
    th[idx] <- a[idx] + bg[g + 1] * vapply(z[idx], bcm_t, numeric(1),
                                           lam = lg[g + 1])
  }
  k   <- length(alpha)
  Pm  <- plogis(outer(th, delta, "-") * rep(alpha, each = n))
  Y   <- matrix(rbinom(n * k, 1L, Pm), n, k)
  sw  <- as.numeric(Y %*% alpha)   # weighted scores
  sr  <- rowSums(Y)                # raw scores
  list(x1 = x1, theta = th, Y = Y, sw = sw, sr = sr)
}

## frequency list (scores 0..sum(alpha)) by group from simulated data
freq_by_group <- function(sw, x1, smax) {
  lapply(0:1, function(g) {
    tab <- tabulate(sw[x1 == g] + 1L, nbins = smax + 1L)
    tab
  })
}
