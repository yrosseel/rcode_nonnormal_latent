## glmm_epilepsy.R ----------------------------------------------------------
## Paper §6.2: Poisson GLMM with a BCM (shifted log-normal) random intercept,
## applied to the epilepsy data of Thall & Vail (1990), MASS::epil.
##
## Model:  y_ij | theta_i ~ Poisson( exp(x_ij' beta + theta_i) )
##         theta_i = b * t(z_i; lambda),  z_i ~ N(0,1)   [a = 0: the fixed
##         intercept absorbs location; theta_i has median 0]
## Covariates (Breslow & Clayton 1993 model IV without visit trend):
##         lbase*trt + lage + V4
##
## Integration: ADAPTIVE Gauss-Hermite per subject (posterior mode +
## curvature on the z-scale), Q = 25 nodes -- necessary here because the
## large counts concentrate the posterior of z far from 0.
##
## M0: lambda = 0 (normal random intercept; cross-checked against lme4
##     glmer nAGQ = 25).  M_skew: lambda free.  1-df LR test.
## Run: Rscript glmm_epilepsy.R           (writes epil_results.rds)
## --------------------------------------------------------------------------
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "bcm.R"))

library(MASS)
data(epil)

subj  <- unique(epil$subject)
Xf    <- model.matrix(~ lbase * trt + lage + V4, data = epil)
y     <- epil$y
id    <- match(epil$subject, subj)
s_i   <- as.numeric(tapply(y, id, sum))
lfac  <- as.numeric(tapply(lgamma(y + 1), id, sum))
nsubj <- length(subj)

gh <- gauss_hermite(25L)

## theta(z) = b t(z; lam) and derivatives wrt z
th_f   <- function(z, b, lam) b * bcm_t(z, lam)
th_d1  <- function(z, b, lam) if (abs(lam) < 1e-8) rep(b, length(z)) else b * exp(lam * z)
th_d2  <- function(z, b, lam) if (abs(lam) < 1e-8) rep(0, length(z)) else b * lam * exp(lam * z)

## negative log marginal likelihood by adaptive GH
## par = c(beta, lnb [, lambda])
negll_glmm <- function(par, est_lambda) {
  p    <- ncol(Xf)
  beta <- par[seq_len(p)]
  b    <- exp(par[p + 1L])
  lam  <- if (est_lambda) par[p + 2L] else 0
  eta  <- as.numeric(Xf %*% beta)
  C_i  <- as.numeric(tapply(exp(eta), id, sum))
  A_i  <- as.numeric(tapply(y * eta, id, sum)) - lfac
  ## posterior mode per subject: maximize g(z) = s th(z) - C e^{th(z)} - z^2/2
  zm <- numeric(nsubj)
  for (it in 1:50) {
    th  <- th_f(zm, b, lam); d1 <- th_d1(zm, b, lam); d2 <- th_d2(zm, b, lam)
    Ce  <- C_i * exp(th)
    g1  <- (s_i - Ce) * d1 - zm
    g2  <- -Ce * d1^2 + (s_i - Ce) * d2 - 1
    g2  <- pmin(g2, -1e-8)                       # keep concave direction
    step <- g1 / g2
    step <- sign(step) * pmin(abs(step), 1)      # damp large steps
    zm   <- zm - step
    if (max(abs(g1)) < 1e-10) break
  }
  th  <- th_f(zm, b, lam); d1 <- th_d1(zm, b, lam); d2 <- th_d2(zm, b, lam)
  Ce  <- C_i * exp(th)
  sig <- 1 / sqrt(pmax(Ce * d1^2 - (s_i - Ce) * d2 + 1, 1e-8))
  ## adaptive nodes: z = zm_i + sig_i * z_q.  With g(z) = s th(z) - C e^th
  ## - z^2/2:  log L_i = A_i + log sig_i + lse_q[log w_q + z_q^2/2 + g(node)]
  ## (the sqrt(2pi) factors of prior and node kernel cancel).
  ZN  <- zm + outer(sig, gh$z)                 # nsubj x Q node matrix
  TH  <- b * (if (abs(lam) < 1e-8) ZN + 0.5*lam*ZN^2 else (exp(lam*ZN)-1)/lam)
  G   <- s_i * TH - C_i * exp(TH) - 0.5 * ZN^2
  M   <- sweep(G, 2, log(gh$w) + 0.5 * gh$z^2, "+")
  rmax <- apply(M, 1, max)
  -sum(rmax + log(rowSums(exp(M - rmax))) + log(sig) + A_i)
}

fit_glmm <- function(est_lambda, start) {
  op <- optim(start, negll_glmm, est_lambda = est_lambda, method = "BFGS",
              control = list(maxit = 1000, reltol = 1e-12))
  H  <- numDeriv::hessian(negll_glmm, op$par, est_lambda = est_lambda)
  V  <- tryCatch(solve(H), error = function(e) matrix(NA, length(op$par),
                                                      length(op$par)))
  nm <- c(colnames(Xf), "lnb", if (est_lambda) "lambda")
  list(par = setNames(op$par, nm), se = setNames(sqrt(pmax(diag(V), 0)), nm),
       vcov = V, pcor = tryCatch(cov2cor(V), error = function(e) NULL),
       loglik = -op$value, npar = length(op$par),
       converged = (op$convergence == 0))
}

## ---- verify the AGQ likelihood against a fine trapezoid grid ---------------
## (integrate() is not used here: it misses the narrow posterior peaks of
## the high-count subjects and is off by several hundredths.)
glm0   <- glm(y ~ lbase * trt + lage + V4, family = poisson, data = epil)
start0 <- c(coef(glm0), log(0.5))
chk_par <- c(coef(glm0), log(0.5), -0.2)
ll_agq  <- -negll_glmm(chk_par, TRUE)
eta  <- as.numeric(Xf %*% chk_par[1:ncol(Xf)])
C_i  <- as.numeric(tapply(exp(eta), id, sum))
A_i  <- as.numeric(tapply(y * eta, id, sum)) - lfac
zt   <- seq(-10, 10, by = 1e-3)
tht  <- exp(chk_par[ncol(Xf)+1]) * bcm_t(zt, chk_par[ncol(Xf)+2])
ll_int <- sum(vapply(seq_len(nsubj), function(i) {
  g <- s_i[i]*tht - C_i[i]*exp(tht) + dnorm(zt, log = TRUE)
  max(g) + log(sum(exp(g - max(g))) * 1e-3) + A_i[i]
}, numeric(1)))
cat(sprintf("AGQ check: agq %.6f vs trapezoid %.6f (diff %.2e)\n",
            ll_agq, ll_int, abs(ll_agq - ll_int)))

cat("\n== M0: normal random intercept ==\n")
f0 <- fit_glmm(FALSE, start0)
print(round(rbind(est = f0$par, se = f0$se), 4))
cat(sprintf("logLik = %.3f\n", f0$loglik))

if (requireNamespace("lme4", quietly = TRUE)) {
  fm <- lme4::glmer(y ~ lbase * trt + lage + V4 + (1 | subject),
                    data = epil, family = poisson, nAGQ = 25)
  cat("lme4 cross-check:\n")
  print(round(rbind(bcm_M0 = f0$par[seq_len(ncol(Xf))],
                    lme4   = lme4::fixef(fm)), 4))
  cat(sprintf("sd(RE): bcm %.4f vs lme4 %.4f\n",
              exp(f0$par["lnb"]), sqrt(unlist(lme4::VarCorr(fm))[1])))
}

cat("\n== M_skew: BCM random intercept ==\n")
grid <- seq(-0.9, 0.9, by = 0.1); grid <- grid[abs(grid) > 1e-9]
vals <- vapply(grid, function(l) negll_glmm(c(f0$par, l), TRUE), numeric(1))
l0   <- grid[which.min(vals)]
cat(sprintf("lambda grid minimum at %.2f\n", l0))
f1 <- fit_glmm(TRUE, c(f0$par, l0))
print(round(rbind(est = f1$par, se = f1$se), 4))
cat(sprintf("logLik = %.3f\n", f1$loglik))
if (!is.null(f1$pcor)) {
  cat(sprintf("max |parameter correlation| = %.3f\n",
              max(abs(f1$pcor[upper.tri(f1$pcor)]))))
}

lt <- c(LR = 2 * (f1$loglik - f0$loglik), df = 1,
        p = pchisq(2 * (f1$loglik - f0$loglik), 1, lower.tail = FALSE))
cat(sprintf("\nLR test of latent normality: LR = %.2f, df = 1, p = %.4f\n",
            lt["LR"], lt["p"]))

## random-effect densities for the figure
thgrid <- seq(-3.5, 3.5, length.out = 400)
dens_normal <- dnorm(thgrid, 0, exp(f0$par["lnb"]))
b1 <- exp(f1$par["lnb"]); l1 <- f1$par["lambda"]
u1 <- -b1 / l1
dens_skew <- ifelse(l1 * (thgrid - u1) > 0,
  1 / (sqrt(2 * pi) * abs(l1) * abs(thgrid - u1)) *
    exp(-(log(pmax(l1 * (thgrid - u1) / b1, 1e-300)))^2 / (2 * l1^2)), 0)

saveRDS(list(f0 = f0, f1 = f1, lr = lt,
             dens = data.frame(theta = thgrid, normal = dens_normal,
                               skew = dens_skew)),
        file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))),
          "epil_results.rds"))
cat("DONE\n")
