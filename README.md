# R code: non-normal latent variables and random effects

R code for the paper

> N. D. Verhelst, Y. Rosseel and M. Vidal.
> *Non-normal latent variables and random effects: transformations,
> testability and likelihood inference.*

The latent variable is a Box–Cox–Manly transformation of a standard normal
variable, `theta = a + b * (exp(lambda * z) - 1) / lambda`, which gives a
shifted log-normal latent distribution with the normal at `lambda = 0`.
This repository contains everything needed to reproduce the computed
results of the paper: the simulation studies, the synthetic replica of the
first application, the second application, the inflation factor of the
two-stage likelihood-ratio test, and the figures.

## Requirements

R (the results were produced with R 4.6.0) and the packages `numDeriv`,
`MASS` and `lme4` (the last only for the cross-check in
`glmm_epilepsy.R`). The simulations use `parallel::mclapply`, which runs
in parallel on Linux and macOS only.

Every script locates its own directory, so run each one with `Rscript`
from anywhere, for example `Rscript sim_S3.R`.

## Scripts

| script | what it does | paper |
|---|---|---|
| `bcm.R` | core module: the transformation, the weighted-score kernel of the restricted 2PLM, models M0–M5 with analytic gradient, the three-step fit, LR test, latent moments, predicted score distribution, data simulation | Sections 3–4 |
| `validate_bcm.R` | checks of `bcm.R`: analytic gradient against `numDeriv`, parameter recovery, LR test | — |
| `replica_entreetoets.R` | synthetic replica of the Entreetoets latent regression (16,000 students, 2 × 60 items), models M0–M5 on both item blocks | Section 6.1, Figures 2–3 |
| `glmm_epilepsy.R` | Poisson mixed model with a BCM random intercept for the epilepsy data (`MASS::epil`), adaptive Gauss–Hermite quadrature | Section 6.2, Table 4 |
| `compute_c.R` | exact inflation factor *c* of the plug-in LR test with CML-estimated item parameters, and the efficient information for lambda | Corollary 1, Sections 4.5, 6.1, 7 |
| `validate_c_mc.R` | Monte Carlo validation of `compute_c.R` | Section 6.1 |
| `sim_S1S2.R` | simulation studies S1 (size, power, recovery) and S2 (misspecified normal model) | Section 7, Table 5 |
| `sim_S3.R` | simulation study S3: comparison with the three-parameter skew-normal latent distribution | Section 7, Table 6 |
| `figures.R` | all figures, from the saved results | Figures 1–3 |

## Saved results

The output of each script is included, so that the tables and figures can
be checked without re-running the simulations.

| file | written by |
|---|---|
| `replica_results.rds`, `replica_log.txt` | `replica_entreetoets.R` |
| `entreetoets_fake.csv`, `entreetoets_fake_items.csv` | `replica_entreetoets.R` (the synthetic data set itself: 16,000 rows × id, girl, 120 item responses, generating `theta_true`; and the item parameters per block) |
| `epil_results.rds` | `glmm_epilepsy.R` |
| `compute_c_results.rds`, `validate_c_mc_log.txt` | `compute_c.R`, `validate_c_mc.R` |
| `sim_S1S2.rds`, `sim_S1S2_summary.csv`, `sim_S1S2_log.txt` | `sim_S1S2.R` |
| `sim_S3.rds`, `sim_S3_summary.csv`, `sim_S3_log.txt` | `sim_S3.R` |

## Running the simulations

The simulation scripts read their settings from environment variables:

| variable | meaning | default |
|---|---|---|
| `REPS` | replications per cell | 500 |
| `NQ` | Gauss–Hermite nodes | 96 (`sim_S1S2.R`), 192 (`sim_S3.R`) |
| `CORES` | parallel workers | 18 |
| `NFINE` | `sim_S3.R` only: replications per cell refitted with `2 * NQ` nodes | 25 |

For example, a quick trial run: `REPS=5 CORES=4 Rscript sim_S3.R`.

With the defaults on 18 cores, `sim_S1S2.R` took about 40 minutes and
`sim_S3.R` about 75 minutes. Seeds are fixed per replication, so a run
with the same settings regenerates the same data sets.

## Numerical integration: practical notes

`fit_bcm()` integrates over the basic normal variable with a fixed Gauss–Hermite rule
(`nq` nodes). Default node counts of other software (20–50) are **not** enough for this
model: on the synthetic replica, model M5 with 48 nodes has the boys' scale and skewness off
by 3.6 and 4.4 standard errors and the LR statistic against M0 off by 14%; 96 nodes still
leave λ one SE short; 180 and 360 nodes agree to three decimals. The error oscillates in the
number of nodes, so "add nodes until it stabilizes" is a weak check.

Why: of a 99-point rule only 45 nodes carry weight above 1e-12 (spacing 0.32 at the centre,
the rest beyond |z| ≈ 7), while the posterior of a unit is narrowest in the steep tail of the
transformation (precision 1 + b²I(θ)e^{2λz}) and the information about λ sits in the tails.

Recommended checks and alternatives:

- Refit once with a rule of a *different type* and confirm the estimates do not move.
- Uniform (trapezoid) grid on the z-scale: nodes `z = seq(-L, L, by = h)`, weights
  `dnorm(z) * h` normalized to sum to one, with L = 8–10 and h = 0.05 (321–401 nodes).
  The integrand is bounded by φ(z) for discrete data, so truncation at L = 8 costs < 1e-14;
  for a posterior of SD s the relative error is of order exp(−2π²s²/h²). A drop-in
  replacement for the Gauss–Hermite nodes/weights (see `make_grid_unif()` in
  `compute_c.R`); the finite-mixture likelihood keeps the analytic gradient exact.
- Adaptive Gauss–Hermite (as in `glmm_epilepsy.R`, 25 nodes) when unit-level information is
  large, e.g. individual-level counts.

## Data

- The Entreetoets data of the first application are proprietary to the Dutch
  National Institute of Educational Measurement (Cito) and are not included.
  `replica_entreetoets.R` simulates a data set with the same design and known
  parameter values and saves it as `entreetoets_fake.csv` (one row per
  simulated student: `id`, `girl`, items `i001`–`i120`, and the generating
  latent value `theta_true`) with the item parameters in
  `entreetoets_fake_items.csv`; the figures of the paper are drawn from this
  replica. The numbers of the original analysis quoted in the paper cannot be
  reproduced from this repository.
- The epilepsy data are the data set `epil` in the R package `MASS`.

## License

GPL-3; see `LICENSE`.
