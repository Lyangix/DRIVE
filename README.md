# DRIVE

[![R-CMD-check](https://github.com/Lyangix/DRIVE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Lyangix/DRIVE/actions/workflows/R-CMD-check.yaml)

Doubly Robust Instrumental Variable Estimation for survival data with treatment
switching. The package provides simulation tooling and a set of estimators for
additive-hazards models when patients may switch treatment during follow-up.

## Installation

```r
# install.packages("remotes")
remotes::install_github("lyangix/drive")
```

The estimating equations are implemented in C++ (Rcpp / RcppArmadillo), so a
working C++ toolchain is required to build from source.

## Estimation methods

`DataFitting()` / `TRTSWE()` / `SimuRun()` dispatch on the method name:

| Method | Description |
|---|---|
| `ITT` | Intention-to-treat additive-hazards fit (ignores switching). |
| `remove` | Drop person-time after switching. |
| `recensor` | Censor at the switching time. |
| `TimeVar` | Time-varying treatment (Aalen additive model). |
| `DRIV.s` | Doubly robust IV estimator (parametric nuisances). |
| `DRIV.cf.hz.ml.est` | Cross-fitted DRIV with user-supplied machine-learning nuisance estimators. |

## Quick start

```r
library(DRIVE)

# Single-dataset estimation
fit <- TRTSWE(dat_DRIV, max_t = 5,
              methods = c("ITT", "remove", "recensor", "TimeVar",
                          "DRIV.s", "DRIV.cf.hz.ml.est"),
              ml_fitting_surv = my_surv_fitter,
              ml_fitting_propensity = my_propensity_fitter)
fit
```

See `scripts/SimuArgs.R` for a full simulation workflow and
`scripts/RealData.R` for an applied analysis example. The `scripts/` directory
is excluded from the build (`.Rbuildignore`) and is provided for reference.

## Repository layout

```
R/      package source (constructors, generics, estimators, wrappers)
src/    C++ estimating-equation solvers
scripts/ example / reproduction scripts (not part of the package build)
```

## Development notes

The Rcpp glue in `R/RcppExports.R` and `src/RcppExports.cpp` is hand-written.
After changing any `// [[Rcpp::export]]` signature in `src/`, regenerate it with:

```r
Rcpp::compileAttributes()
```
