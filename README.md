# DRIVE

[![R-CMD-check](https://github.com/Lyangix/DRIVE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Lyangix/DRIVE/actions/workflows/R-CMD-check.yaml)

Doubly Robust Instrumental Variable Estimation for survival data with treatment
switching. The package provides simulation tooling and a set of estimators for
additive-hazards models when patients may switch treatment during follow-up.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Lyangix/DRIVE")
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
| `DRIVE.joint` | **DRIVE joint**: joint estimation of the treatment effect and parametric survival nuisance coefficients, with a logistic propensity model. |
| `DRIVE.ML` | **DRIVE ML**: cross-fitted estimation using user-supplied propensity and survival nuisance learners. |

## Quick start

The package includes a fully synthetic, ready-to-fit dataset of 200 subjects.

```r
library(DRIVE)
data(drive_toy)

head(drive_toy$Covariates)
fit <- TRTSWE(drive_toy, max_t = 5,
              methods = c("ITT", "remove", "recensor", "TimeVar", "DRIVE.joint"))
fit
```

`drive_toy` is a list containing measured `Covariates`, initial treatment `Z`,
switching time `W`, observed time `T_D_c`, and the logical `event` indicator.
`W = Inf` means no switch; a finite switch can also occur after follow-up ends.
The true treatment effect is 0.1 and the follow-up limit is 5. 

Generate another dataset without creating files:

```r
dat <- SimulateData(N = 200, Scenario = "ii", model = "both", seed = 123)
fit <- TRTSWE(dat, max_t = 5, methods = c("ITT", "recensor"))
```

An explicit `seed` makes generation reproducible and preserves the caller's
random-number state. Estimation also uses random draws, so call `set.seed()`
before fitting when reproducible estimates are needed.

## Simulation scenarios

`R/Settings.R` defines the data-generating building blocks. The eight scenarios are available
through `New_SimuScenario()`:

| Argument | Options |
|---|---|
| `Scenario` | `"i"`: scenario i (exogenous switching); `"ii"`: scenario ii (endogenous switching) |
| `model` | `"linear"`, `"propensity"`, `"survival"`, `"both"` |

`model` identifies which components are nonlinear. Each preset generates two
measured uniform covariates and one unmeasured uniform confounder, assigns
initial treatment using a logistic propensity, generates a switching time,
and inverts a piecewise additive cumulative hazard for survival. Treatment
changes from `Z` to `1 - Z` after switching. Observed follow-up is the minimum
of the event time, random censoring, and `max_t`.

In scenario ii, the switching and survival generators share the same
latent exponential draw. The nonlinear-survival presets retain the original
linear reference time in the switching model. Nonpositive switching draws from
the original formula are treated as no switch. The presets retain the original
coefficients and default to smaller `N = 200`, `nrep = 1`; use `N = 1600` or
`3200` and `nrep = 1000` for the original study sizes. Use `New_SimuArg()` and
custom building blocks for other dimensions or coefficient settings.

```r
args <- New_SimuScenario(N = 200, nrep = 3, Scenario = "i", model = "linear")
set.seed(123)
generated <- DataGenerating(args, return_data = TRUE)
generated$summary

results <- SimuRun(args, methods = c("ITT", "recensor"), data = generated$data)
print(results, Comp_parameters = c(0.1, 0.25, 0.25))

# Generate and fit fresh replicates in one step.
set.seed(123)
results <- SimuRun(args, methods = c("ITT", "recensor"))
```

In-memory datasets keep unmeasured covariates separately in `U`, so they are
never included in `Covariates` passed to an estimator. They also retain latent
`T_D`, `T_0`, and `C` for diagnostics. `SimulateData()` and `drive_toy` contain
only the five observed-data components. The original `DataGenerating(args)`
summary-only return remains available.

For disk-based studies, enable JSON saving, generate first, then fit:

```r
args <- New_SimuScenario(
  N = 200, nrep = 2,
  Control = list(json_save = TRUE, save_path = tempdir(), Annotation = "_linear")
)
set.seed(123)
DataGenerating(args)
results <- SimuRun(args, methods = "ITT", sequence = 1:2)
```

JSON files retain the original covariate layout for compatibility; `SimuRun()`
uses only the measured columns. Unselected replicates are `NA` and excluded
from printed summaries. The installed example is available at:

```r
system.file("examples", "simulation.R", package = "DRIVE")
```

## Applying DRIVE ML

DRIVE handles the cross-fitting. Supply two functions that fit on the training
folds and return predictions for the held-out rows. The installed example uses
an `rpart` classification tree and a `randomForestSRC` survival forest:

```r
# Install these optional learner packages once:
# install.packages(c("rpart", "randomForestSRC"))
library(DRIVE)
source(system.file("examples", "drive_ml_learners.R", package = "DRIVE"))
data(drive_toy)

set.seed(123)
fit <- TRTSWE(
  drive_toy, max_t = 5, methods = c("DRIVE.joint", "DRIVE.ML"),
  ml_fitting_propensity = ml_propensity_tree,
  ml_fitting_surv = ml_survival_forest,
  nfolds = 5, seed = 123,
  Control = list(max_iter = 50, tol = 1e-5)
)
fit
fit$DRIVE.ML$Estim$Convergence
```

The complete fitting script is in
[`inst/examples/drive_ml.R`](inst/examples/drive_ml.R), also installed at
`system.file("examples", "drive_ml.R", package = "DRIVE")`.
The `source()` call above loads both callback definitions from
[`inst/examples/drive_ml_learners.R`](inst/examples/drive_ml_learners.R).
For example, this is the propensity callback used in that script:

```r
ml_propensity_tree <- function(data, predictx) {
  if (!requireNamespace("rpart", quietly = TRUE)) stop("Install rpart first.")
  # Train on this fold's training subjects; IV is the target.
  data$IV <- factor(data$IV, levels = c(0, 1))
  if (length(unique(data$IV)) < 2L) {
    stop("Both treatment arms are needed in each training fold.")
  }
  model <- rpart::rpart(
    IV ~ ., data = data, method = "class",
    control = rpart::rpart.control(cp = 0.01, minbucket = 10)
  )
  # predictx has held-out covariates only, with no IV column.
  probability <- stats::predict(model, newdata = predictx, type = "prob")[, "1"]
  # Example clipping; choose stabilization within the training folds.
  pmin(0.99, pmax(0.01, as.numeric(probability)))
}
```

Pass this function as `ml_fitting_propensity = ml_propensity_tree`, as in the
fit call above. The same learner file contains `ml_survival_forest`, including
the conversion of survival-forest predictions to the hazard increments DRIVE
requires.

Use `?DRIVE_ML` for the detailed interface. To substitute your own learners:

| Callback | Training input (`data`) | Prediction input (`predictx`) | Required return |
|---|---|---|---|
| `ml_fitting_propensity(data, predictx)` | Data frame with binary `IV` and propensity covariates | Held-out propensity covariates only | Numeric vector of `P(IV = 1)` in `[0, 1]`, one value per prediction row |
| `ml_fitting_surv(data, predictx)` | List with `time`, `event`, `Covariates`, `Covariates2`, `IV`, `D_status`, and the shared grid `stime` | Held-out survival covariates only | Numeric matrix of untreated conditional **hazard increments**, with prediction subjects in rows and `stime` in columns |

The survival return is a matrix. If a learner
predicts cumulative hazards `H`, evaluate them on `data$stime` and return
`H(t1) - H(0), H(t2) - H(t1), ...`, with `H(0) = 0`. Survival probabilities,
cumulative hazards, and unscaled hazard rates are not interchangeable with
these increments. The example maps the forest's `chf` predictions from
`time.interest` onto this grid before differencing; these outputs are described
in the [randomForestSRC prediction documentation](https://www.randomforestsrc.org/reference/predict.rfsrc.html).

Fit preprocessing, tuning and learners inside each training fold and preserve
prediction-row order. Both callbacks use `dat$Covariates` by default;
`dat$Covariates2` can provide different propensity predictors. The survival
example follows the original scripts by training on subjects whose recorded
treatment history is entirely zero. Adapt this nuisance-training strategy to
your study's assumptions. The small toy example and its 100-tree forest
illustrate the interface.

The same learners work in a simulation:

```r
args <- New_SimuScenario(N = 400, nrep = 3, Scenario = "i")
set.seed(123)
results <- SimuRun(args, methods = c("DRIVE.joint", "DRIVE.ML"),
                   ml_fitting_propensity = ml_propensity_tree,
                   ml_fitting_surv = ml_survival_forest,
                   nfolds = 5, seed = 123)
```

DRIVE joint estimates treatment and parametric covariate coefficients together.
DRIVE ML estimates the treatment effect; nuisance-coefficient entries in a
combined table are `NA`.

## Compatibility

To use an existing ML algorithm, wrap its training and prediction steps in a
function with arguments `data` and `predictx`, as in the `ml_propensity_tree`
definition above. DRIVE calls this function separately for each fold: with
`nfolds = 5`, it supplies four folds as training `data` and the remaining fold's
covariates as `predictx`. The function fits a model on `data` and returns
predictions for `predictx`. Reusing a model trained on the entire analysis
dataset would not provide this cross-fitting, because it has already seen the
held-out subjects.

For the propensity callback, training `data` still includes `IV` as the target;
`predictx` contains only covariates. When updating an older callback, remove any
assumption that `predictx` also contains held-out `IV` labels. Return one
probability per prediction row. The survival callback returns a hazard-increment
matrix directly. Malformed callback outputs produce an explicit error before
the estimating equation is called.

The research scripts remain in `scripts/SimuArgs.R` and `scripts/RealData.R`.
They are excluded from the package build and require additional data or
optional ML packages.
