#' Applying cross-fitted DRIVE ML
#'
#' DRIVE has two estimation versions: `DRIVE.joint` jointly estimates the
#' treatment effect and parametric survival nuisance coefficients, using a
#' logistic propensity model; `DRIVE.ML` estimates the treatment effect using
#' cross-fitted predictions from two user-supplied nuisance learners.
#'
#' @section Fitting:
#' Pass `methods = "DRIVE.ML"`, `ml_fitting_propensity`, `ml_fitting_surv`,
#' `nfolds`, and `seed` to [TRTSWE()] or [SimuRun()]. DRIVE partitions subjects,
#' trains each learner on the other folds, and requests predictions for the
#' held-out subjects. Do preprocessing, tuning and model fitting within each
#' callback using only its training data. Preserve prediction-row order.
#'
#' @section Propensity callback:
#' `ml_fitting_propensity(data, predictx)` receives a training data frame with
#' binary `IV` and propensity covariates. `predictx` contains only the held-out
#' propensity covariates, with the same column names; held-out `IV` is not
#' supplied. Return a finite numeric vector of length `nrow(predictx)` giving
#' `P(IV = 1 | covariates)`, with values in `[0, 1]`. Return probabilities,
#' not class labels, logits or a two-column probability matrix.
#'
#' @section Survival callback:
#' `ml_fitting_surv(data, predictx)` receives a training list containing `time`,
#' `event`, `Covariates`, `Covariates2`, `IV`, and `D_status`, plus the common
#' evaluation grid `stime`. `D_status` has training subjects in rows and grid
#' times in columns. `predictx` is a data frame of held-out survival covariates;
#' held-out outcomes are not supplied. The grid is shared across folds.
#'
#' Return a finite numeric matrix with `nrow(predictx)` rows and
#' `length(data$stime)` columns, representing the untreated conditional hazard
#' increments. If `H(t | X)` is a predicted cumulative hazard, column 1 is
#' `H(stime[1] | X) - H(0 | X)` and column j is
#' `H(stime[j] | X) - H(stime[j-1] | X)`. Use `H(0 | X) = 0`.
#' Evaluate cumulative hazards on DRIVE's grid first, then difference adjacent
#' columns. Return the matrix itself, not a list, survival probabilities,
#' cumulative hazards, or unscaled hazard rates.
#'
#' By default both learners use `dat$Covariates`. Supply `dat$Covariates2` for a
#' different propensity-covariate matrix (or pass `Covariates2` through `...`).
#' DRIVE ML returns a scalar treatment effect and its variance; it does not
#' estimate parametric covariate coefficients. Inspect
#' `fit$DRIVE.ML$Estim$Convergence` after fitting.
#' The default initial treatment effect is 0.1; `Control$init_parameters`
#' overrides it (the first element is used if a vector is shared with DRIVE joint).
#'
#' @section Runnable example:
#' The installed `examples/drive_ml_learners.R` defines a classification-tree
#' propensity learner and a random-survival-forest hazard learner. The latter
#' follows the original scripts by using training subjects whose recorded
#' treatment history is entirely zero; adapt that nuisance-training strategy
#' to your study's assumptions. The example uses 100 trees for a short run.
#' Install the optional `rpart` and `randomForestSRC` packages first.
#'
#' Source `examples/drive_ml.R` to run a complete comparison of DRIVE joint and
#' DRIVE ML on [drive_toy]. The toy sample illustrates the API; use an appropriate
#' sample size and learner tuning for your analysis.
#'
#' @section Compatibility:
#' The old names `DRIV.s` and `DRIV.cf.hz.ml.est` remain aliases for
#' `DRIVE.joint` and `DRIVE.ML`. New results use the canonical names.
#' Existing propensity callbacks should predict from covariates only;
#' `predictx` no longer includes held-out treatment labels.
#'
#' @examples
#' \dontrun{
#' # install.packages(c("rpart", "randomForestSRC"))
#' source(system.file("examples", "drive_ml_learners.R", package = "DRIVE"))
#' data(drive_toy)
#' set.seed(123)
#' fit <- TRTSWE(drive_toy, max_t = 5, methods = "DRIVE.ML",
#'                ml_fitting_propensity = ml_propensity_tree,
#'                ml_fitting_surv = ml_survival_forest,
#'                nfolds = 5, seed = 123)
#' fit
#' fit$DRIVE.ML$Estim$Convergence
#' }
#' @seealso [TRTSWE()], [SimuRun()], [drive_ml_est_cpp()]
#' @name DRIVE_ML
#' @md
NULL
