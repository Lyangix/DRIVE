# R wrappers around the compiled estimating-equation solvers in src/.
# The compiled entry points driv_s_est() and driv_cf_ml_est() are
# made available through R/RcppExports.R (useDynLib).

# Partition indices into `nfolds` cross-fitting groups.
cf_group <- function(nfolds, datasize, seed) {
  if (!is.numeric(nfolds) || length(nfolds) != 1L || !is.finite(nfolds) ||
      nfolds != floor(nfolds) || nfolds < 2L || nfolds > datasize) {
    stop("nfolds must be an integer between 2 and the number of subjects.")
  }
  if (!is.null(seed)) set.seed(seed)
  n <- rep(seq_len(nfolds), ceiling(datasize / nfolds))[seq_len(datasize)]
  temp <- sample(n, datasize)
  x <- seq_len(nfolds)
  dataseq <- seq_len(datasize)
  cvlist <- lapply(x, function(x) dataseq[temp == x])
  return(cvlist)
}


#' DRIVE joint estimating-equation solver (Newton with backtracking line search)
#'
#' Wraps the compiled `driv_s_est()` solver. The instrument residual is formed
#' from a logistic propensity model of `IV` on `Covariates2`.
#'
#' @param init_parameters Initial parameter vector.
#' @param time Observed event/censoring times.
#' @param event Event indicator.
#' @param IV Instrumental variable (treatment assignment).
#' @param Covariates Confounders in the survival model (matrix).
#' @param Covariates2 Confounders in the propensity model.
#' @param D_status Treatment status at each grid time (matrix).
#' @param stime Ordered evaluation-time grid.
#' @param max_iter,tol,contraction,eta Newton / line-search controls.
#' @return A list with the estimate `x`, the variance `var` (the joint
#'   semi-parametric sandwich estimator), `var_orig` (the original
#'   scalar-sandwich estimator, for reference), `var_joint` (alias of `var`),
#'   `Convergence`, and `stime`.
#' @aliases driv_s_est_cpp
#' @export
drive_joint_est_cpp <- function(init_parameters, time, event, IV,
                              Covariates, Covariates2, D_status, stime, max_iter = 50, tol = 1e-5,
                              contraction = 0.5, eta = 1e-4) {
  mod <- glm(IV ~ Covariates2, family = binomial(link = "logit"))
  IV_c <- IV - expit(predict(mod))

  out <- driv_s_est(init_parameters = init_parameters, time = time,
                       event = event, IV = IV, IV_c = IV_c, Covariates = Covariates,
                       D_status = D_status, stime = stime, max_iter = max_iter,
                       tol = tol, contraction = contraction, eta = eta)
  out$stime <- stime
  return(out)
}


#' DRIVE ML cross-fitted estimating-equation solver
#'
#' Cross-fits user-supplied machine-learning nuisance estimators for the
#' propensity score and the untreated conditional hazard, then solves the compiled
#' `driv_cf_ml_est()` estimating equation.
#'
#' @inheritParams drive_joint_est_cpp
#' @param ml_fitting_surv Function called as `ml_fitting_surv(train_list, predictx)`.
#'   Return a numeric matrix of hazard increments with `nrow(predictx)` rows
#'   and `length(train_list$stime)` columns. See [DRIVE_ML] for the full contract.
#' @param ml_fitting_propensity Function fitting the propensity score; called as
#'   `ml_fitting_propensity(train_df, predictx)`. Return one probability of
#'   `IV = 1` per prediction row. Prediction data contain covariates only.
#' @param nfolds Number of cross-fitting folds, an integer from 2 to the number
#'   of subjects.
#' @param seed Random seed for fold assignment.
#' @param init_parameters Initial treatment-effect value; if a shared control
#'   vector is supplied, its first element is used.
#' @seealso [DRIVE_ML]
#' @md
#' @return A list with the estimate `x`, variance `var`, and `Convergence`.
#' @aliases driv_cf_ml_est_cpp
#' @export
drive_ml_est_cpp <- function(init_parameters, time, event, IV,
                                      Covariates, Covariates2, D_status, stime, ml_fitting_surv,
                                      ml_fitting_propensity,
                                      max_iter = 20, tol = 1e-5,
                                      contraction = 0.5, eta = 1e-4, nfolds = 10, seed = 5884419) {
  N <- length(time)
  if (!is.function(ml_fitting_surv) || !is.function(ml_fitting_propensity)) {
    stop("DRIVE.ML requires both nuisance-fitting functions; see ?DRIVE_ML.")
  }
  if (!is.numeric(stime) || !length(stime) || any(!is.finite(stime)) ||
      any(stime < 0) || is.unsorted(stime, strictly = TRUE)) {
    stop("stime must contain strictly increasing, finite, nonnegative times.")
  }
  if (!is.matrix(Covariates) || nrow(Covariates) != N ||
      !is.matrix(Covariates2) || nrow(Covariates2) != N ||
      !is.matrix(D_status) || !identical(dim(D_status), c(N, length(stime)))) {
    stop("Covariates, Covariates2 and D_status must be matrices with one row per subject; D_status must have one column per stime.")
  }
  if (!is.numeric(init_parameters) || !length(init_parameters) ||
      !is.finite(init_parameters[1L])) stop("init_parameters must start with a finite numeric value.")
  cflist <- cf_group(nfolds = nfolds, datasize = N, seed = seed)
  cf_IV_c <- rep(0, N)
  cf_surv <- matrix(0, nrow = N, ncol = length(stime))
  for (i in seq_along(cflist)) {
    cat("Fold ", i, "\n")
    ind <- cflist[[i]]
    tmp_df <- data.frame(IV = IV[-ind], Covariates2[-ind, , drop = FALSE])
    pred_df <- data.frame(Covariates2[ind, , drop = FALSE])
    mod <- ml_fitting_propensity(tmp_df, predictx = pred_df)
    if (!is.numeric(mod) || !is.null(dim(mod)) || length(mod) != length(ind) ||
        any(!is.finite(mod)) || any(mod < 0 | mod > 1)) {
      stop("ml_fitting_propensity must return a finite probability vector of length nrow(predictx) in [0, 1] (fold ", i, ").")
    }

    tmp_df_surv <- list(time = time[-ind],
                        event = event[-ind],
                        Covariates = Covariates[-ind, , drop = FALSE],
                        Covariates2 = Covariates2[-ind, , drop = FALSE],
                        IV = IV[-ind],
                        D_status = D_status[-ind, , drop = FALSE],
                        stime = stime)
    pred_df <- data.frame(Covariates[ind, , drop = FALSE])
    out <- ml_fitting_surv(tmp_df_surv, predictx = pred_df)
    if (!is.matrix(out) || !is.numeric(out) ||
        !identical(dim(out), c(length(ind), length(stime))) || any(!is.finite(out))) {
      stop("ml_fitting_surv must return a finite numeric matrix of hazard increments: nrow(predictx) by length(data$stime) (fold ", i, ").")
    }

    cf_IV_c[ind] <- IV[ind] - mod
    cf_surv[ind, ] <- out
  }

  out <- driv_cf_ml_est(init_parameters[1L], time = time, event = event, IV = IV, IV_c = cf_IV_c,
                                 D_status = D_status, stime = stime,
                                 ConfoundingPart = cf_surv, max_iter = max_iter,
                                 tol = tol, contraction = contraction, eta = eta)

  return(out)
}

# Compatibility aliases for the original exported R wrappers.
driv_s_est_cpp <- drive_joint_est_cpp
driv_cf_ml_est_cpp <- drive_ml_est_cpp
