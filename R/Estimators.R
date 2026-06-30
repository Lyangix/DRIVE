# R wrappers around the compiled estimating-equation solvers in src/.
# The compiled entry points driv_s_est() and driv_cf_ml_est() are
# made available through R/RcppExports.R (useDynLib).

# Partition indices into `nfolds` cross-fitting groups.
cf_group <- function(nfolds, datasize, seed) {
  cvlist <- list()
  if (!is.null(seed)) set.seed(seed)
  n <- rep(1:nfolds, ceiling(datasize / nfolds))[1:datasize]
  temp <- sample(n, datasize)
  x <- 1:nfolds
  dataseq <- 1:datasize
  cvlist <- lapply(x, function(x) dataseq[temp == x])
  return(cvlist)
}


#' DRIV.s estimating-equation solver (Newton with backtracking line search)
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
#' @param stime Grid of event times.
#' @param max_iter,tol,contraction,eta Newton / line-search controls.
#' @return A list with the estimate `x`, the variance `var` (the joint
#'   semi-parametric sandwich estimator), `var_orig` (the original
#'   scalar-sandwich estimator, for reference), `var_joint` (alias of `var`),
#'   `Convergence`, and `stime`.
#' @export
driv_s_est_cpp <- function(init_parameters, time, event, IV,
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


#' DRIV.cf.hz.ml.est cross-fitted, ML-nuisance estimating-equation solver
#'
#' Cross-fits user-supplied machine-learning nuisance estimators for the
#' propensity score and the conditional survival, then solves the compiled
#' `driv_cf_ml_est()` estimating equation.
#'
#' @inheritParams driv_s_est_cpp
#' @param ml_fitting_surv Function fitting the conditional survival; called as
#'   `ml_fitting_surv(train_list, predictx)`.
#' @param ml_fitting_propensity Function fitting the propensity score; called as
#'   `ml_fitting_propensity(train_df, predictx)`.
#' @param nfolds Number of cross-fitting folds.
#' @param seed Random seed for fold assignment.
#' @return A list with the estimate `x`, variance `var`, and `Convergence`.
#' @export
driv_cf_ml_est_cpp <- function(init_parameters, time, event, IV,
                                      Covariates, Covariates2, D_status, stime, ml_fitting_surv,
                                      ml_fitting_propensity,
                                      max_iter = 20, tol = 1e-5,
                                      contraction = 0.5, eta = 1e-4, nfolds = 10, seed = 5884419) {
  N <- length(time)
  cflist <- cf_group(nfolds = nfolds, datasize = N, seed = seed)
  cf_IV_c <- rep(0, N)
  cf_surv <- matrix(0, nrow = N, ncol = length(stime))
  for (i in 1:length(cflist)) {
    cat("Folder ", i, "\n")
    ind <- cflist[[i]]
    tmp_df <- data.frame(IV = IV[-ind], Covariates2[-ind, , drop = FALSE])
    pred_df <- data.frame(IV = IV[ind], Covariates2[ind, , drop = FALSE])
    mod <- ml_fitting_propensity(tmp_df, predictx = pred_df)

    tmp_df_surv <- list(time = time[-ind],
                        event = event[-ind],
                        Covariates = Covariates[-ind, , drop = FALSE],
                        Covariates2 = Covariates2[-ind, , drop = FALSE],
                        IV = IV[-ind],
                        D_status = D_status[-ind, ],
                        stime = stime)
    pred_df <- data.frame(Covariates[ind, , drop = FALSE])
    out <- ml_fitting_surv(tmp_df_surv, predictx = pred_df)

    cf_IV_c[ind] <- IV[ind] - mod
    cf_surv[ind, ] <- out
  }

  out <- driv_cf_ml_est(0.1, time = time, event = event, IV = IV, IV_c = cf_IV_c,
                                 D_status = D_status, stime = stime,
                                 ConfoundingPart = cf_surv, max_iter = max_iter,
                                 tol = tol, contraction = contraction, eta = eta)

  return(out)
}
