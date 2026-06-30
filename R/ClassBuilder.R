#' Construct a simulation-argument object
#'
#' Builds a `SimuArg.<Scenario>` object bundling the data-generating functions
#' and parameters used by [DataGenerating()].
#'
#' @param nrep Number of replicates.
#' @param N Sample size.
#' @param p Number of measured covariates.
#' @param p_U Number of unmeasured confounders.
#' @param Scenario Either `"exogenous"` or `"endogenous"`.
#' @param max_t Administrative censoring time.
#' @param theta True treatment effect.
#' @param unmeasured_Confounding,InitCovariates,InitAssignment,SurvTime,SwitchingTime,CensoringTime,TransferX
#'   Data-generating building-block functions.
#' @param Control List of generation controls (e.g. `json_save`, `save_path`, `Annotation`).
#' @param ... Additional scenario parameters (e.g. `gamma`, `alpha`, `beta`).
#' @return An object of class `paste0("SimuArg.", Scenario)`.
#' @export
New_SimuArg <- function(nrep, N, p, p_U, Scenario, max_t, theta,
                        unmeasured_Confounding = function(N, p_U, ...) NULL,
                        InitCovariates = function(N, p, ...) NULL,
                        InitAssignment = function(N, p, Covariates, ...) NULL,
                        SurvTime = function(N, p, Covariates, W, theta, ...) NULL,
                        SwitchingTime = function(N, p, Covariates, Z, ...) NULL,
                        CensoringTime = function(N, p, Covariates, ...) NULL,
                        TransferX = function(N, p, Covariates, ...) NULL,
                        Control = list(json_save = FALSE,
                                       save_path = "",
                                       Annotation = ""), ...)
{
  stopifnot(is.function(InitCovariates))
  stopifnot(is.function(InitAssignment))
  stopifnot(is.function(SurvTime))
  stopifnot(is.function(SwitchingTime))
  stopifnot(is.function(CensoringTime))
  Validate_scenario(Scenario)
  Control <- rlang::dots_list(!!!Control, json_save = FALSE, save_path = "",
                              Annotation = "", .homonyms = "first")
  data <- rlang::dots_list(initials = list(nrep = nrep, N = N, p = p, p_U = p_U,
                                           max_t = max_t, theta = theta),
                           unmeasured_Confounding = unmeasured_Confounding,
                           InitCovariates = InitCovariates,
                           InitAssignment = InitAssignment,
                           SurvTime = SurvTime,
                           SwitchingTime = SwitchingTime,
                           CensoringTime = CensoringTime,
                           TransferX = TransferX,
                           Control = Control,
                           parameters = list(...), .homonyms = "first")
  structure(data, class = paste0("SimuArg.", Scenario))
}


#' Construct a model-parameter object
#'
#' Builds a `ModelPar.<method>` object dispatched on by [DataFitting()].
#'
#' @param N Sample size.
#' @param p Number of measured covariates.
#' @param max_t Administrative censoring time.
#' @param dat A list with the observed data (covariates, treatment, times, etc.).
#' @param method Estimation method name (see [Validate_method()]).
#' @param Control List of estimation controls (grid, iterations, tolerances, etc.).
#' @param ... Additional method arguments (e.g. `ml_fitting_surv`, `ml_fitting_propensity`).
#' @return An object of class `paste0("ModelPar.", method[1])`.
#' @export
New_ModelPar <- function(N, p, max_t, dat, method,
                         Control = list(grid = 100,
                                        max_iter = 20,
                                        tol = 1e-5,
                                        contraction = 0.5,
                                        eta = 1e-4,
                                        init_parameters = runif(p + 1),
                                        learning_rate = 0.1),
                         ...)
{
  Validate_method(method)
  Control <- rlang::dots_list(!!!Control, grid = 100,
                              max_iter = 20,
                              tol = 1e-5,
                              contraction = 0.5,
                              eta = 1e-4,
                              learning_rate = 0.1,
                              init_parameters = runif(p + 1), .homonyms = "first")
  data <- rlang::dots_list(N = N, p = p, max_t = max_t, dat = dat, Control = Control,
                           !!!list(...), .homonyms = "first")
  structure(data, class = paste0("ModelPar.", method[1]))
}
