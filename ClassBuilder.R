New_SimuArg = function(nrep, N, p, p_U, Scenario, max_t, theta,
                       unmeasured_Confounding = function(N, p_U, ...) NULL,
                       InitCovariates = function(N, p, ...) NULL,
                       InitAssignment = function(N, p, Covariates, ...) NULL,
                       SurvTime = function(N, p, Covariates, W, theta,...) NULL,
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
  Control = rlang::dots_list(!!!Control, json_save = FALSE, save_path = "",
                             Annotation = "", .homonyms = "first")
  data = rlang::dots_list(initials = list(nrep = nrep, N = N, p = p, p_U = p_U, 
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


New_ModelPar = function(N, p, max_t, dat, method, 
                        Control = list(grid = 100,
                                       max_iter = 20,
                                       tol = 1e-5,
                                       contraction = 0.5,
                                       eta = 1e-4,
                                       init_parameters = runif(p+1),
                                       learning_rate = 0.1), 
                        ...)
{
  Validate_method(method)
  Control = rlang::dots_list(!!!Control, grid = 100,
                             max_iter = 20,
                             tol = 1e-5,
                             contraction = 0.5,
                             eta = 1e-4,
                             learning_rate = 0.1,
                             init_parameters = runif(p+1), .homonyms = "first")
  data = rlang::dots_list(N = N, p = p, max_t = max_t, dat = dat, Control = Control,
                          !!!list(...), .homonyms = "first")
  structure(data, class = paste0("ModelPar.", method[1]))
}

