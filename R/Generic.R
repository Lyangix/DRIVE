#' Generate simulation data
#'
#' Generate replicate datasets using the supplied building blocks. By default,
#' return the original summary rates; optionally retain datasets in memory.
#'
#' @param SimuArg A simulation argument object from [New_SimuArg()] or
#'   [New_SimuScenario()].
#' @param return_data Logical; return datasets as well as summary rates?
#' @return With `return_data = FALSE`, a list of mean treatment, switching and
#'   censoring rates, plus `anyNegetive` (legacy spelling), the total number of
#'   negative switching draws across replicates. With `return_data = TRUE`, a
#'   list with `summary` and `data` (one dataset per replicate). Each in-memory
#'   dataset contains only measured `Covariates`; unmeasured covariates are in
#'   `U`. Latent `T_D`, `T_0`, and `C` are retained for simulation diagnostics.
#' @details
#' `Control$json_save = TRUE` also saves datasets under
#' `file.path(save_path, "DataGenerated", paste0(N, Annotation))`.
#' For compatibility, saved JSON files retain the original layout with measured
#' and unmeasured covariates together; [SimuRun()] removes the unmeasured columns.
#' `adcensoring_rate` retains the original definition `mean(T_D > max_t)`,
#' regardless of whether random censoring occurs first.
#' @examples
#' set.seed(123)
#' generated <- DataGenerating(New_SimuScenario(N = 100), return_data = TRUE)
#' generated$summary
#' head(generated$data[[1]]$Covariates)
#' @md
#' @export
DataGenerating <- function(SimuArg, return_data = FALSE) {
  UseMethod("DataGenerating")
}

#' @rdname DataGenerating
#' @export
DataGenerating.SimuArg.i <- function(SimuArg, return_data = FALSE) {
  generate_simulation(SimuArg, return_data, endogenous = FALSE)
}

#' @rdname DataGenerating
#' @export
DataGenerating.SimuArg.ii <- function(SimuArg, return_data = FALSE) {
  generate_simulation(SimuArg, return_data, endogenous = TRUE)
}

# Compatibility for previously saved simulation argument objects.
DataGenerating.SimuArg.exogenous <- DataGenerating.SimuArg.i
DataGenerating.SimuArg.endogenous <- DataGenerating.SimuArg.ii


#' Fit a treatment-switching model
#'
#' Generic that dispatches on the `ModelPar.<method>` class. Supported methods
#' are `ITT`, `remove`, `recensor`, `TimeVar`, `DRIVE.joint` and `DRIVE.ML`.
#' See [DRIVE_ML] for the user-supplied learner interface.
#'
#' @param ModelPar A `ModelPar` object from [New_ModelPar()].
#' @return A list with at least `Coef` and `Var`.
#' @export
DataFitting <- function(ModelPar) {
  UseMethod("DataFitting")
}


#' @rdname DataFitting
#' @export
DataFitting.ModelPar.ITT <- function(ModelPar) {
  event <- ModelPar$dat$event
  surv <- survival::Surv(ModelPar$dat$T_D_c + runif(ModelPar$N, 0, 0.01),
                         event = event,
                         type = "right")
  mod <- ahaz::ahaz(surv, cbind(ModelPar$dat$Z,
                                ModelPar$dat$Covariates[, ]))
  return(list(Coef = coef(mod),
              Var = diag(vcov(mod))))
}


#' @rdname DataFitting
#' @export
DataFitting.ModelPar.remove <- function(ModelPar) {
  event <- ModelPar$dat$event
  event_remove <- ModelPar$dat$T_D_c < ModelPar$dat$W
  surv <- survival::Surv(ModelPar$dat$T_D_c[event_remove] + runif(sum(event_remove), 0, 0.01),
                         event = event[event_remove], type = "right")
  mod <- ahaz::ahaz(surv, cbind(ModelPar$dat$Z[event_remove],
                                ModelPar$dat$Covariates[event_remove, ]))
  return(list(Coef = coef(mod),
              Var = diag(vcov(mod))))
}


#' @rdname DataFitting
#' @export
DataFitting.ModelPar.recensor <- function(ModelPar) {
  event <- ModelPar$dat$event
  event_w <- event & ModelPar$dat$T_D_c < ModelPar$dat$W
  T_D_c_w <- ifelse(ModelPar$dat$T_D_c < ModelPar$dat$W, ModelPar$dat$T_D_c, ModelPar$dat$W)
  surv <- survival::Surv(T_D_c_w + runif(ModelPar$N, 0, 0.01), event = event_w, type = "right")
  mod <- ahaz::ahaz(surv, cbind(ModelPar$dat$Z,
                                ModelPar$dat$Covariates))
  return(list(Coef = coef(mod),
              Var = diag(vcov(mod))))
}


#' @rdname DataFitting
#' @export
DataFitting.ModelPar.TimeVar <- function(ModelPar) {
  # The formula is evaluated in this environment, including in a clean session
  # where survival has been loaded as a dependency but has not been attached.
  Surv <- survival::Surv
  const <- timereg::const
  event <- ModelPar$dat$event
  event_w <- ModelPar$dat$T_D_c > ModelPar$dat$W
  p <- ncol(ModelPar$dat$Covariates)
  colnames(ModelPar$dat$Covariates) <- paste0("X", 1:p)
  tvdat <- data.frame(id = 1:ModelPar$N, treatment = ModelPar$dat$Z,
                      ModelPar$dat$Covariates, event = event,
                      start_time = 0, end_time = ModelPar$dat$T_D_c)
  swdat <- tvdat[event_w, ]
  tvdat$end_time[event_w] <- ModelPar$dat$W[event_w]
  tvdat$event[event_w] <- FALSE
  swdat$start_time <- ModelPar$dat$W[event_w]
  swdat$treatment <- 1 - swdat$treatment
  adat <- rbind(tvdat, swdat)
  adat <- adat[order(adat$id, adat$start_time), ]
  if (any(adat$end_time == 0)) {
    nn <- length(adat$end_time[adat$end_time == 0])
    adat$end_time[adat$end_time == 0] <- runif(nn, 0, 0.01)
  }
  str_formula_aalen <- stringr::str_c(paste0("const(X", 1:(ModelPar$p), ")"), collapse = "+")
  str_formula_aalen <- paste0("Surv(start_time, end_time, event) ~ const(treatment) + ",
                              str_formula_aalen)
  mod <- timereg::aalen(formula(str_formula_aalen), data = adat,
                        max.time = ModelPar$max_t, id = adat$id)

  return(list(Coef = coef(mod)[, 1],
              Var = coef(mod)[, 2]^2))
}


#' @rdname DataFitting
#' @export
DataFitting.ModelPar.DRIVE.joint <- function(ModelPar) {
  ModelPar$Control <- rlang::dots_list(!!!ModelPar$Control,
                                       init_parameters = rep(0, ncol(ModelPar$dat$Covariates) + 1),
                                       .homonyms = "first")
  event <- ModelPar$dat$event
  T_D_c <- ModelPar$dat$T_D_c
  if (is.null(ModelPar$dat$stime)) {
    stime <- sort(T_D_c)
    stime <- unique(stime)
  } else {
    stime <- ModelPar$dat$stime
  }
  k <- length(stime)
  if (is.null(ModelPar$dat$D_status)) {
    D_status <- matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if (T_D_c[i] > ModelPar$dat$W[i]) {
        D_status[i, which(stime <= ModelPar$dat$W[i])] <- ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] <- 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] <- ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status <- ModelPar$dat$D_status
  }

  if (is.null(ModelPar$Covariates2)) ModelPar$Covariates2 <- ModelPar$dat$Covariates2
  if (is.null(ModelPar$Covariates2)) {
    args <- rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                             IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates,
                             Covariates2 = ModelPar$dat$Covariates,
                             D_status = D_status, stime = stime)
  } else {
    args <- rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                             IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates,
                             Covariates2 = ModelPar$Covariates2,
                             D_status = D_status, stime = stime)
  }
  mod <- easy_call(drive_joint_est_cpp, args)
  return(list(Coef = mod$x,
              Var = mod$var,            # joint semi-parametric sandwich variance
              Convergence = mod$Convergence,
              dLam = mod$dLam))
}


#' @rdname DataFitting
#' @export
DataFitting.ModelPar.DRIVE.ML <- function(ModelPar) {
  event <- ModelPar$dat$event
  T_D_c <- ModelPar$dat$T_D_c
  if (is.null(ModelPar$dat$stime)) {
    stime <- sort(T_D_c)
    stime <- unique(stime)
  } else {
    stime <- ModelPar$dat$stime
  }
  k <- length(stime)
  if (is.null(ModelPar$dat$D_status)) {
    D_status <- matrix(nrow = ModelPar$N, ncol = k)
    for (i in 1:ModelPar$N) {
      if (T_D_c[i] > ModelPar$dat$W[i]) {
        D_status[i, which(stime <= ModelPar$dat$W[i])] <- ModelPar$dat$Z[i]
        D_status[i, which(stime > ModelPar$dat$W[i])] <- 1 - ModelPar$dat$Z[i]
      } else {
        D_status[i, ] <- ModelPar$dat$Z[i]
      }
    }
  } else {
    D_status <- ModelPar$dat$D_status
  }

  if (!is.function(ModelPar$ml_fitting_surv)) stop("DRIVE.ML requires ml_fitting_surv; see ?DRIVE_ML.")
  if (!is.function(ModelPar$ml_fitting_propensity)) stop("DRIVE.ML requires ml_fitting_propensity; see ?DRIVE_ML.")
  if (is.null(ModelPar$nfolds)) ModelPar$nfolds <- 10
  if (is.null(ModelPar$seed)) ModelPar$seed <- 5884419

  if (is.null(ModelPar$Covariates2)) ModelPar$Covariates2 <- ModelPar$dat$Covariates2
  if (is.null(ModelPar$Covariates2)) {
    args <- rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                             IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates,
                             ml_fitting_surv = ModelPar$ml_fitting_surv,
                             ml_fitting_propensity = ModelPar$ml_fitting_propensity,
                             Covariates2 = ModelPar$dat$Covariates,
                             D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed, .homonyms = "first")
  } else {
    args <- rlang::dots_list(!!!ModelPar$Control, time = T_D_c, event = event,
                             IV = ModelPar$dat$Z, Covariates = ModelPar$dat$Covariates,
                             ml_fitting_surv = ModelPar$ml_fitting_surv,
                             ml_fitting_propensity = ModelPar$ml_fitting_propensity,
                             Covariates2 = ModelPar$Covariates2,
                             D_status = D_status, stime = stime, nfolds = ModelPar$nfolds, seed = ModelPar$seed, .homonyms = "first")
  }
  mod <- do.call(drive_ml_est_cpp, arg_filter(args, drive_ml_est_cpp))
  return(list(Coef = mod$x,
              Var = mod$var,
              Convergence = mod$Convergence))
}

# Compatibility for previously saved model objects.
DataFitting.ModelPar.DRIV.s <- DataFitting.ModelPar.DRIVE.joint
DataFitting.ModelPar.DRIV.cf.hz.ml.est <- DataFitting.ModelPar.DRIVE.ML


#' Print simulation results
#' @method print SimuResults
#' @param x A `SimuResults` object from [SimuRun()].
#' @param ... Optional `Comp_parameters` for bias computation.
#' @return Invisibly `NULL`; prints bias, SD and mean SE tables.
#' @export
print.SimuResults <- function(x, ...) {
  SimuResults <- x
  results <- rlang::dots_list(!!!SimuResults, !!!list(...),
                              Comp_parameters = rep(0, SimuResults$initials$p + 1),
                              .homonyms = "first")
  cat("Simulation Results for ")
  cat(normalize_methods(results$methods), ":\n")
  tb <- NULL
  tb2 <- NULL
  tb3 <- NULL
  selected <- if (is.null(results$sequence)) seq_len(results$initials$nrep) else results$sequence
  for (j in results$methods) {
    method <- normalize_methods(j)
    tb <- rbind(tb, apply(results$SimuResults[[j]]$Coef[, selected, drop = FALSE], 1, mean) - results$Comp_parameters)
    tb2 <- rbind(tb2, apply(results$SimuResults[[j]]$Coef[, selected, drop = FALSE], 1, sd))
    if (method %in% c("DRIVE.joint", "DRIVE.ML")) {
      tb3 <- rbind(tb3, c(mean(sqrt(results$SimuResults[[j]]$Var[selected])), rep(NA_real_, nrow(results$SimuResults[[j]]$Coef) - 1)))
    } else {
      tb3 <- rbind(tb3, apply(sqrt(results$SimuResults[[j]]$Var[, selected, drop = FALSE]), 1, mean))
    }
  }
  rownames(tb) <- normalize_methods(results$methods)
  colnames(tb) <- c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  rownames(tb2) <- normalize_methods(results$methods)
  colnames(tb2) <- c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  rownames(tb3) <- normalize_methods(results$methods)
  colnames(tb3) <- c("theta", paste0("alpha", 1:(SimuResults$initials$p)))
  cat("\t Mean bias or sampling mean: ", "\n")
  print.default(round(tb, 4), print.gap = 2L)
  cat("\n")
  cat("\t Standard deviation: ", "\n")
  print.default(round(tb2, 4), print.gap = 2L)
  cat("\n")
  cat("\t Mean standard error: ", "\n")
  print.default(round(tb3, 4), print.gap = 2L)
  cat("\n")
}


#' Print treatment-switching estimation results
#' @method print TRTSWE
#' @param x A `TRTSWE` object from [TRTSWE()].
#' @param all Logical; if `TRUE` print full coefficient / SE / p-value tables.
#' @param ... Unused.
#' @return Invisibly `NULL`; prints estimate tables.
#' @export
print.TRTSWE <- function(x, all = FALSE, ...) {
  Results <- x
  results <- rlang::dots_list(!!!Results, !!!list(...),
                              .homonyms = "first")

  cat("Results for ")
  cat(normalize_methods(names(Results)), ":\n")
  tb <- NULL
  tb3 <- NULL
  tb4 <- NULL
  p <- 0
  for (j in names(Results)) {
    if (length(results[[j]]$Estim$Coef) > p) p <- length(results[[j]]$Estim$Coef)
  }

  for (j in names(Results)) {
    method <- normalize_methods(j)
    if (method %in% c("DRIVE.joint", "DRIVE.ML")) {
      if (method == "DRIVE.joint") {
        tb <- rbind(tb, as.vector(results[[j]]$Estim$Coef))
      } else {
        tb <- rbind(tb, c(results[[j]]$Estim$Coef, rep(NA, p - 1)))
      }
      tb3 <- rbind(tb3, c(sqrt(results[[j]]$Estim$Var), rep(NA, p - 1)))
    } else {
      tb <- rbind(tb, results[[j]]$Estim$Coef)
      tb3 <- rbind(tb3, sqrt(results[[j]]$Estim$Var))
    }
  }
  rownames(tb) <- normalize_methods(names(Results))
  colnames(tb) <- c("theta", if (p > 1L) paste0("alpha", seq_len(p - 1L)))
  rownames(tb3) <- normalize_methods(names(Results))
  colnames(tb3) <- colnames(tb)

  if (all) {
    cat("\t Coef: ", "\n")
    print.default(round(tb, 4), print.gap = 2L)
    cat("\n")
    cat("\t Standard error: ", "\n")
    print.default(round(tb3, 4), print.gap = 2L)
    cat("\n")
    tb4 <- 2 * (1 - pnorm(abs(tb / tb3)))
    cat("\t P-value: ", "\n")
    print.default(round(tb4, 4), print.gap = 2L)
  } else {
    tb4 <- 2 * (1 - pnorm(abs(tb / tb3)))
    tb <- cbind(tb[, 1], tb3[, 1], tb4[, 1])
    colnames(tb) <- c("theta", "SE", "Pval")
    print.default(round(tb, 4), print.gap = 2L)
  }
}
