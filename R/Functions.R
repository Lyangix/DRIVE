#' Run a simulation study
#'
#' Reads previously generated replicate datasets and fits each requested method,
#' collecting coefficients and variances into a `SimuResults` object.
#'
#' @param SimuArg A `SimuArg` object from [New_SimuArg()].
#' @param methods Character vector of method names (see [Validate_method()]).
#' @param Control List of estimation controls passed to [New_ModelPar()].
#' @param sequence Optional integer vector of replicate indices; defaults to all.
#' @param ... Additional arguments forwarded to [New_ModelPar()]
#'   (e.g. `ml_fitting_surv`, `ml_fitting_propensity`).
#' @return A `SimuResults` object.
#' @export
SimuRun <- function(SimuArg, methods,
                    Control = list(grid = 100,
                                   max_iter = 20,
                                   tol = 1e-5,
                                   contraction = 0.5,
                                   eta = 1e-4), sequence = NULL, ...) {
  Validate_method(methods)
  args <- rlang::dots_list(N = SimuArg$initials$N,
                           p = SimuArg$initials$p,
                           max_t = SimuArg$initials$max_t,
                           Control = Control, !!!list(...), .homonyms = "first")
  out <- vector("list", length(methods))
  names(out) <- methods
  templist <- list(Coef = matrix(0, nrow = args$p + 1, ncol = SimuArg$initials$nrep),
                   Var = matrix(0, nrow = args$p + 1, ncol = SimuArg$initials$nrep))
  for (kk in methods) {
    if (kk %in% c("DRIV.s", "DRIV.cf.hz.ml.est")) {
      out[[kk]] <- list(Coef = matrix(0, nrow = args$p + 1, ncol = SimuArg$initials$nrep),
                        Var = rep(0, SimuArg$initials$nrep),
                        Convergence = rep(FALSE, SimuArg$initials$nrep))
      next
    }
    out[[kk]] <- templist
  }
  if (is.null(sequence)) sequence <- 1:SimuArg$initials$nrep

  for (i in sequence) {
    data <- jsonlite::read_json(paste0(SimuArg$Control$save_path,
                                       "DataGenerated/",
                                       SimuArg$initials$N,
                                       SimuArg$Control$Annotation,
                                       "/", i, ".json"), simplifyVector = TRUE)
    data$Covariates <- data$Covariates[, 1:args$p, drop = FALSE]
    colnames(data$Covariates) <- paste0("X", 1:args$p)
    args$dat <- data
    cat("[[rep", i)
    for (kk in methods) {
      args$method <- kk
      ModelPar <- do.call(New_ModelPar, args)
      mod <- DataFitting(ModelPar)
      out[[kk]]$Coef[, i] <- mod$Coef
      cat("\t", out[[kk]]$Coef[1, i])
      if (kk %in% c("DRIV.s", "DRIV.cf.hz.ml.est")) {
        if (is.null(mod$Var)) next
        out[[kk]]$Var[i] <- mod$Var
        out[[kk]]$Convergence[i] <- mod$Convergence
        cat("\t", out[[kk]]$Var[i])
        cat("\t", out[[kk]]$Convergence[i])
        next
      }
      out[[kk]]$Var[, i] <- mod$Var
    }
    cat("]]\n")
  }
  SimuArg$methods <- methods
  SimuArg$SimuResults <- out
  structure(SimuArg, class = "SimuResults")
}


#' Estimate treatment-switching effects on a single dataset
#'
#' Applies each requested method to one observed dataset and returns a `TRTSWE`
#' object holding the per-method estimates.
#'
#' @param dat A list with the observed data (`Covariates`, `Z`, `W`, `T_D_c`,
#'   `event`, and optionally `stime`, `D_status`, `Covariates2`).
#' @param max_t Administrative censoring time.
#' @param methods Character vector of method names (see [Validate_method()]).
#' @param Control List of estimation controls passed to [New_ModelPar()].
#' @param ... Additional arguments forwarded to [New_ModelPar()].
#' @return A `TRTSWE` object.
#' @export
TRTSWE <- function(dat, max_t, methods,
                   Control = list(), ...) {
  N <- nrow(dat$Covariates)
  p <- ncol(dat$Covariates)
  Validate_method(methods)
  args <- rlang::dots_list(N = N, p = p, max_t = max_t,
                           dat = dat, Control = Control,
                           !!!list(...), .homonyms = "first")
  ModelPar <- vector("list", length(methods))
  names(ModelPar) <- methods

  for (kk in methods) {
    args$method <- kk
    ModelPar[[kk]] <- do.call(New_ModelPar, args)
    mod <- DataFitting(ModelPar[[kk]])
    ModelPar[[kk]]$Estim <- mod
  }
  structure(ModelPar, class = "TRTSWE")
}
