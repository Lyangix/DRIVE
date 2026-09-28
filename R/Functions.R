#' Run a simulation study
#'
#' Generate replicates in memory, use supplied datasets, or read previously
#' saved JSON datasets, then fit each requested method,
#' collecting coefficients and variances into a `SimuResults` object.
#'
#' @param SimuArg A `SimuArg` object from [New_SimuArg()].
#' @param methods Character vector of method names (see [Validate_method()]).
#' @param Control List of estimation controls passed to [New_ModelPar()].
#' @param sequence Optional integer vector of replicate indices; defaults to all.
#'   Unselected replicates are stored as `NA`.
#' @param data Optional list of replicate datasets, such as the `data` component
#'   returned by `DataGenerating(SimuArg, return_data = TRUE)`. Must have `nrep`
#'   entries. When omitted, read JSON files if `Control$json_save` in `SimuArg`
#'   is `TRUE`; otherwise generate replicates in memory.
#' @param ... Additional arguments forwarded to [New_ModelPar()]
#'   (e.g. `ml_fitting_surv`, `ml_fitting_propensity`, `nfolds`, `seed`).
#'   See [DRIVE_ML] for the ML callback interface.
#' @return A `SimuResults` object.
#' @md
#' @export
SimuRun <- function(SimuArg, methods,
                    Control = list(grid = 100,
                                   max_iter = 20,
                                   tol = 1e-5,
                                   contraction = 0.5,
                                   eta = 1e-4), sequence = NULL, ..., data = NULL) {
  methods <- unique(normalize_methods(methods))
  check_simulation_count(SimuArg$initials$nrep, "nrep")
  nrep <- SimuArg$initials$nrep
  if (!length(methods)) stop("Specify at least one estimation method.")
  if (is.null(sequence)) sequence <- seq_len(nrep)
  if (!is.numeric(sequence) || !length(sequence) || anyNA(sequence) ||
      any(!is.finite(sequence)) || any(sequence != floor(sequence)) ||
      any(sequence < 1 | sequence > nrep) || anyDuplicated(sequence)) {
    stop("sequence must contain unique replicate indices between 1 and nrep.")
  }
  if (is.null(data) && !isTRUE(SimuArg$Control$json_save)) {
    data <- DataGenerating(SimuArg, return_data = TRUE)$data
  }
  if (!is.null(data) && (!is.list(data) || length(data) != nrep)) {
    stop("data must be a list with one dataset per replicate (nrep entries).")
  }
  args <- rlang::dots_list(N = SimuArg$initials$N,
                           p = SimuArg$initials$p,
                           max_t = SimuArg$initials$max_t,
                           Control = Control, !!!list(...), .homonyms = "first")
  out <- vector("list", length(methods))
  names(out) <- methods
  templist <- list(Coef = matrix(NA_real_, nrow = args$p + 1, ncol = SimuArg$initials$nrep),
                   Var = matrix(NA_real_, nrow = args$p + 1, ncol = SimuArg$initials$nrep))
  for (kk in methods) {
    if (kk %in% c("DRIVE.joint", "DRIVE.ML")) {
      out[[kk]] <- list(Coef = matrix(NA_real_, nrow = args$p + 1, ncol = SimuArg$initials$nrep),
                        Var = rep(NA_real_, SimuArg$initials$nrep),
                        Convergence = rep(NA, SimuArg$initials$nrep))
      next
    }
    out[[kk]] <- templist
  }
  for (i in sequence) {
    if (is.null(data)) {
      dat <- jsonlite::read_json(
        file.path(simulation_directory(SimuArg), paste0(i, ".json")),
        simplifyVector = TRUE)
      # jsonlite represents infinite switching times as strings.
      dat$W <- as.numeric(dat$W)
    } else {
      dat <- data[[i]]
    }
    args$dat <- observed_simulation_data(dat, args$p)
    cat("[[rep", i)
    for (kk in methods) {
      args$method <- kk
      ModelPar <- do.call(New_ModelPar, args)
      mod <- DataFitting(ModelPar)
      # Cross-fitted DRIVE ML estimates theta only; do not recycle it into the
      # nuisance-coefficient rows used by the other estimators.
      out[[kk]]$Coef[seq_along(mod$Coef), i] <- mod$Coef
      cat("\t", out[[kk]]$Coef[1, i])
      if (kk %in% c("DRIVE.joint", "DRIVE.ML")) {
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
  SimuArg$sequence <- sequence
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
#' @param ... Additional arguments forwarded to [New_ModelPar()]. For
#'   `DRIVE.ML`, supply `ml_fitting_surv`, `ml_fitting_propensity`, and optionally
#'   `nfolds` and `seed`; see [DRIVE_ML] for the callback interface.
#' @return A `TRTSWE` object.
#' @seealso [DRIVE_ML]
#' @md
#' @export
TRTSWE <- function(dat, max_t, methods,
                   Control = list(), ...) {
  N <- nrow(dat$Covariates)
  p <- ncol(dat$Covariates)
  methods <- unique(normalize_methods(methods))
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
