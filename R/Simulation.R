#' Construct one of the original simulation scenarios
#'
#' Presets from the original `SimuArgs.R`, with smaller defaults for examples.
#' Each preset uses two measured uniform covariates and one independent
#' unmeasured uniform confounder. Customize other designs with [New_SimuArg()].
#'
#' @param N Number of subjects per replicate.
#' @param nrep Number of replicates.
#' @param Scenario `"i"` for scenario i (exogenous switching), or `"ii"` for
#'   scenario ii (endogenous switching). The former names remain accepted aliases.
#' @param model Nonlinear component: `"linear"` (neither), `"propensity"`,
#'   `"survival"`, or `"both"`.
#' @param max_t Administrative follow-up limit.
#' @param theta Treatment effect on the additive hazard.
#' @param Control Generation controls passed to [New_SimuArg()]. JSON saving
#'   is disabled by default. The default `Annotation` distinguishes the eight
#'   presets when saving; an explicitly supplied annotation takes precedence.
#' @details
#' The original study used `nrep = 1000` and `N = 1600` or `3200`.
#' These presets retain its coefficients: `gamma = c(1, -1)`,
#' `alpha = beta = rep(0.25, 3)`, `diffcoef = 0.5`,
#' `censoring_par = c(0.01, 0.01)`, and `censoring_intercept = 0.1`.
#'
#' Endogenous switching depends on the same exponential draw used to generate
#' survival. For nonlinear survival, the switching mechanism still uses the
#' original linear latent reference time. Nonpositive switching times from the
#' original formula denote no switch and are stored as `Inf`.
#' @return A simulation argument object accepted by [DataGenerating()] and
#'   [SimuRun()].
#' @examples
#' args <- New_SimuScenario(N = 200, nrep = 2, Scenario = "ii",
#'                         model = "both")
#' set.seed(123)
#' generated <- DataGenerating(args, return_data = TRUE)
#' generated$summary
#' @md
#' @export
New_SimuScenario <- function(N = 200, nrep = 1,
                             Scenario = "i",
                             model = c("linear", "propensity", "survival", "both"),
                             max_t = 5, theta = 0.1, Control = list()) {
  Scenario <- normalize_scenario(Scenario)
  model <- match.arg(model)
  check_simulation_count(N, "N")
  check_simulation_count(nrep, "nrep")
  if (!is.numeric(max_t) || length(max_t) != 1L ||
      !is.finite(max_t) || max_t <= 0) {
    stop("max_t must be a positive finite number.")
  }
  nonlinear_survival <- model %in% c("survival", "both")
  baseline <- if (nonlinear_survival) 0.1 else 0.25
  if (!is.numeric(theta) || length(theta) != 1L ||
      !is.finite(theta) || theta <= -baseline) {
    stop("theta must be finite and greater than ", -baseline,
         " to keep the preset hazards positive.")
  }
  survival_function <- if (Scenario == "ii") {
    if (nonlinear_survival) SurvTime_endogenous_Nonlinear else SurvTime_endogenous
  } else {
    if (nonlinear_survival) SurvTime_Nonlinear else SurvTime
  }
  suffix <- c(linear = "", propensity = "prop", survival = "surv", both = "both")[[model]]
  annotation <- if (Scenario == "ii") {
    if (nzchar(suffix)) paste0("end_", suffix) else "end"
  } else suffix
  Control <- rlang::dots_list(!!!Control, Annotation = annotation, .homonyms = "first")
  New_SimuArg(
    nrep = nrep, N = N, p = 2, p_U = 1, Scenario = Scenario,
    max_t = max_t, theta = theta,
    unmeasured_Confounding = unmeasured_Confounding,
    InitCovariates = InitCovariates,
    InitAssignment = if (model %in% c("propensity", "both")) {
      InitAssignment_Nonlinear
    } else InitAssignment,
    SurvTime = survival_function,
    SwitchingTime = if (Scenario == "ii") {
      SwitchingTime_endogenous
    } else SwitchingTime,
    CensoringTime = CensoringTime, Control = Control,
    gamma = c(1, -1), alpha = rep(0.25, 3), beta = rep(0.25, 3),
    diffcoef = 0.5, censoring_par = c(0.01, 0.01),
    censoring_intercept = 0.1
  )
}

#' Simulate one ready-to-fit treatment-switching dataset
#'
#' Generate observed data using the original simulation presets, without
#' creating files. Only measured covariates are included in `Covariates`.
#'
#' @inheritParams New_SimuScenario
#' @param seed Optional nonnegative integer seed. If supplied, the caller's
#'   random-number state is restored on exit. If `NULL`, use the current state.
#' @return A list containing `Covariates` (an `N` by 2 matrix), initial treatment
#'   `Z`, switching time `W` (`Inf` means no switch), observed time `T_D_c`, and
#'   logical event indicator `event`. Pass it directly to [TRTSWE()].
#' @examples
#' dat <- SimulateData(N = 200, seed = 2025)
#' head(dat$Covariates)
#' fit <- TRTSWE(dat, max_t = 5, methods = c("ITT", "recensor"))
#' fit
#' @md
#' @export
SimulateData <- function(N = 200, Scenario = "i",
                         model = c("linear", "propensity", "survival", "both"),
                         max_t = 5, theta = 0.1, seed = NULL) {
  args <- New_SimuScenario(N = N, Scenario = normalize_scenario(Scenario),
                          model = match.arg(model), max_t = max_t, theta = theta)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
      stop("seed must be NULL or a nonnegative integer.")
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }
  dat <- DataGenerating(args, return_data = TRUE)$data[[1L]]
  dat[c("Covariates", "Z", "W", "T_D_c", "event")]
}

check_simulation_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop(name, " must be a positive integer.")
  }
}

simulation_directory <- function(SimuArg) {
  root <- SimuArg$Control$save_path
  if (identical(root, "")) root <- "."
  file.path(root, "DataGenerated",
            paste0(SimuArg$initials$N, SimuArg$Control$Annotation))
}

# Preserve the legacy JSON layout, but never expose U as a fitting covariate.
observed_simulation_data <- function(dat, p) {
  if (ncol(dat$Covariates) > p) {
    dat$U <- dat$Covariates[, seq.int(p + 1L, ncol(dat$Covariates)), drop = FALSE]
  }
  dat$Covariates <- dat$Covariates[, seq_len(p), drop = FALSE]
  colnames(dat$Covariates) <- paste0("X", seq_len(p))
  dat
}

generate_simulation <- function(SimuArg, return_data, endogenous) {
  check_simulation_count(SimuArg$initials$nrep, "nrep")
  check_simulation_count(SimuArg$initials$N, "N")
  if (!is.logical(return_data) || length(return_data) != 1L || is.na(return_data)) {
    stop("return_data must be TRUE or FALSE.")
  }
  args <- rlang::dots_list(!!!SimuArg$parameters, !!!SimuArg$initials)
  rate_names <- c("Z_proportion", "switching_rate_overall", "switching_rate_from_0",
                  "switching_rate_from_1", "censoring_rate", "adcensoring_rate")
  rates <- matrix(NA_real_, nrow = args$nrep, ncol = length(rate_names),
                  dimnames = list(NULL, rate_names))
  datasets <- if (return_data) vector("list", args$nrep) else NULL
  negative_count <- 0L
  if (isTRUE(SimuArg$Control$json_save)) {
    directory <- simulation_directory(SimuArg)
    if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
      stop("Could not create simulation directory: ", directory)
    }
  }
  for (kk in seq_len(args$nrep)) {
    Covariates <- easy_call(SimuArg$InitCovariates,
      c(args, list(unmeasured_Confounding = SimuArg$unmeasured_Confounding)))
    Z <- easy_call(SimuArg$InitAssignment, c(args, list(Covariates = Covariates)))
    switching_args <- c(args, list(Covariates = Covariates, Z = Z))
    if (endogenous) switching_args$T <- rexp(args$N)
    W <- easy_call(SimuArg$SwitchingTime, switching_args)
    T <- easy_call(SimuArg$SurvTime, c(switching_args, list(W = W)))
    C <- as.vector(easy_call(SimuArg$CensoringTime,
                             c(args, list(Covariates = Covariates))))
    T_D <- as.vector(T$T_D)
    T_D_c <- pmin(T_D, C, args$max_t)
    event <- T_D <= C & T_D <= args$max_t
    negative_count <- negative_count + sum(W < 0)
    W[W <= 0] <- Inf
    if (any(T_D <= 0)) warning("Nonpositive time-to-event outcome occurs")
    dat <- list(Covariates = Covariates, Z = Z, W = as.vector(W),
                T_D_c = T_D_c, T_D = T_D, T_0 = as.vector(T$T_0),
                C = C, event = event)
    if (isTRUE(SimuArg$Control$json_save)) {
      # String encoding preserves Inf; full precision avoids introducing ties.
      jsonlite::write_json(dat, file.path(directory, paste0(kk, ".json")),
                            digits = NA)
    }
    if (return_data) datasets[[kk]] <- observed_simulation_data(dat, args$p)
    switched <- T_D_c > W
    rates[kk, ] <- c(mean(Z), mean(switched), mean(switched[Z == 0]),
                     mean(switched[Z == 1]), mean(!event), mean(T_D > args$max_t))
  }
  summary <- as.list(colMeans(rates))
  # Keep the original spelling for compatibility; now counts all replicates.
  summary$anyNegetive <- negative_count
  if (return_data) list(summary = summary, data = datasets) else summary
}
