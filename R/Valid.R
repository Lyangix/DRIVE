#' Validate a simulation scenario name
#'
#' @param x Scenario label: `"i"` (exogenous switching) or `"ii"` (endogenous
#'   switching). The old names `"exogenous"` and `"endogenous"` are aliases.
#' @return Invisibly `NULL`; called for its side effect of stopping on invalid input.
#' @export
Validate_scenario <- function(x) {
  normalize_scenario(x)
  invisible(NULL)
}


#' Validate estimation method name(s)
#'
#' @param x Character vector of method names. Supported methods are
#'   `"ITT"`, `"remove"`, `"recensor"`, `"TimeVar"`, `"DRIVE.joint"` and
#'   `"DRIVE.ML"`. `"DRIV.s"` and `"DRIV.cf.hz.ml.est"` remain accepted aliases
#'   for the joint and ML estimators, respectively.
#' @return Invisibly `NULL`; called for its side effect of stopping on invalid input.
#' @export
Validate_method <- function(x) {
  normalize_methods(x)
  invisible(NULL)
}

normalize_scenario <- function(x) {
  aliases <- c(i = "i", ii = "ii", exogenous = "i", endogenous = "ii")
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !tolower(x) %in% names(aliases)) {
    stop('Scenario must be "i" or "ii" (aliases: "exogenous", "endogenous").')
  }
  unname(aliases[tolower(x)])
}

normalize_methods <- function(x) {
  aliases <- c(ITT = "ITT", remove = "remove", recensor = "recensor",
               TimeVar = "TimeVar", DRIVE.joint = "DRIVE.joint", DRIVE.ML = "DRIVE.ML",
               DRIV.s = "DRIVE.joint", DRIV.cf.hz.ml.est = "DRIVE.ML")
  if (!is.character(x) || !length(x) || anyNA(x) || any(!x %in% names(aliases))) {
    stop('methods must use ITT, remove, recensor, TimeVar, DRIVE.joint or DRIVE.ML.')
  }
  unname(aliases[x])
}
