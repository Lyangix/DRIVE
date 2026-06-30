#' Validate a simulation scenario name
#'
#' @param x Scenario label; must be one of `"exogenous"` or `"endogenous"`.
#' @return Invisibly `NULL`; called for its side effect of stopping on invalid input.
#' @export
Validate_scenario <- function(x) {
  values <- unclass(x)
  stopifnot(all(values %in% c("exogenous", "endogenous")))
  stopifnot(length(x) <= 1)
}


#' Validate estimation method name(s)
#'
#' @param x Character vector of method names. Supported methods are
#'   `"ITT"`, `"remove"`, `"recensor"`, `"TimeVar"`, `"DRIV.s"` and
#'   `"DRIV.cf.hz.ml.est"`.
#' @return Invisibly `NULL`; called for its side effect of stopping on invalid input.
#' @export
Validate_method <- function(x) {
  values <- unclass(x)
  stopifnot(all(values %in% c("ITT", "remove", "recensor", "TimeVar",
                              "DRIV.s", "DRIV.cf.hz.ml.est")))
}
