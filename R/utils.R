#' Internal helpers
#'
#' @keywords internal
#' @noRd

# Keep only the arguments of `args` that match the formals of `fn`.
arg_filter <- function(args, fn) {
  args[names(formals(fn))]
}

# Call `fn` with the subset of `args` it accepts, reporting missing arguments.
easy_call <- function(fn, args) {
  out <- try(do.call(fn, arg_filter(args, fn)), FALSE)
  if (inherits(out, "try-error")) {
    para <- names(formals(fn))
    stop("These arguments are needed: ", para[!(para %in% names(args))], "\n")
  } else {
    return(out)
  }
}

#' Expit (inverse logit) transform
#'
#' @param d Numeric vector on the linear-predictor scale.
#' @return Numeric vector of probabilities `1 / (1 + exp(-d))`.
#' @export
expit <- function(d) return(1 / (1 + exp(-d)))
