# Default data-generating building blocks for the simulation scenarios.
# Each function is passed by name into New_SimuArg() and called by DataGenerating().

#' Generate unmeasured confounders (independent)
#' @param N Sample size.
#' @param p_U Number of unmeasured confounders.
#' @return An `N x p_U` matrix of uniform draws.
#' @export
unmeasured_Confounding <- function(N, p_U) return(matrix(runif(N * p_U), nrow = N, ncol = p_U))

#' Generate unmeasured confounders dependent on a measured covariate
#' @param N Sample size.
#' @param p Covariate index (column of `Covariates` used as the base).
#' @param p_U Number of unmeasured confounders.
#' @param Covariates Measured covariate matrix.
#' @return A length-`N` vector.
#' @export
unmeasured_Confounding_dependent <- function(N, p, p_U, Covariates) return(Covariates[, p] * runif(N, -1, 1))


#' Initialise measured (and optionally unmeasured) covariates
#' @param N Sample size.
#' @param p Number of measured covariates.
#' @param p_U Number of unmeasured confounders.
#' @param unmeasured_Confounding Function generating unmeasured confounders.
#' @return Covariate matrix.
#' @export
InitCovariates <- function(N, p, p_U, unmeasured_Confounding) {
  L <- matrix(runif(N * (p)), ncol = p, nrow = N)
  args <- as.list(match.call())[-1]
  args <- rlang::dots_list(!!!args,
                           Covariates = L)
  if (p_U >= 1) {
    U <- easy_call(unmeasured_Confounding, args)
    return(cbind(L, U))
  } else {
    return(L)
  }
}


#' Assign treatment via a linear propensity score
#' @param N Sample size.
#' @param p Number of covariates.
#' @param Covariates Covariate matrix.
#' @param gamma Coefficient vector for the propensity score linear predictor.
#' @return Binary treatment vector.
#' @export
InitAssignment <- function(N, p, Covariates, gamma) {
  Z_p <- expit(Covariates[, 1:p, drop = FALSE] %*% gamma[1:p] - mean(Covariates[, 1:p, drop = FALSE] %*% gamma[1:p]))
  return(Z = rbinom(N, 1, Z_p))
}


#' Assign treatment via a nonlinear propensity score
#' @param N Sample size.
#' @param p Number of covariates.
#' @param Covariates Covariate matrix.
#' @return Binary treatment vector.
#' @export
InitAssignment_Nonlinear <- function(N, p, Covariates) {
  L <- Covariates
  L <- cbind(exp(L[, 1] / 2))
  L <- ifelse(L[, 1] > exp(0.25), exp(2) - L[, 1], 0.25 * L[, 1])
  Z_p <- expit(-0.5 * mean(L) + L * 0.5)
  return(Z = rbinom(N, 1, Z_p))
}


#' Survival time under a linear structural model
#' @param N Sample size.
#' @param p Number of covariates.
#' @param Covariates Covariate matrix.
#' @param W Switching-time vector.
#' @param theta True treatment effect.
#' @param alpha Covariate effect vector.
#' @param Z Treatment assignment vector.
#' @return List with `T_D` (observed) and `T_0` (control) event times.
#' @export
SurvTime <- function(N, p, Covariates, W, theta, alpha, Z) {
  W_copy <- W
  W[W <= 0] <- Inf
  W_copy[W_copy <= 0] <- 0
  T <- rexp(N)
  T_0 <- T / (0.25 + Covariates %*% alpha)
  T_D <- T / (0.25 + theta * Z + Covariates %*% alpha)
  T_D_ind <- T_D >= W
  if (any(T_D_ind)) {
    T_D[T_D_ind] <- ((T_D * (0.25 + theta * Z + Covariates %*% alpha) -
                        (theta * Z * W_copy - theta * (1 - Z) * W_copy)) /
                       (0.25 + theta * (1 - Z) + Covariates %*% alpha))[T_D_ind]
  }
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}


#' Survival time under a linear model with endogenous switching
#' @inheritParams SurvTime
#' @param T Pre-drawn baseline event time.
#' @return List with `T_D` and `T_0`.
#' @export
SurvTime_endogenous <- function(N, p, Covariates, W, theta, alpha, Z, T) {
  W_copy <- W
  W[W <= 0] <- Inf
  W_copy[W_copy <= 0] <- 0
  T_0 <- T / (0.25 + Covariates %*% alpha)
  T_D <- T / (0.25 + theta * Z + Covariates %*% alpha)
  T_D_ind <- T_D >= W
  if (any(T_D_ind)) {
    T_D[T_D_ind] <- ((T_D * (0.25 + theta * Z + Covariates %*% alpha) -
                        (theta * Z * W_copy - theta * (1 - Z) * W_copy)) /
                       (0.25 + theta * (1 - Z) + Covariates %*% alpha))[T_D_ind]
  }
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}


#' Survival time under a nonlinear structural model
#' @param N Sample size.
#' @param p Number of covariates (must be 2).
#' @param Covariates Covariate matrix.
#' @param W Switching-time vector.
#' @param theta True treatment effect.
#' @param Z Treatment assignment vector.
#' @return List with `T_D` and `T_0`.
#' @export
SurvTime_Nonlinear <- function(N, p, Covariates, W, theta, Z) {
  if (p != 2) stop("p must be 2")
  L <- Covariates
  L[, 1] <- cbind(exp(L[, 1] / 2))
  L[, 1] <- ifelse(L[, 1] > exp(0.25), exp(3) - L[, 1], 0.25 * L[, 1])
  L[, 2] <- Covariates[, 2] / (1 + exp(Covariates[, 1])) + 1
  Covariates <- cbind(L, Covariates[, p + 1])
  W_copy <- W
  W[W <= 0] <- Inf
  W_copy[W_copy <= 0] <- 0
  T <- rexp(N)
  Covbeta <- 0.1 * abs(Covariates[, 1]) + 0.1 * Covariates[, 2] + 0.1 * Covariates[, p + 1]
  T_0 <- T / (0.1 + Covbeta)
  T_D <- T / (0.1 + theta * Z + Covbeta)
  T_D_ind <- T_D >= W
  if (any(T_D_ind)) {
    T_D[T_D_ind] <- ((T_D * (0.1 + theta * Z + Covbeta) -
                        (theta * Z * W_copy - theta * (1 - Z) * W_copy)) /
                       (0.1 + theta * (1 - Z) + Covbeta))[T_D_ind]
  }
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}


#' Survival time under a nonlinear model with endogenous switching
#' @inheritParams SurvTime_Nonlinear
#' @param T Pre-drawn baseline event time.
#' @return List with `T_D` and `T_0`.
#' @export
SurvTime_endogenous_Nonlinear <- function(N, p, Covariates, W, theta, Z, T) {
  if (p != 2) stop("p must be 2")
  L <- Covariates
  L[, 1] <- cbind(exp(L[, 1] / 2))
  L[, 1] <- ifelse(L[, 1] > exp(0.25), exp(3) - L[, 1], 0.25 * L[, 1])
  L[, 2] <- Covariates[, 2] / (1 + exp(Covariates[, 1])) + 1
  Covariates <- cbind(L, Covariates[, p + 1])
  W_copy <- W
  W[W <= 0] <- Inf
  W_copy[W_copy <= 0] <- 0
  Covbeta <- 0.1 * abs(Covariates[, 1]) + 0.1 * Covariates[, 2] + 0.1 * Covariates[, p + 1]
  T_0 <- T / (0.1 + Covbeta)
  T_D <- T / (0.1 + theta * Z + Covbeta)
  T_D_ind <- T_D >= W
  if (any(T_D_ind)) {
    T_D[T_D_ind] <- ((T_D * (0.1 + theta * Z + Covbeta) -
                        (theta * Z * W_copy - theta * (1 - Z) * W_copy)) /
                       (0.1 + theta * (1 - Z) + Covbeta))[T_D_ind]
  }
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}


#' Treatment-switching time (exogenous)
#' @param N Sample size.
#' @param p Number of covariates.
#' @param Covariates Covariate matrix.
#' @param Z Treatment assignment vector.
#' @param beta Covariate effect vector on switching rate.
#' @param diffcoef Treatment-arm differential for switching rate.
#' @return Switching-time vector.
#' @export
SwitchingTime <- function(N, p, Covariates, Z, beta, diffcoef) {
  W <- ifelse(as.logical(Z), rexp(N) / (0.5 * (0.1 + diffcoef * Z + Covariates %*% beta)),
              rexp(N) / (0.5 * (0.1 - diffcoef * (1 - Z) + Covariates %*% beta)))
  return(W)
}

#' Treatment-switching time (endogenous)
#' @inheritParams SwitchingTime
#' @param T Pre-drawn baseline event time.
#' @param alpha Covariate effect vector on the baseline hazard.
#' @return Switching-time vector.
#' @export
SwitchingTime_endogenous <- function(N, p, Covariates, Z, beta, diffcoef, T, alpha) {
  T_0 <- T / (0.25 + Covariates %*% alpha)
  W <- ifelse(as.logical(Z), rexp(N) / (exp(-T_0) * 2 * (0.1 + diffcoef * Z + Covariates %*% beta)),
              rexp(N) / (exp(-T_0) * 2 * (0.1 - diffcoef * (1 - Z) + Covariates %*% beta)))
  return(W)
}


#' Censoring time (linear)
#' @param N Sample size.
#' @param p Number of covariates.
#' @param Covariates Covariate matrix.
#' @param censoring_par Covariate effect vector on censoring rate.
#' @param censoring_intercept Baseline (intercept) censoring rate.
#' @return Censoring-time vector.
#' @export
CensoringTime <- function(N, p, Covariates, censoring_par, censoring_intercept) {
  C <- rexp(N) / (censoring_intercept + Covariates[, 1:(p), drop = FALSE] %*% censoring_par[1:p])
}

#' Censoring time (nonlinear)
#' @param N Sample size.
#' @param p Number of covariates.
#' @param Covariates Covariate matrix.
#' @return Censoring-time vector.
#' @export
CensoringTime_nonlinear <- function(N, p, Covariates) {
  C <- rexp(N) / (0.1 + Covariates[, 1:(p), drop = FALSE] %*% rep(0.8, p))
}
