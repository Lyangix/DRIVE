source("ClassBuilder.R")
source("Generic.R")
source("Functions.R")
source("Valid.R")
source("Function/integral2_est.R")
source("Function/integral2.R")
source("Function/integral.R")
source("Function/cf_integral.R")
source("Function/cf_integral_est.R")
source("Function/integral_customized.R")
source("Function/cf_ml_integral_est.R")
source("Function/cf_ml_integral_est_rateCal_t.R")
source("Function/Strcf_ml_integral_est.R")
library(RcppArmadillo)


# Unmeasured Confounding --------------------------------------------------

unmeasured_Confounding = function(N, p_U) return(matrix(runif(N*p_U), nrow = N, ncol = p_U))
unmeasured_Confounding_dependent = function(N, p, p_U, Covariates) return(Covariates[, p]*runif(N, -1, 1))


# Covariates --------------------------------------------------------------

InitCovariates = function(N, p, p_U, unmeasured_Confounding){
  L = matrix(runif(N*(p)), ncol = p, nrow = N)
  args = as.list(match.call())[-1]
  args = rlang::dots_list(!!!args,
                          Covariates = L)
  if(p_U >= 1){
    U = easy_call(unmeasured_Confounding, args)
    return(cbind(L, U))
  } else {
    return(L)
  }
}



# Propensity Score --------------------------------------------------------

InitAssignment = function(N, p, Covariates, gamma){
  Z_p = expit(Covariates[, 1:p, drop = FALSE] %*% gamma[1:p] - mean(Covariates[, 1:p, drop = FALSE] %*% gamma[1:p]))
  return(Z = rbinom(N, 1, Z_p))
}


InitAssignment_Nonlinear = function(N, p, Covariates) {
  # if(p != 2) stop("p must be 2")
  L = Covariates
  L = cbind(exp(L[, 1]/2))
  L = ifelse(L[, 1] > exp(0.25), exp(2) - L[, 1], 0.25 * L[, 1])
  Z_p = expit(-0.5*mean(L)+ L*0.5)
  return(Z = rbinom(N, 1, Z_p))
}


# Covariates = matrix(runif(1000), nrow = 1000)
# plot(Covariates[, 1], L)
# plot(Covariates[, 2], L[, 2])


# Survival Model ----------------------------------------------------------

SurvTime = function(N, p, Covariates, W, theta, alpha, Z){
  W_copy = W
  W[W <= 0] = Inf
  W_copy[W_copy <= 0] = 0
  T = rexp(N)
  T_0 = T/(0.25 + Covariates %*% alpha)
  T_D = T/(0.25 + theta*Z + Covariates %*% alpha)
  T_D_ind = T_D >= W
  if(any(T_D_ind)){
    T_D[T_D_ind] = ((T_D*(0.25 + theta*Z + Covariates %*% alpha)- 
                       (theta * Z * W_copy - theta* (1-Z)*W_copy))/
                      (0.25 + theta*(1-Z) + Covariates %*% alpha))[T_D_ind]
  }
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}


SurvTime_endogenous = function(N, p, Covariates, W, theta, alpha, Z, T){
  W_copy = W
  W[W <= 0] = Inf
  W_copy[W_copy <= 0] = 0
  T_0 = T/(0.25 + Covariates %*% alpha)
  T_D = T/(0.25 + theta*Z + Covariates %*% alpha)
  T_D_ind = T_D >= W
  if(any(T_D_ind)){
    T_D[T_D_ind] = ((T_D*(0.25 + theta*Z + Covariates %*% alpha)- 
                       (theta * Z * W_copy - theta* (1-Z)*W_copy))/
                      (0.25 + theta*(1-Z) + Covariates %*% alpha))[T_D_ind]
  }
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}


SurvTime_Nonlinear = function(N, p, Covariates, W, theta, Z){
  if(p != 2) stop("p must be 2")
  L = Covariates
  L[, 1] = cbind(exp(L[, 1]/2))
  L[, 1] = ifelse(L[, 1] > exp(0.25), exp(3) - L[, 1], 0.25 * L[, 1])
  L[, 2] = Covariates[, 2]/(1 + exp(Covariates[, 1])) + 1
  Covariates = cbind(L, Covariates[, p+1])
  W_copy = W
  W[W <= 0] = Inf
  W_copy[W_copy <= 0] = 0
  T = rexp(N)
  Covbeta = 0.1 * abs(Covariates[, 1]) + 0.1*Covariates[, 2] + 0.1 * Covariates[, p+1]
  T_0 = T/(0.1 + Covbeta)
  T_D = T/(0.1 + theta*Z + Covbeta)
  T_D_ind = T_D >= W
  if(any(T_D_ind)){
    T_D[T_D_ind] = ((T_D*(0.1 + theta*Z + Covbeta)- 
                       (theta * Z * W_copy - theta* (1-Z)*W_copy))/
                      (0.1 + theta*(1-Z) + Covbeta))[T_D_ind]
  }
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}


SurvTime_endogenous_Nonlinear = function(N, p, Covariates, W, theta, Z, T){
  if(p != 2) stop("p must be 2")
  L = Covariates
  L[, 1] = cbind(exp(L[, 1]/2))
  L[, 1] = ifelse(L[, 1] > exp(0.25), exp(3) - L[, 1], 0.25 * L[, 1])
  L[, 2] = Covariates[, 2]/(1 + exp(Covariates[, 1])) + 1
  Covariates = cbind(L, Covariates[, p+1])
  W_copy = W
  W[W <= 0] = Inf
  W_copy[W_copy <= 0] = 0
  # T = rexp(N)
  Covbeta = 0.1 * abs(Covariates[, 1]) + 0.1*Covariates[, 2] + 0.1 * Covariates[, p+1]
  T_0 = T/(0.1 + Covbeta)
  T_D = T/(0.1 + theta*Z + Covbeta)
  T_D_ind = T_D >= W
  if(any(T_D_ind)){
    T_D[T_D_ind] = ((T_D*(0.1 + theta*Z + Covbeta)- 
                       (theta * Z * W_copy - theta* (1-Z)*W_copy))/
                      (0.1 + theta*(1-Z) + Covbeta))[T_D_ind]
  }
  return(list(T_D = as.vector(T_D),
              T_0 = T_0))
}


# Censoring and Switching -------------------------------------------------

SwitchingTime = function(N, p, Covariates, Z, beta, diffcoef){
  W = ifelse(as.logical(Z), rexp(N)/(0.5*(0.1 + diffcoef*Z + Covariates %*% beta)),
             rexp(N)/(0.5*(0.1 - diffcoef*(1 - Z) + Covariates %*% beta)))
  return(W)
}

SwitchingTime_endogenous = function(N, p, Covariates, Z, beta, diffcoef, T, alpha) {
  T_0 = T/(0.25 + Covariates %*% alpha)
  W = ifelse(as.logical(Z), rexp(N)/(exp(-T_0)*2*(0.1 + diffcoef*Z + Covariates %*% beta)),
             rexp(N)/(exp(-T_0)*2*(0.1 - diffcoef*(1 - Z) + Covariates %*% beta)))
  return(W)
}


CensoringTime = function(N, p, Covariates, censoring_par, censoring_intercept) {
  C = rexp(N)/(censoring_intercept + Covariates[, 1:(p), drop = FALSE] %*% censoring_par[1:p])
}

CensoringTime_nonlinear = function(N, p, Covariates) {
  C = rexp(N)/(0.1 + Covariates[, 1:(p), drop = FALSE] %*% rep(0.8, p))
}


  
