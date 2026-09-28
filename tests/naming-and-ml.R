library(DRIVE)

expect_error <- function(expr, pattern) {
  err <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err), fixed = TRUE))
}

# Canonical names, legacy input names and previously saved S3 objects.
for (pair in list(c("i", "exogenous"), c("ii", "endogenous"))) {
  current <- New_SimuScenario(N = 20, Scenario = pair[1])
  alias <- New_SimuScenario(N = 20, Scenario = pair[2])
  stopifnot(identical(class(current), paste0("SimuArg.", pair[1])),
            identical(class(current), class(alias)),
            identical(current$initials, alias$initials),
            identical(current$parameters, alias$parameters),
            identical(current$Control, alias$Control))
  set.seed(34)
  new_data <- DataGenerating(current, return_data = TRUE)
  class(alias) <- paste0("SimuArg.", pair[2])
  set.seed(34)
  stopifnot(identical(new_data, DataGenerating(alias, return_data = TRUE)))
}
stopifnot(inherits(New_SimuScenario(Scenario = "II"), "SimuArg.ii"))
dat <- SimulateData(N = 80, seed = 17)
for (pair in list(c("DRIVE.joint", "DRIV.s"), c("DRIVE.ML", "DRIV.cf.hz.ml.est"))) {
  model <- New_ModelPar(80, 2, 5, dat, pair[2])
  stopifnot(identical(class(model), paste0("ModelPar.", pair[1])))
}
stopifnot(identical(drive_joint_est_cpp, driv_s_est_cpp),
          identical(drive_ml_est_cpp, driv_cf_ml_est_cpp))
expect_error(Validate_scenario(character()), "Scenario")
expect_error(Validate_method(NA_character_), "methods")
expect_error(TRTSWE(dat, 5, "DRIVE.ML"), "requires ml_fitting_surv")

# Reach the callback through public S3 dispatch, including the legacy name,
# and verify that dat$Covariates2 supplies the propensity predictors.
dat$Covariates2 <- cbind(propensity_only = dat$Covariates[, 1])
check_predictors <- function(data, predictx) {
  stopifnot(identical(names(predictx), "propensity_only"), nrow(data) == 64L)
  stop("propensity predictors verified")
}
for (name in c("DRIVE.ML", "DRIV.cf.hz.ml.est")) {
  expect_error(TRTSWE(dat, 5, name, nfolds = 5,
                      ml_fitting_propensity = check_predictors,
                      ml_fitting_surv = function(data, predictx) NULL),
               "propensity predictors verified")
}

# Exercise the real R cross-fitting wrapper with a spy for the final native
# call. This isolates fold boundaries, argument forwarding and matrix assembly
# from the separate numerical solver tests.
wrapper <- drive_ml_est_cpp
spy_env <- new.env(parent = environment(wrapper))
spy_env$captured <- NULL
spy_env$driv_cf_ml_est <- function(init_parameters, ...) {
  spy_env$captured <- c(list(init_parameters = init_parameters), list(...))
  list(x = 0.2, var = 0.04, Convergence = TRUE)
}
environment(wrapper) <- spy_env
N <- 8L
grid <- c(0.5, 1, 2)
X <- cbind(subject = seq_len(N), x = seq_len(N) / N)
X2 <- cbind(subject = seq_len(N), other = 2 * seq_len(N))
IV <- rep(0:1, N / 2)
status <- matrix(IV, nrow = N, ncol = length(grid))
seen <- integer()
propensity <- function(data, predictx) {
  stopifnot(!"IV" %in% names(predictx),
            identical(names(predictx), c("subject", "other")),
            !any(predictx$subject %in% data$subject),
            nrow(data) + nrow(predictx) == N)
  predictx$subject / (N + 1)
}
hazard <- function(data, predictx) {
  stopifnot(!any(predictx$subject %in% data$Covariates[, "subject"]),
            is.matrix(data$D_status),
            identical(ncol(data$D_status), length(grid)),
            identical(data$stime, grid))
  seen <<- c(seen, predictx$subject)
  outer(predictx$subject / 100, c(grid[1], diff(grid)))
}
args <- list(init_parameters = 0.07, time = rep(2, N), event = rep(TRUE, N),
             IV = IV, Covariates = X, Covariates2 = X2,
             D_status = status, stime = grid,
             ml_fitting_propensity = propensity, ml_fitting_surv = hazard,
             nfolds = 4, seed = 42)
do.call(wrapper, args)
stopifnot(identical(as.integer(sort(seen)), seq_len(N)),
          identical(spy_env$captured$init_parameters, 0.07),
          isTRUE(all.equal(spy_env$captured$IV_c, IV - seq_len(N) / (N + 1))),
          isTRUE(all.equal(spy_env$captured$ConfoundingPart,
                            outer(seq_len(N) / 100, c(grid[1], diff(grid))))))

# Leave-one-out folds and a single grid point must retain matrix dimensions.
one <- args
one$nfolds <- N
one$stime <- 1
one$D_status <- matrix(IV, nrow = N, ncol = 1)
one$ml_fitting_surv <- function(data, predictx) {
  stopifnot(identical(dim(data$D_status), c(N - 1L, 1L)))
  matrix(0.1, nrow = nrow(predictx), ncol = 1)
}
do.call(wrapper, one)
for (folds in list(0, 1, 2.5, N + 1L, NA_real_)) {
  bad <- args
  bad$nfolds <- folds
  expect_error(do.call(wrapper, bad), "nfolds")
}
bad <- args
bad$ml_fitting_propensity <- function(data, predictx) 0.5
expect_error(do.call(wrapper, bad), "probability vector")
bad$ml_fitting_propensity <- function(data, predictx) rep(NA_real_, nrow(predictx))
expect_error(do.call(wrapper, bad), "probability vector")
bad <- args
bad$ml_fitting_surv <- function(data, predictx) list(surv = matrix(0, 2, 3))
expect_error(do.call(wrapper, bad), "matrix of hazard increments")
bad$ml_fitting_surv <- function(data, predictx) matrix(0, nrow(predictx), 1)
expect_error(do.call(wrapper, bad), "matrix of hazard increments")

# ML-only print used to create bogus alpha1/alpha0 column names for a scalar.
fit <- structure(list(DRIVE.ML = list(Estim = list(Coef = 0.2, Var = 0.04))),
                  class = "TRTSWE")
printed <- capture.output(print(fit, all = TRUE))
stopifnot(any(grepl("theta", printed)), !any(grepl("alpha", printed)))
names(fit) <- "DRIV.cf.hz.ml.est"
stopifnot(identical(printed, capture.output(print(fit, all = TRUE))))

# Correct cumulative-hazard grid alignment, including time before the first
# forest event, between events, and beyond the last training event.
source(system.file("examples", "drive_ml_learners.R", package = "DRIVE"))
H <- rbind(c(0.2, 0.5, 0.9), c(0.1, 0.4, 0.8))
increments <- drive_grid_increments(H, c(1, 2, 4), c(0.5, 1, 1.5, 3, 5))
stopifnot(isTRUE(all.equal(increments,
  rbind(c(0, 0.2, 0, 0.3, 0.4), c(0, 0.1, 0, 0.3, 0.4)))))
stopifnot(identical(dim(drive_grid_increments(H[1, , drop = FALSE], c(1, 2, 4), 1)),
                    c(1L, 1L)))
