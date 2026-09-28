library(DRIVE)

expect_error <- function(expr) {
  stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
}

# All eight original combinations generate valid observed follow-up and keep
# the latent confounder out of fitting covariates.
for (scenario in c("i", "ii")) {
  for (model in c("linear", "propensity", "survival", "both")) {
    args <- New_SimuScenario(N = 100, nrep = 2, Scenario = scenario, model = model)
    set.seed(2025)
    generated <- DataGenerating(args, return_data = TRUE)
    stopifnot(length(generated$data) == 2L)
    for (dat in generated$data) {
      stopifnot(identical(dim(dat$Covariates), c(100L, 2L)),
                identical(dim(dat$U), c(100L, 1L)),
                identical(colnames(dat$Covariates), c("X1", "X2")),
                all(dat$Z %in% 0:1), is.logical(dat$event),
                all(dat$W > 0), all(is.finite(dat$T_D)),
                all(dat$T_D_c > 0 & dat$T_D_c <= 5),
                identical(dat$T_D_c, pmin(dat$T_D, dat$C, 5)),
                identical(dat$event, dat$T_D <= dat$C & dat$T_D <= 5))
      if (model %in% c("linear", "propensity")) {
        baseline <- 0.25 + 0.25 * rowSums(cbind(dat$Covariates, dat$U))
        treated_time <- dat$Z * pmin(dat$T_D, dat$W) +
          (1 - dat$Z) * pmax(dat$T_D - dat$W, 0)
        stopifnot(isTRUE(all.equal(baseline * dat$T_D + 0.1 * treated_time,
                                   baseline * dat$T_0)))
      }
    }
    set.seed(2025)
    stopifnot(identical(DataGenerating(args), generated$summary))
  }
}

# Seeded convenience calls reproduce the data without disturbing caller RNG.
set.seed(12)
before <- .Random.seed
dat <- SimulateData(seed = 2025)
stopifnot(identical(before, .Random.seed),
          identical(dat, SimulateData(seed = 2025)),
          !identical(dat, SimulateData(seed = 2026)))
rm(.Random.seed, envir = .GlobalEnv)
invisible(SimulateData(N = 10, seed = 12))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
invisible(SimulateData(N = 10))
stopifnot(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

# The checked-in fixture is exactly reproducible and directly fits all of the
# estimators that do not require user-supplied machine-learning functions.
data(drive_toy)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
stopifnot(identical(drive_toy, SimulateData(N = 200, seed = 2025)))
set.seed(123)
fits <- TRTSWE(drive_toy, max_t = 5,
               methods = c("ITT", "remove", "recensor", "TimeVar", "DRIVE.joint"))
for (fit in fits) {
  stopifnot(length(fit$Estim$Coef) == 3L, all(is.finite(fit$Estim$Coef)))
}

# In-memory, explicit-data and JSON workflows fit identical data. A nested
# save_path needs no trailing slash, and Inf must survive the JSON round trip.
args <- New_SimuScenario(N = 200, nrep = 3)
set.seed(42)
generated <- DataGenerating(args, return_data = TRUE)
set.seed(100)
explicit <- SimuRun(args, methods = "ITT", data = generated$data, sequence = c(1, 3))
stopifnot(all(is.na(explicit$SimuResults$ITT$Coef[, 2])),
          all(is.na(explicit$SimuResults$ITT$Var[, 2])))
set.seed(42)
automatic <- SimuRun(args, methods = "ITT")
set.seed(42)
generated <- DataGenerating(args, return_data = TRUE)
manual <- SimuRun(args, methods = "ITT", data = generated$data)
stopifnot(identical(automatic$SimuResults, manual$SimuResults))

saved_args <- New_SimuScenario(N = 200, nrep = 3,
  Control = list(json_save = TRUE,
                 save_path = file.path(tempdir(), "drive-tests", "nested"),
                 Annotation = "_roundtrip"))
set.seed(42)
saved <- DataGenerating(saved_args, return_data = TRUE)
stopifnot(identical(generated$data, saved$data),
          any(is.infinite(saved$data[[1]]$W)))
set.seed(100)
disk <- SimuRun(saved_args, methods = "ITT", sequence = c(1, 3))
stopifnot(isTRUE(all.equal(explicit$SimuResults, disk$SimuResults, tolerance = 1e-12)))

# The default empty save_path means the working directory; presets use distinct
# annotations unless the caller explicitly overrides one.
previous <- getwd()
scratch <- tempfile("drive-relative-path-")
dir.create(scratch)
setwd(scratch)
tryCatch(DataGenerating(New_SimuScenario(N = 5, Control = list(json_save = TRUE))),
         finally = setwd(previous))
stopifnot(file.exists(file.path(scratch, "DataGenerated", "5", "1.json")))
annotations <- c()
for (scenario in c("i", "ii")) {
  for (model in c("linear", "propensity", "survival", "both")) {
    annotations <- c(annotations,
      New_SimuScenario(Scenario = scenario, model = model)$Control$Annotation)
  }
}
stopifnot(!anyDuplicated(annotations),
          New_SimuScenario(Scenario = "ii",
            Control = list(Annotation = "custom"))$Control$Annotation == "custom")

# Printing a subset has the same summary as an equivalent compact result.
compact <- explicit
compact$initials$nrep <- 2
compact$sequence <- 1:2
compact$SimuResults$ITT$Coef <- explicit$SimuResults$ITT$Coef[, c(1, 3)]
compact$SimuResults$ITT$Var <- explicit$SimuResults$ITT$Var[, c(1, 3)]
stopifnot(identical(capture.output(print(compact)), capture.output(print(explicit))))

# Validate public preset inputs and replicate selection before fitting.
for (n in list(0, -1, 1.5, NA_real_, Inf, numeric(), "100")) {
  expect_error(New_SimuScenario(N = n))
  expect_error(New_SimuScenario(nrep = n))
}
expect_error(New_SimuScenario(Scenario = "unknown"))
expect_error(New_SimuScenario(model = "unknown"))
expect_error(New_SimuScenario(max_t = 0))
expect_error(New_SimuScenario(theta = -0.25))
expect_error(SimulateData(seed = -1))
expect_error(DataGenerating(args, return_data = NA))
expect_error(SimuRun(args, "ITT", sequence = c(1, 1)))
expect_error(SimuRun(args, "ITT", sequence = 4))
expect_error(SimuRun(args, "ITT", data = list(dat)))

# Degenerate no-switch/no-event custom generators retain matrix dimensions.
custom <- New_SimuScenario(N = 1, nrep = 2)
custom$SwitchingTime <- function(N) rep(-1, N)
custom$SurvTime <- function(N) list(T_D = rep(10, N), T_0 = rep(10, N))
custom$CensoringTime <- function(N) rep(20, N)
generated <- DataGenerating(custom, return_data = TRUE)
stopifnot(generated$summary$anyNegetive == 2,
          identical(dim(generated$data[[1]]$Covariates), c(1L, 2L)),
          identical(generated$data[[1]]$event, FALSE),
          generated$data[[1]]$W == Inf)
