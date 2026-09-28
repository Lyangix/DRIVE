# This example is installed with DRIVE. Find it with:
# system.file("examples", "simulation.R", package = "DRIVE")
library(DRIVE)

# Use the bundled dataset immediately.
data(drive_toy)
fit <- TRTSWE(drive_toy, max_t = 5,
              methods = c("ITT", "remove", "recensor", "TimeVar", "DRIVE.joint"))
print(fit)

# Or generate new data from any of the eight original presets.
dat <- SimulateData(N = 200, Scenario = "ii", model = "both", seed = 123)
head(dat$Covariates)

# A small reproducible study; increase N/nrep for a research simulation.
args <- New_SimuScenario(N = 200, nrep = 3)
set.seed(123)
generated <- DataGenerating(args, return_data = TRUE)
print(generated$summary)
results <- SimuRun(args, methods = c("ITT", "recensor"), data = generated$data)
print(results, Comp_parameters = c(0.1, 0.25, 0.25))

# SimuRun can also generate fresh replicates directly, without any files.
set.seed(123)
results <- SimuRun(args, methods = "ITT")

# For larger studies, retain the original disk-backed workflow.
saved_args <- New_SimuScenario(
  N = 200, nrep = 2,
  Control = list(json_save = TRUE, save_path = tempdir(), Annotation = "_example")
)
set.seed(123)
DataGenerating(saved_args)
saved_results <- SimuRun(saved_args, methods = "ITT")
print(saved_results)
