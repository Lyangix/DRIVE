# A complete application of cross-fitted DRIVE ML to the bundled toy data.
# install.packages(c("rpart", "randomForestSRC"))
library(DRIVE)
source(system.file("examples", "drive_ml_learners.R", package = "DRIVE"))
data(drive_toy)

set.seed(123)
fit <- TRTSWE(
  drive_toy, max_t = 5, methods = c("DRIVE.joint", "DRIVE.ML"),
  ml_fitting_propensity = ml_propensity_tree,
  ml_fitting_surv = ml_survival_forest,
  nfolds = 5, seed = 123,
  Control = list(max_iter = 50, tol = 1e-5)
)
print(fit)
fit$DRIVE.ML$Estim$Convergence

# The same callbacks work with repeated simulations. Uncomment to run.
# args <- New_SimuScenario(N = 400, nrep = 3, Scenario = "i")
# set.seed(123)
# results <- SimuRun(args, methods = c("DRIVE.joint", "DRIVE.ML"),
#                    ml_fitting_propensity = ml_propensity_tree,
#                    ml_fitting_surv = ml_survival_forest,
#                    nfolds = 5, seed = 123)
# print(results, Comp_parameters = c(0.1, 0.25, 0.25))
