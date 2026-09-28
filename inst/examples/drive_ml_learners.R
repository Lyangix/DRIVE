# Example nuisance learners for DRIVE.ML. Sourcing defines functions only.
# Install the optional packages first: install.packages(c("rpart", "randomForestSRC"))

ml_propensity_tree <- function(data, predictx) {
  if (!requireNamespace("rpart", quietly = TRUE)) stop("Install rpart first.")
  # data contains training IV and covariates; predictx has covariates only.
  data$IV <- factor(data$IV, levels = c(0, 1))
  if (length(unique(data$IV)) < 2L) stop("Both treatment arms are needed in each training fold.")
  model <- rpart::rpart(IV ~ ., data = data, method = "class",
                        control = rpart::rpart.control(cp = 0.01, minbucket = 10))
  probability <- stats::predict(model, newdata = predictx, type = "prob")[, "1"]
  # Example stabilization; choose clipping/tuning for your study within training folds.
  pmin(0.99, pmax(0.01, as.numeric(probability)))
}

# Map a cumulative hazard onto DRIVE's grid, then difference it. The first
# column is H(stime[1]) - H(0). Do not interpolate increments or repeat them
# at all later times.
drive_grid_increments <- function(chf, time.interest, stime) {
  H <- cbind(0, chf)[, findInterval(stime, time.interest) + 1L, drop = FALSE]
  previous <- cbind(0, H)[, seq_along(stime), drop = FALSE]
  H - previous
}

ml_survival_forest <- function(data, predictx) {
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    stop("Install randomForestSRC first.")
  }
  # This follows the control-only learner in the original research scripts.
  # Use only TRAINING subjects whose recorded treatment history is entirely 0.
  # Adapt the nuisance-training strategy to the assumptions of your study.
  untreated <- rowSums(data$D_status != 0) == 0
  if (sum(untreated) < 10L || sum(data$event[untreated]) < 2L) {
    stop("Too few untreated training subjects/events; use more data or adapt the nuisance learner.")
  }
  train <- data.frame(time = data$time[untreated],
                       event = as.integer(data$event[untreated]),
                       data$Covariates[untreated, , drop = FALSE])
  xnames <- paste0("X", seq_len(ncol(data$Covariates)))
  names(train) <- c("time", "event", xnames)
  names(predictx) <- xnames
  Surv <- survival::Surv
  model <- randomForestSRC::rfsrc(Surv(time, event) ~ ., data = train,
                                   ntree = 100, nodesize = 10, ntime = 0)
  prediction <- stats::predict(model, newdata = predictx)
  # randomForestSRC returns cumulative hazards in chf, on time.interest:
  # https://www.randomforestsrc.org/reference/predict.rfsrc.html
  drive_grid_increments(prediction$chf, prediction$time.interest, data$stime)
}
