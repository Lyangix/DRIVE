library(DRIVE)

# Check the compiled Newton update against an independent finite difference of
# the score, with fixed nuisance predictions. There is no treatment switching
# in this fixture, so cumulative treatment exposure is simply Z_i * t.
dat <- list(time = c(1, 1, 2, 2), event = c(1, 0, 1, 0),
            IV = c(0, 1, 0, 1), IV_c = c(-0.4, 0.6, -0.4, 0.6),
            ConfoundingPart = matrix(c(0.1, 0.2, 0.05, 0.15), 4, 2),
            D_status = matrix(c(0, 1, 0, 1), 4, 2), stime = c(1, 2))

score <- function(beta, dat) {
  previous_time <- c(0, head(dat$stime, -1L))
  increments <- dat$stime - previous_time
  event <- outer(dat$time, dat$stime, "==") * dat$event
  risk <- outer(dat$time, dat$stime, ">=")
  previous_risk <- outer(dat$time, previous_time, ">=")
  weight <- exp(beta * outer(dat$IV, dat$stime))
  previous_weight <- exp(beta * outer(dat$IV, previous_time))
  baseline <- colSums(weight * risk *
    (event - beta * outer(dat$IV, increments) - dat$ConfoundingPart)) /
    colSums(weight * risk)
  residual <- weight * (event - risk * sweep(dat$ConfoundingPart, 2, baseline, "+")) -
    previous_risk * (weight - previous_weight)
  sum(dat$IV_c * rowSums(residual))
}

native_step <- function(beta, dat) {
  do.call(getFromNamespace("driv_cf_ml_est", "DRIVE"),
          c(list(init_parameters = beta), dat,
            list(max_iter = 1L, tol = 1e-12, eta = 1e-4, contraction = 0.5)))$x
}

# Putting a different treatment arm first must not change the update. The old
# IV[j] code uses IV[0] for everyone in the first interval and fails this test.
permuted <- dat
order <- c(2, 1, 3, 4)
for (name in c("time", "event", "IV", "IV_c")) {
  permuted[[name]] <- dat[[name]][order]
}
for (name in c("ConfoundingPart", "D_status")) {
  permuted[[name]] <- dat[[name]][order, , drop = FALSE]
}

for (beta in c(-0.2, 0, 0.2)) {
  h <- 1e-6
  derivative <- (score(beta + h, dat) - score(beta - h, dat)) / (2 * h)
  expected <- beta - score(beta, dat) / derivative
  actual <- native_step(beta, dat)
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-8)),
            isTRUE(all.equal(actual, native_step(beta, permuted), tolerance = 1e-12)))
}
