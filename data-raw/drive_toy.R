# Run from the repository root after installing the development package.
# Use the standard R uniform/exponential generators for a reproducible fixture.
library(DRIVE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
drive_toy <- SimulateData(N = 200, seed = 2025)
dir.create("data", showWarnings = FALSE)
save(drive_toy, file = "data/drive_toy.rda", compress = "xz", version = 2)
