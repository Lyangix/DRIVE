# DRIV.s variance-estimator validation

Standalone Monte-Carlo check comparing the two variance estimators for the
`DRIV.s` treatment-effect estimator (`src/driv_s_est.cpp`):

- **original** — the simplified scalar sandwich `Asyvar = hatsigma / hatphi^2`,
  which treats the baseline hazard as fixed and ignores estimation of the
  covariate effect `alpha` and the propensity coefficient `gamma`;
- **joint** — the full multivariate sandwich over `beta = (theta, alpha, gamma)`
  described in the supplementary section *"Variance estimator under
  semi-parametric specification"*, including the baseline-profiling corrections
  (`B_Psi`, `B_Phi`) in the meat and the profiled derivatives plus the
  propensity block in the bread. This is the estimator now returned as `var`.

The harness embeds a faithful port of the estimator core (with the `IV[j]`→`IV[i]`
`j==0` indexing fix), and additionally pastes the **verbatim package block** to
assert the integrated code reproduces the validated computation bit-for-bit.

## Build & run

Requires Armadillo (`apt-get install libarmadillo-dev`).

```bash
g++ -O2 validate_var.cpp -o validate_var -larmadillo
# args: R(replicates) n(sample size) scen(0=no switch,1=switch) p(#covariates)
./validate_var 600 800 1 2
```

## Findings

Across `n = 150…800`, with and without treatment switching, and `p = 1, 2, 3`:

- the **joint** estimator's mean SE tracks the Monte-Carlo empirical SD of
  `theta_hat` (ratios ≈ 0.95–1.05) with ~94–96% CI coverage — i.e. it is valid;
- it is **essentially equivalent to the original** (agreeing to ~1% per
  replicate). This is expected: the propensity and baseline-profiling
  corrections are first-order *orthogonal*, so omitting them is asymptotically
  harmless. The joint estimator is the more complete finite-sample form.

The gold standard in each run is the empirical SD of `theta_hat` across
replicates; a valid SE estimator should match it on average.
