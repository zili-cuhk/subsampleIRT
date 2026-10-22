# subsampleIRT

Target-specific optimal subsampling for penalized marginal maximum likelihood
in multidimensional, multi-group mixed 2PL/3PL item response models.

The computations are those used for the simulation studies and the TIMSS 2023
application of the accompanying paper (*Target-Specific Optimal Subsampling
for Penalized Marginal-Likelihood Estimation with Numerical Integration*).

## Installation

```sh
R CMD build subsampleIRT
R CMD INSTALL subsampleIRT_0.1.1.tar.gz
```

Only base R and the recommended package Matrix are needed.

## Example

```r
library(subsampleIRT)
sim <- irt_simulate(20000, difficulty = seq(-2, 2, length.out = 20),
                    discrimination = rep(1.2, 20), guessing = rep(c(0.2, 0), 10),
                    seed = 1)
dat <- irt_data(sim$responses)
mod <- irt_model(dat, guessing = rep(c(TRUE, FALSE), 10))

# pilot fit and probabilities for two targets
des <- sub_design(dat, mod, targets = c("difficulty", "guessing"), pilot_size = 2000, seed = 2)
des

# subsample estimate of the difficulties with conditional sandwich SEs
fit <- sub_fit(des, size = 3000, target = "difficulty", seed = 3)
summary(fit)

# subsample size for a standard-error tolerance of 0.02
fit2 <- sub_size(des, target = "difficulty", tolerance = 0.02, stage1 = 1000)
fit2$size
```

Multi-group data with booklets and survey weights:

```r
dat <- irt_data(Y, group = country, booklet = booklet, weights = senwgt)
mod <- irt_model(dat, dimension = domain, guessing = item_type == "MC", means = "sum_zero")
```

`inst/examples/timss2023.R` runs the TIMSS 2023 analysis of the paper.

## Main functions

| Function | Purpose |
|---|---|
| `irt_data()` | responses, groups, booklets and base weights |
| `irt_model()`, `irt_penalty()` | model and penalty |
| `irt_fit()` | full-data penalized MML |
| `sub_design()` | pilot fit, score pass, target-specific probabilities |
| `sub_fit()` | subsample fit with sandwich standard errors |
| `sub_size()` | two-stage choice of the subsample size |
| `sub_path()` | penalty path on one subsample |
| `irt_simulate()` | simulated responses |

## Tests

`tests/test-core.R` checks analytic gradients against finite differences and
summed person scores, the sparse and dense evaluation engines, the group-mean
constraints, the correlation parameterization and the sampling probabilities.
Set `SUBSAMPLEIRT_LONG_TESTS=true` to add parameter recovery and sandwich
calibration checks. `tests/test-api.R` runs the full workflow.
