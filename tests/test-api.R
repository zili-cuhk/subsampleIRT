# End-to-end checks of the user interface on a small multi-group,
# matrix-sampled 2PL/3PL data set.
library(subsampleIRT)
set.seed(1)

J <- 12L
country <- sample(c("A", "B", "C"), 4000, replace = TRUE)
sim <- irt_simulate(4000, difficulty = rep(seq(-1.5, 1.5, length.out = 6), 2),
                    discrimination = rep(seq(0.8, 1.8, length.out = 6), 2),
                    guessing = rep(c(0.2, 0.2, 0, 0, 0, 0), 2), dimension = rep(1:2, each = 6),
                    correlation = matrix(c(1, 0.5, 0.5, 1), 2), group = country,
                    means = rbind(A = c(0, 0), B = c(0.4, 0.2), C = c(-0.3, 0.1)),
                    booklets = list(c(1:4, 7:10), c(3:6, 9:12), c(1, 2, 5, 6, 7, 8, 11, 12)))
Y <- sim$responses

dat <- irt_data(Y, group = country, booklet = sim$booklet, weights = runif(4000, 0.5, 1.5))
print(dat)
mod <- irt_model(dat, dimension = rep(1:2, each = 6), guessing = rep(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), 2),
                 nodes = 7)
print(mod)
stopifnot(mod$spec$p == 2 * J + 4 + 1 + 4)

full <- irt_fit(dat, mod)
print(full)
stopifnot(full$fit$converged, full$fit$projected_gradient_inf < 1e-4)
cf <- coef(full)
stopifnot(nrow(cf$items) == J, abs(sum(cf$means[, 1] * tapply(dat$weights, dat$group, sum))) < 1e-8)

des <- sub_design(dat, mod, targets = c("difficulty", "joint", "mean"), pilot_size = 1500, seed = 2)
print(des)
pr <- des$design$designs$difficulty$probability
stopifnot(abs(sum(pr) - 1) < 1e-12, min(pr) >= 0.1 / dat$n - 1e-15)

sf <- sub_fit(des, size = 1500, target = "difficulty", seed = 3)
s <- summary(sf)
print(head(s))
stopifnot(nrow(s) == J, all(s$se > 0), all(is.finite(s$estimate)))
V <- vcov(sf, "mean")
stopifnot(nrow(V) == 6, isSymmetric(V))
# the subsample estimate is close to the full-data estimate
err <- target_estimate(sf, "difficulty") - target_estimate(full, "difficulty")
stopifnot(max(abs(err) / s$se) < 6)

su <- sub_fit(des, size = 1500, target = "difficulty", method = "stratified", se = FALSE, seed = 3)
stopifnot(is.null(su$sandwich))

ts <- sub_size(des, target = "difficulty", tolerance = 0.08, stage1 = 800, criterion = "mean", seed = 4)
print(ts)
stopifnot(ts$size >= 800)

path <- sub_path(des, size = 1500, target = "joint", seed = 5,
                 penalties = list(irt_penalty(guessing = 0.002), irt_penalty(guessing = 0.02)))
stopifnot(length(path) == 2L)
g <- sapply(path, function(f) target_estimate(f, "guessing"))
print(g)
cat("API checks passed\n")
