# Checks of the computational core (internal functions). Base R only.
# Set SUBSAMPLEIRT_LONG_TESTS=true to add the recovery and sandwich checks
# (a few minutes).
library(subsampleIRT)
ns <- asNamespace("subsampleIRT")
for (f in ls(ns)) assign(f, get(f, envir = ns))
set.seed(20261021)

failures <- 0L
check <- function(name, value, tol) {
  ok <- is.finite(value) && value <= tol
  cat(sprintf("%-48s %11.3e  %s\n", name, value, if (ok) "ok" else "FAILED"))
  if (!ok) failures <<- failures + 1L
}

test_items <- function(D, per_dim = 6L) {
  n_free <- round(2 / 3 * per_dim)
  free <- rep(c(rep(TRUE, n_free), rep(FALSE, per_dim - n_free)), D)
  data.frame(dim = rep(seq_len(D), each = per_dim),
             difficulty = rep(seq(-2, 2, length.out = per_dim), D),
             discrimination = rep(seq(0.8, 1.8, length.out = per_dim), D),
             guessing = ifelse(free, rep(seq(0.10, 0.25, length.out = per_dim), D), 0),
             free = free)
}
test_corr <- function(D) switch(D, matrix(1), matrix(c(1, 0.5, 0.5, 1), 2),
                                matrix(c(1, 0.3, 0.5, 0.3, 1, 0.6, 0.5, 0.6, 1), 3))
test_booklets <- function(D, per_dim = 6L) {
  third <- split(seq_len(per_dim), rep(1:3, each = per_dim / 3))
  lapply(1:3, function(b) unlist(lapply(seq_len(D), function(k) (k - 1L) * per_dim + unlist(third[-b]))))
}
pen <- mirt_penalty()
fd_grad <- function(eta, fn) vapply(seq_along(eta), function(k) {
  e <- numeric(length(eta)); e[k] <- 1e-6
  (fn(eta + e) - fn(eta - e)) / 2e-6
}, numeric(1L))

# 1. Gradient vs finite differences and vs summed person scores (D = 1, 2, 3)
for (D in 1:3) {
  items <- test_items(D); R <- test_corr(D)
  sim <- simulate_mirt(1500L, items, R, if (D == 2L) test_booklets(D) else NULL)
  resp <- prepare_responses(sim$Y)
  spec <- mirt_spec(items$dim, items$free, D = D)
  quad <- gh_normal_product(5L, D)
  eta <- true_eta(items, R, spec) + rnorm(spec$p, 0, 0.05)
  omega <- runif(resp$n, 0.5, 2); nrm <- 1.7 * sum(omega)
  ev <- evaluate_mirt(eta, resp, omega, nrm, spec, quad, pen)
  fd <- fd_grad(eta, function(e) evaluate_mirt(e, resp, omega, nrm, spec, quad, pen)$objective)
  U <- person_scores(eta, resp, spec, quad, chunk = 400L)
  g2 <- -colSums(U * omega) / nrm + penalty_value(eta, spec, pen)$gradient
  check(sprintf("D=%d gradient vs finite differences", D), max(abs(ev$gradient - fd)), 1e-6)
  check(sprintf("D=%d gradient vs summed scores", D), max(abs(ev$gradient - g2)), 1e-10)
}

# 2. Multi-group means: reference and sum-zero constraints, sparse = dense engine
items <- test_items(2L); R <- test_corr(2L)
sim <- simulate_mirt(2000L, items, R, test_booklets(2L))
grp <- sample.int(3L, 2000L, replace = TRUE)
quad <- gh_normal_product(5L, 2L)
eta1 <- true_eta(items, R, mirt_spec(items$dim, items$free, D = 2L))
for (mc in c("reference", "sum_zero")) {
  resp <- prepare_responses(sim$Y, grp, block = paste(grp, sim$form))
  spec <- mirt_spec(items$dim, items$free, D = 2L, G = 3L, mean_constraint = mc,
                    group_mass = c(0.5, 0.3, 0.2))
  eta <- c(eta1, 0.3, -0.2, -0.4, 0.1) + rnorm(spec$p, 0, 0.05)
  omega <- runif(resp$n, 0.5, 2); nrm <- 1.3 * sum(omega)
  run <- function(engine, rsp, om, e = eta) {
    old <- options(subsampleIRT.engine = engine); on.exit(options(old))
    evaluate_mirt(e, rsp, om, 1.3 * sum(om), spec, quad, pen, chunk = 97L)
  }
  ev <- run("dense", resp, omega)
  fd <- fd_grad(eta, function(e) run("dense", resp, omega, e)$objective)
  U <- person_scores(eta, resp, spec, quad, chunk = 300L)
  g2 <- -colSums(U * omega) / nrm + penalty_value(eta, spec, pen)$gradient
  check(sprintf("G=3 %s: gradient vs finite differences", mc), max(abs(ev$gradient - fd)), 1e-6)
  check(sprintf("G=3 %s: gradient vs summed scores", mc), max(abs(ev$gradient - g2)), 1e-10)
  es <- run("sparse", resp, omega)
  check(sprintf("G=3 %s: sparse vs dense engine", mc),
        max(abs(c(ev$objective - es$objective, ev$gradient - es$gradient))), 1e-10)
  idx <- sample.int(resp$n, 600L, replace = TRUE)
  sub <- subset_resp(resp, idx)
  a <- run("dense", sub, omega[idx]); b <- run("sparse", sub, omega[idx]); cc <- run("auto", sub, omega[idx])
  check(sprintf("G=3 %s: engines on an IPW subsample", mc),
        max(abs(c(a$gradient - b$gradient, a$gradient - cc$gradient))), 1e-10)
  Gm <- target_jacobian(eta, spec, "mean")
  fdm <- vapply(seq_len(spec$p), function(k) {
    e <- numeric(spec$p); e[k] <- 1e-6
    (target_value(eta + e, spec, "mean") - target_value(eta - e, spec, "mean")) / 2e-6
  }, numeric(nrow(Gm)))
  check(sprintf("G=3 %s: mean-target Jacobian", mc), max(abs(Gm - fdm)), 1e-8)
  if (mc == "sum_zero") {
    mu <- parse_eta(eta, spec)$mu
    check("sum-zero: weighted average of the means", max(abs(colSums(mu * c(0.5, 0.3, 0.2)))), 1e-12)
  }
}

# 3. Correlation parameterization
phi <- c(0.3, -0.4, 0.7)
R3 <- corr_from_phi(phi, 3L)
check("CPC: unit diagonal", max(abs(diag(R3) - 1)), 1e-12)
it <- test_items(3L, per_dim = 3L)
sp3 <- mirt_spec(it$dim, rep(FALSE, 9L), D = 3L)
check("CPC: inverse map", max(abs(true_eta(transform(it, guessing = 0), R3, sp3)[sp3$idx_phi] - phi)), 1e-10)
cs <- complex_step_jacobian(function(p) corr_from_phi(p, 3L), phi)
fdR <- lapply(1:3, function(m) {
  e <- numeric(3); e[m] <- 1e-6
  (corr_from_phi(phi + e, 3L) - corr_from_phi(phi - e, 3L)) / 2e-6
})
check("CPC: complex-step derivative", max(mapply(function(a, b) max(abs(a - b)), cs$derivative, fdR)), 1e-8)

# 4. Probabilities: lower bound and stratum masses
raw <- rexp(1000L); st <- sample(letters[1:4], 1000L, TRUE); w <- runif(1000L)
pr <- stratified_probabilities(raw, st, w, 0.1 / 1000)
check("probabilities: total mass", abs(sum(pr) - 1), 1e-12)
check("probabilities: stratum masses", max(abs(tapply(pr, st, sum) - tapply(w, st, sum) / sum(w))), 1e-12)
check("probabilities: lower bound", max(0, 0.1 / 1000 - min(pr)), 1e-15)

# 5. Optional: recovery and sandwich calibration
if (identical(tolower(Sys.getenv("SUBSAMPLEIRT_LONG_TESTS")), "true")) {
  items <- test_items(2L, per_dim = 9L); R <- test_corr(2L)
  sim <- simulate_mirt(20000L, items, R, test_booklets(2L, per_dim = 9L))
  spec <- mirt_spec(items$dim, items$free, D = 2L)
  resp <- prepare_responses(sim$Y)
  quad <- gh_normal_product(11L, 2L)
  fit <- fit_mirt(resp, spec, quad, pen, factr = 1e5, pgtol = 1e-8, maxit = 2000L)
  check("recovery: converged (0 = yes)", 1 - fit$converged, 0)
  check("recovery: projected gradient", fit$projected_gradient_inf, 1e-5)
  check("recovery: |rho - 0.5|", abs(fit$parsed$correlation[2, 1] - 0.5), 0.03)
  des <- design_from_pilot(fit$eta, resp, spec, quad, pen, targets = "correlation")
  pr <- des$designs$correlation$probability
  draws <- replicate(30L, {
    f <- fit_ipw(resp, draw_indices(3000L, pr), pr, rep(1, resp$n), spec, quad, pen, fit$eta,
                 factr = 1e5, pgtol = 1e-8, maxit = 2000L)
    sw <- sandwich_ipw(f, resp, spec, quad, pen)
    c(est = f$parsed$correlation[2, 1], se = sqrt(target_covariance(f, sw, spec, "correlation")[1, 1]))
  })
  check("sandwich: |empirical SD / mean SE - 1|", abs(sd(draws["est", ]) / mean(draws["se", ]) - 1), 0.35)
}

if (failures) stop(failures, " check(s) failed")
cat("all checks passed\n")
