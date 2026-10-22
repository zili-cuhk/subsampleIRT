# User-facing interface. The computations are those of core.R.

.ctrl <- function(control) {
  def <- list(maxit = 2000L, factr = 1e5, pgtol = 1e-8, chunk = 20000L)
  def[names(control)] <- control
  def
}

.target_names <- c("difficulty", "discrimination", "guessing", "correlation",
                   "mean", "joint", "native")

#' Prepare item responses
#'
#' Collects dichotomous item responses, optional group (for example country)
#' and booklet (test form) identifiers and base weights for the other
#' functions of the package.
#'
#' @param responses A matrix or data frame of 0/1 item responses with
#'   \code{NA} for items that were not administered or not reached. Rows are
#'   persons and columns are items.
#' @param group Optional vector of group labels (for example countries), one
#'   per person. Groups have their own latent means in \code{irt_model()}.
#' @param booklet Optional vector of booklet (test form) labels, one per person.
#'   Persons with the same group and booklet are processed together, which is
#'   much faster for matrix-sampled designs.
#' @param weights Optional positive base weights (for example survey weights).
#' @return An object of class \code{"irt_data"}.
#' @examples
#' sim <- irt_simulate(500, difficulty = seq(-1, 1, length.out = 6),
#'                     discrimination = rep(1.2, 6), seed = 1)
#' dat <- irt_data(sim$responses)
#' dat
#' @export
irt_data <- function(responses, group = NULL, booklet = NULL, weights = NULL) {
  Y <- as.matrix(responses)
  storage.mode(Y) <- "double"
  ok <- is.na(Y) | Y == 0 | Y == 1
  if (!all(ok)) stop("responses must be 0, 1 or NA.", call. = FALSE)
  n <- nrow(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("item", seq_len(ncol(Y)))
  if (any(rowSums(!is.na(Y)) == 0)) stop("every person needs at least one response.", call. = FALSE)
  if (any(colSums(!is.na(Y)) == 0)) stop("every item needs at least one response.", call. = FALSE)
  if (is.null(weights)) weights <- rep(1, n)
  stopifnot(length(weights) == n, all(is.finite(weights)), all(weights > 0))
  if (is.null(group)) {
    gf <- factor(rep("all", n))
  } else {
    stopifnot(length(group) == n, !anyNA(group))
    gf <- factor(group)
  }
  bk <- NULL
  if (!is.null(booklet)) {
    stopifnot(length(booklet) == n, !anyNA(booklet))
    bk <- interaction(gf, factor(booklet), drop = TRUE)
  }
  resp <- prepare_responses(Y, group = if (nlevels(gf) > 1L) as.integer(gf) else NULL,
                            block = if (is.null(bk)) NULL else as.integer(bk))
  structure(list(resp = resp, weights = as.numeric(weights), group = gf,
                 booklet = bk, items = colnames(Y), n = n, J = ncol(Y)),
            class = "irt_data")
}

#' @export
print.irt_data <- function(x, ...) {
  cat(sprintf("Item response data: %d persons, %d items, %d observed responses\n",
              x$n, x$J, sum(x$resp$M)))
  if (nlevels(x$group) > 1L) cat(sprintf("  groups: %d\n", nlevels(x$group)))
  if (!is.null(x$booklet)) cat(sprintf("  group-by-booklet blocks: %d\n", nlevels(x$booklet)))
  invisible(x)
}

#' Penalty for discrimination and guessing parameters
#'
#' A log-normal penalty on the discrimination parameters and a Beta penalty
#' on the free guessing parameters, both on the per-person (mean) scale of
#' the criterion, so that their influence does not change with the sample
#' size or the subsample size.
#'
#' @param guessing Strength of the Beta penalty on guessing parameters
#'   (0 for none).
#' @param guessing_shape The two Beta shape parameters.
#' @param log_slope Strength of the normal penalty on log-discriminations.
#' @param log_slope_sd Its standard deviation.
#' @return A list used by \code{irt_model()}.
#' @examples
#' irt_penalty(guessing = 0.01)
#' @export
irt_penalty <- function(guessing = 0.005, guessing_shape = c(2.5, 8.5),
                        log_slope = 0.001, log_slope_sd = 0.5) {
  structure(mirt_penalty(guessing_strength = guessing, guessing_alpha = guessing_shape[1],
                         guessing_beta = guessing_shape[2],
                         log_discrimination_strength = log_slope,
                         log_discrimination_sd = log_slope_sd),
            class = c("irt_penalty", "list"))
}

#' Specify a (multidimensional, multi-group) item response model
#'
#' Between-item model in which item \eqn{j} measures dimension \eqn{k(j)}:
#' \eqn{P(Y_{ij}=1\mid\theta)=c_j+(1-c_j)\,\mathrm{logit}^{-1}(d_j+a_j\theta_{k(j)})},
#' with free guessing parameters \eqn{c_j} for the items in \code{guessing}
#' (3PL) and \eqn{c_j=0} otherwise (2PL). The latent vector is normal with
#' unit variances, a correlation matrix shared by all groups and
#' group-specific means. Integration uses a tensor-product Gauss--Hermite
#' rule.
#'
#' @param data An \code{irt_data} object.
#' @param dimension Integer vector giving the dimension of each item
#'   (default: all items on one dimension).
#' @param guessing Logical vector (or a single value) marking the items with a
#'   free guessing parameter.
#' @param means Identification of the group means: \code{"sum_zero"} (the
#'   weighted average of the group means is zero) or \code{"reference"} (the
#'   first group has mean zero).
#' @param nodes Number of quadrature nodes per dimension.
#' @param penalty An \code{irt_penalty()} object.
#' @return An object of class \code{"irt_model"}.
#' @examples
#' sim <- irt_simulate(500, difficulty = seq(-1, 1, length.out = 6),
#'                     discrimination = rep(1.2, 6), seed = 1)
#' mod <- irt_model(irt_data(sim$responses))
#' mod
#' @export
irt_model <- function(data, dimension = NULL, guessing = FALSE,
                      means = c("sum_zero", "reference"), nodes = 11L,
                      penalty = irt_penalty()) {
  stopifnot(inherits(data, "irt_data"))
  means <- match.arg(means)
  J <- data$J
  if (is.null(dimension)) dimension <- rep(1L, J)
  dimension <- as.integer(dimension)
  if (length(guessing) == 1L) guessing <- rep(guessing, J)
  stopifnot(length(dimension) == J, length(guessing) == J, all(dimension >= 1L))
  D <- max(dimension)
  if (!all(seq_len(D) %in% dimension)) stop("every dimension needs at least one item.", call. = FALSE)
  G <- nlevels(data$group)
  mass <- as.numeric(tapply(data$weights, data$group, sum))
  spec <- mirt_spec(dimension, as.logical(guessing), D = D, G = G,
                    mean_constraint = means, group_mass = mass)
  structure(list(spec = spec, quad = gh_normal_product(as.integer(nodes), D),
                 penalty = penalty, nodes = as.integer(nodes), means = means,
                 items = data$items, groups = levels(data$group)),
            class = "irt_model")
}

#' @export
print.irt_model <- function(x, ...) {
  s <- x$spec
  cat(sprintf("IRT model: %d items (%d with free guessing), %d dimension(s), %d group(s)\n",
              s$J, s$nF, s$D, s$G))
  cat(sprintf("  %d parameters, %d quadrature nodes per dimension (%d in total)\n",
              s$p, x$nodes, length(x$quad$weight)))
  invisible(x)
}

# Readable parameter tables from eta.
.param_table <- function(eta, model) {
  s <- model$spec
  par <- parse_eta(eta, s)
  items <- data.frame(item = model$items, dimension = s$dim, a = par$alpha,
                      b = par$difficulty, d = par$d, c = par$guessing,
                      row.names = NULL)
  R <- par$correlation
  dimnames(R) <- list(paste0("F", seq_len(s$D)), paste0("F", seq_len(s$D)))
  mu <- par$mu
  dimnames(mu) <- list(model$groups, paste0("F", seq_len(s$D)))
  list(items = items, correlation = R, means = mu)
}

.target_labels <- function(target, model) {
  s <- model$spec
  pairs <- if (s$nR) corr_pairs(s$D) else matrix(0L, 0, 2)
  lab <- list(
    difficulty = paste0("b:", model$items),
    discrimination = paste0("a:", model$items),
    guessing = paste0("c:", model$items[s$free_index]),
    correlation = if (s$nR) sprintf("rho:F%d,F%d", pairs[, 2], pairs[, 1]) else character(0),
    mean = if (s$G > 1L) {
      g <- rep(model$groups, each = s$D); k <- rep(seq_len(s$D), s$G)
      paste0("mu:", g, ",F", k)[mean_rows(s)]
    } else character(0))
  if (target == "joint") return(c(lab$difficulty, lab$discrimination, lab$guessing, lab$correlation))
  if (target == "native") return(paste0("eta", seq_len(s$p)))
  lab[[target]]
}

#' Full-data penalized marginal maximum likelihood
#'
#' Fits the model to all persons by minimizing the weighted penalized
#' negative marginal log likelihood with L-BFGS-B and analytic gradients.
#'
#' @param data An \code{irt_data} object.
#' @param model An \code{irt_model} object.
#' @param start Optional starting values on the optimization scale.
#' @param control List with \code{maxit}, \code{factr}, \code{pgtol} and
#'   \code{chunk} (rows processed at a time).
#' @return An object of class \code{"irt_fit"}; use \code{coef()} for the
#'   parameter tables.
#' @examples
#' sim <- irt_simulate(1000, difficulty = seq(-1, 1, length.out = 6),
#'                     discrimination = rep(1.2, 6), seed = 1)
#' dat <- irt_data(sim$responses)
#' fit <- irt_fit(dat, irt_model(dat, nodes = 7))
#' coef(fit)$items
#' @export
irt_fit <- function(data, model, start = NULL, control = list()) {
  stopifnot(inherits(data, "irt_data"), inherits(model, "irt_model"))
  ctrl <- .ctrl(control)
  f <- fit_mirt(data$resp, model$spec, model$quad, model$penalty, omega = data$weights,
                start = start, maxit = ctrl$maxit, factr = ctrl$factr, pgtol = ctrl$pgtol,
                chunk = ctrl$chunk)
  structure(list(fit = f, eta = f$eta, model = model, n = data$n, size = NA_integer_,
                 kind = "full data"), class = "irt_fit")
}

#' Parameter tables of a fit
#'
#' @param object An \code{irt_fit} or \code{sub_fit} object.
#' @param ... Not used.
#' @return A list with \code{items} (discrimination \code{a}, difficulty
#'   \code{b}, intercept \code{d} and guessing \code{c} of each item), the
#'   latent \code{correlation} matrix and the group \code{means}.
#' @export
coef.irt_fit <- function(object, ...) .param_table(object$eta, object$model)

#' @export
print.irt_fit <- function(x, ...) {
  f <- x$fit
  cat(sprintf("Penalized MML fit (%s): objective %.6f, %s after %d evaluations (%.1f s)\n",
              x$kind, f$objective, if (isTRUE(f$converged)) "converged" else "NOT converged",
              as.integer(f$evaluations), f$elapsed))
  if (f$boundary_contacts > 0) cat(sprintf("  %d parameter(s) at a bound\n", f$boundary_contacts))
  invisible(x)
}

#' Target value
#'
#' Evaluates a target function of the parameters at a fit.
#'
#' @param object An \code{irt_fit} or \code{sub_fit} object.
#' @param target One of \code{"difficulty"}, \code{"discrimination"},
#'   \code{"guessing"}, \code{"correlation"}, \code{"mean"} or \code{"joint"}.
#' @return A named numeric vector.
#' @export
target_estimate <- function(object, target = "difficulty") {
  target <- match.arg(target, .target_names)
  v <- target_value(object$eta, object$model$spec, target)
  stats::setNames(v, .target_labels(target, object$model))
}

#' Pilot fit and target-specific sampling probabilities
#'
#' Fits the model to a uniform pilot sample, evaluates the person-level
#' scores of all persons once at the pilot estimate and computes
#' trace-optimal sampling probabilities for each target, subject to a lower
#' bound and to fixed stratum masses.
#'
#' @param data An \code{irt_data} object.
#' @param model An \code{irt_model} object.
#' @param targets Targets for which probabilities are computed: any of
#'   \code{"difficulty"}, \code{"discrimination"}, \code{"guessing"},
#'   \code{"correlation"}, \code{"mean"} and \code{"joint"} (difficulty,
#'   discrimination, guessing and correlations on scales 1, 0.5, 0.1 and 0.1).
#' @param pilot_size Number of persons drawn uniformly (with replacement) for
#'   the pilot fit.
#' @param strata Optional stratum labels whose weight shares are kept fixed.
#'   By default the groups when there are several, otherwise the booklets.
#' @param lower_bound Probabilities are at least \code{lower_bound / N}.
#' @param seed Optional random seed.
#' @param control Optimizer control, see \code{irt_fit()}.
#' @return An object of class \code{"sub_design"}.
#' @examples
#' sim <- irt_simulate(2000, difficulty = seq(-1.5, 1.5, length.out = 8),
#'                     discrimination = rep(1.2, 8), seed = 1)
#' dat <- irt_data(sim$responses)
#' des <- sub_design(dat, irt_model(dat, nodes = 7), targets = "difficulty",
#'                   pilot_size = 500, seed = 2)
#' des
#' @export
sub_design <- function(data, model, targets = "joint", pilot_size = 2000L,
                       strata = NULL, lower_bound = 0.10, seed = NULL,
                       control = list()) {
  stopifnot(inherits(data, "irt_data"), inherits(model, "irt_model"))
  targets <- match.arg(targets, setdiff(.target_names, "native"), several.ok = TRUE)
  ctrl <- .ctrl(control)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(strata)) {
    strata <- if (nlevels(data$group) > 1L) data$group else data$booklet
  }
  t0 <- proc.time()[["elapsed"]]
  idx0 <- sample.int(data$n, pilot_size, replace = TRUE)
  pilot <- fit_mirt(subset_resp(data$resp, idx0), model$spec, model$quad, model$penalty,
                    omega = data$weights[idx0], maxit = ctrl$maxit, factr = ctrl$factr,
                    pgtol = ctrl$pgtol, chunk = ctrl$chunk)
  des <- design_from_pilot(pilot$eta, data$resp, model$spec, model$quad, model$penalty,
                           base_weights = data$weights, targets = targets, strata = strata,
                           lower_bound_factor = lower_bound, chunk = ctrl$chunk)
  structure(list(data = data, model = model, pilot = pilot, design = des, targets = targets,
                 strata = strata, pilot_size = as.integer(pilot_size), seconds = proc.time()[["elapsed"]] - t0, control = ctrl),
            class = "sub_design")
}

#' @export
print.sub_design <- function(x, ...) {
  cat(sprintf("Subsampling design for N = %d: pilot of %d persons (%.1f s), score pass and probabilities %.1f s\n",
              x$data$n, x$pilot_size, x$pilot$elapsed,
              x$design$score_seconds + x$design$probability_seconds))
  for (tg in x$targets) {
    d <- x$design$designs[[tg]]
    if (is.null(d)) { cat(sprintf("  %-15s (no parameters of this kind)\n", tg)); next }
    cat(sprintf("  %-15s planning variance (sum over coordinates): optimal %.4g, stratified uniform %.4g (predicted reduction %.0f%%)\n",
                tg, sum(d$planning_nu2), sum(d$planning_nu2_stratified_uniform),
                100 * (1 - sum(d$planning_nu2) / sum(d$planning_nu2_stratified_uniform))))
  }
  invisible(x)
}

.probabilities <- function(design, target, method) {
  switch(method,
         optimal = {
           if (!target %in% names(design$design$designs))
             stop("no probabilities for target '", target, "'; add it to sub_design(targets = ).", call. = FALSE)
           design$design$designs[[target]]$probability
         },
         uniform = design$design$uniform,
         stratified = design$design$stratified_uniform)
}

#' Subsample fit with conditional sandwich standard errors
#'
#' Draws a subsample with replacement from the design probabilities, fits the
#' inverse-probability weighted criterion with a fixed Horvitz--Thompson
#' normalizer (so that the penalty keeps its full-data scale), and computes
#' the centered sandwich covariance of the subsampling error.
#'
#' @param design A \code{sub_design} object.
#' @param size Subsample size \eqn{r}.
#' @param target The target whose probabilities are used.
#' @param method \code{"optimal"} (target-specific probabilities),
#'   \code{"uniform"} or \code{"stratified"} (uniform within strata).
#' @param se Compute the sandwich covariance?
#' @param hessian Finite-difference scheme for the Hessian of the sandwich.
#' @param seed Optional random seed.
#' @return An object of class \code{"sub_fit"}; see \code{summary()},
#'   \code{vcov()} and \code{target_estimate()}.
#' @examples
#' sim <- irt_simulate(2000, difficulty = seq(-1.5, 1.5, length.out = 8),
#'                     discrimination = rep(1.2, 8), seed = 1)
#' dat <- irt_data(sim$responses)
#' des <- sub_design(dat, irt_model(dat, nodes = 7), targets = "difficulty",
#'                   pilot_size = 500, seed = 2)
#' sf <- sub_fit(des, size = 800, target = "difficulty", seed = 3)
#' head(summary(sf))
#' @export
sub_fit <- function(design, size, target = design$targets[1], method = c("optimal", "uniform", "stratified"),
                    se = TRUE, hessian = c("central", "forward"), seed = NULL) {
  stopifnot(inherits(design, "sub_design"))
  method <- match.arg(method); hessian <- match.arg(hessian)
  target <- match.arg(target, setdiff(.target_names, "native"))
  m <- design$model; d <- design$data; ctrl <- design$control
  prob <- .probabilities(design, target, method)
  if (!is.null(seed)) set.seed(seed)
  idx <- draw_indices(size, prob)
  f <- fit_ipw(d$resp, idx, prob, d$weights, m$spec, m$quad, m$penalty, design$pilot$eta,
               maxit = ctrl$maxit, factr = ctrl$factr, pgtol = ctrl$pgtol, chunk = ctrl$chunk)
  sw <- if (se) sandwich_ipw(f, d$resp, m$spec, m$quad, m$penalty, chunk = ctrl$chunk,
                             hessian = hessian) else NULL
  structure(list(fit = f, eta = f$eta, sandwich = sw, model = m, size = size, target = target,
                 method = method, n = d$n, kind = sprintf("subsample r = %d, %s", size, method)),
            class = c("sub_fit", "irt_fit"))
}

#' Covariance and summary of a subsample fit
#'
#' \code{vcov()} returns the conditional sandwich covariance of the
#' subsampling error of a target; \code{summary()} returns estimates,
#' standard errors and Wald intervals. Both describe the error of the
#' subsample estimate relative to the full-data estimate.
#'
#' @param object A \code{sub_fit} object fitted with \code{se = TRUE}.
#' @param target Target (default: the target of the probabilities).
#' @param level Confidence level of the intervals.
#' @param ... Not used.
#' @return A covariance matrix (\code{vcov}) or a data frame (\code{summary}).
#' @export
vcov.sub_fit <- function(object, target = object$target, ...) {
  if (is.null(object$sandwich)) stop("refit with se = TRUE.", call. = FALSE)
  target <- match.arg(target, .target_names)
  V <- target_covariance(object$fit, object$sandwich, object$model$spec, target)
  lab <- .target_labels(target, object$model)
  dimnames(V) <- list(lab, lab)
  V
}

#' @rdname vcov.sub_fit
#' @export
summary.sub_fit <- function(object, target = object$target, level = 0.95, ...) {
  est <- target_estimate(object, target)
  out <- data.frame(parameter = names(est), estimate = unname(est), row.names = NULL)
  if (!is.null(object$sandwich)) {
    se <- sqrt(pmax(diag(vcov(object, target)), 0))
    z <- stats::qnorm(1 - (1 - level) / 2)
    out$se <- se; out$lower <- out$estimate - z * se; out$upper <- out$estimate + z * se
  }
  out
}

#' Two-stage choice of the subsample size
#'
#' Fits a stage-one subsample, estimates the per-draw variance of the target
#' with the sandwich, chooses the size \eqn{\hat r=\max\{r_1,\lceil
#' \hat\nu^2/\varepsilon^2\rceil\}} that meets the standard-error tolerance,
#' draws the remaining indices and refits.
#'
#' @param design A \code{sub_design} object.
#' @param target Target whose standard errors should meet the tolerance.
#' @param tolerance Absolute standard-error tolerance \eqn{\varepsilon}.
#' @param stage1 Stage-one size \eqn{r_1}.
#' @param criterion \code{"mean"} (average variance over the target
#'   coordinates) or \code{"max"} (largest variance).
#' @param max_size Upper limit on the subsample size.
#' @param method Probabilities, as in \code{sub_fit()}.
#' @param seed Optional random seed.
#' @return A \code{sub_fit} object with the chosen size in \code{$size}.
#' @export
sub_size <- function(design, target = design$targets[1], tolerance, stage1 = 1000L,
                     criterion = c("mean", "max"), max_size = Inf,
                     method = c("optimal", "uniform", "stratified"), seed = NULL) {
  stopifnot(inherits(design, "sub_design"))
  criterion <- match.arg(criterion); method <- match.arg(method)
  target <- match.arg(target, setdiff(.target_names, "native"))
  m <- design$model; d <- design$data; ctrl <- design$control
  prob <- .probabilities(design, target, method)
  if (!is.null(seed)) set.seed(seed)
  ts <- two_stage_fit(d$resp, prob, d$weights, m$spec, m$quad, m$penalty, target,
                      epsilon = tolerance, r1 = stage1, start = design$pilot$eta,
                      r_max = max_size, criterion = criterion, final_bread = "stage1",
                      maxit = ctrl$maxit, factr = ctrl$factr, pgtol = ctrl$pgtol, chunk = ctrl$chunk)
  structure(list(fit = ts$final, eta = ts$final$eta, sandwich = ts$sandwich, model = m,
                 size = ts$r_hat, target = target, method = method, n = d$n,
                 stage1_variance = ts$nu2_stage1, seconds = ts$seconds,
                 kind = sprintf("two-stage, r = %d", ts$r_hat)),
            class = c("sub_fit", "irt_fit"))
}

#' Amortized fitting along a penalty path
#'
#' Fits a sequence of penalties to one subsample drawn from the design,
#' warm-starting each fit from the previous one.
#'
#' @param design A \code{sub_design} object.
#' @param size Subsample size.
#' @param penalties A list of \code{irt_penalty()} objects, fitted in order.
#' @param target,method Probabilities, as in \code{sub_fit()}.
#' @param seed Optional random seed.
#' @return A list of \code{sub_fit} objects (without standard errors).
#' @export
sub_path <- function(design, size, penalties, target = design$targets[1],
                     method = c("optimal", "uniform", "stratified"), seed = NULL) {
  stopifnot(inherits(design, "sub_design"), is.list(penalties), length(penalties) >= 1L)
  method <- match.arg(method)
  target <- match.arg(target, setdiff(.target_names, "native"))
  m <- design$model; d <- design$data; ctrl <- design$control
  prob <- .probabilities(design, target, method)
  if (!is.null(seed)) set.seed(seed)
  idx <- draw_indices(size, prob)
  steps <- lapply(penalties, function(p) list(penalty = p, spec = m$spec))
  fits <- fit_path(d$resp, idx, prob, d$weights, m$quad, steps, design$pilot$eta, m$spec,
                   maxit = ctrl$maxit, factr = ctrl$factr, pgtol = ctrl$pgtol, chunk = ctrl$chunk)
  lapply(seq_along(fits), function(l) {
    mm <- m; mm$penalty <- penalties[[l]]
    structure(list(fit = fits[[l]], eta = fits[[l]]$eta, sandwich = NULL, model = mm, size = size,
                   target = target, method = method, n = d$n,
                   kind = sprintf("path step %d, r = %d", l, size)),
              class = c("sub_fit", "irt_fit"))
  })
}

#' Simulate item responses
#'
#' Simulates responses from the between-item model with normal abilities,
#' optional group means and optional booklets (each booklet administers a
#' subset of the items).
#'
#' @param n Number of persons.
#' @param difficulty,discrimination Item difficulties \eqn{b_j} and
#'   discriminations \eqn{a_j}; the intercept is \eqn{d_j=-a_jb_j}.
#' @param guessing Guessing parameters (0 for 2PL items).
#' @param dimension Dimension of each item.
#' @param correlation Latent correlation matrix.
#' @param group Optional group label of each person.
#' @param means Optional matrix of latent means with one row per group (row
#'   names matching the labels in \code{group}) and one column per dimension.
#' @param booklets Optional list of integer vectors of administered items;
#'   booklets are assigned uniformly.
#' @param seed Optional random seed.
#' @return A list with \code{responses}, \code{booklet}, \code{group} and
#'   \code{theta}.
#' @examples
#' sim <- irt_simulate(100, difficulty = c(-1, 0, 1), discrimination = c(1, 1.5, 2),
#'                     guessing = c(0.2, 0, 0))
#' colMeans(sim$responses)
#' @export
irt_simulate <- function(n, difficulty, discrimination, guessing = 0, dimension = 1L,
                         correlation = NULL, group = NULL, means = NULL,
                         booklets = NULL, seed = NULL) {
  J <- length(difficulty)
  if (length(guessing) == 1L) guessing <- rep(guessing, J)
  if (length(dimension) == 1L) dimension <- rep(as.integer(dimension), J)
  stopifnot(length(discrimination) == J, length(guessing) == J, length(dimension) == J)
  D <- max(dimension)
  if (is.null(correlation)) correlation <- diag(D)
  correlation <- as.matrix(correlation)
  stopifnot(nrow(correlation) == D)
  if (!is.null(seed)) set.seed(seed)
  theta <- matrix(stats::rnorm(n * D), n, D) %*% chol(correlation)
  if (!is.null(means)) {
    stopifnot(!is.null(group), length(group) == n)
    means <- matrix(means, ncol = D, dimnames = list(rownames(means), NULL))
    gi <- if (is.null(rownames(means))) as.integer(factor(group)) else match(as.character(group), rownames(means))
    if (anyNA(gi)) stop("every group needs a row of 'means'.", call. = FALSE)
    theta <- theta + means[gi, , drop = FALSE]
  }
  x <- sweep(sweep(theta[, dimension, drop = FALSE], 2L, discrimination, "*"), 2L,
             discrimination * difficulty, "-")
  P <- sweep(sweep(stats::plogis(x), 2L, 1 - guessing, "*"), 2L, guessing, "+")
  Y <- matrix(stats::rbinom(length(P), 1L, P), n, J)
  form <- rep(1L, n)
  if (!is.null(booklets)) {
    form <- sample.int(length(booklets), n, replace = TRUE)
    for (b in seq_along(booklets)) Y[form == b, setdiff(seq_len(J), booklets[[b]])] <- NA
  }
  colnames(Y) <- paste0("item", seq_len(J))
  list(responses = Y, booklet = form, group = group, theta = theta)
}
