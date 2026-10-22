# Internal computational core of subsampleIRT. This file is the code used for
# the simulation studies and the TIMSS application of the accompanying paper;
# the user-facing functions are in api.R.
# Multidimensional penalized mixed 2PL/3PL marginal maximum likelihood with
# fixed tensor-product Gauss-Hermite quadrature, target-specific optimal
# subsampling, conditional sandwich inference, a two-stage precision rule for
# the subsample size, and amortized fitting along a penalty/model path.
#
# Base R plus the recommended Matrix package (sparse engine). The
# unidimensional case (D = 1) reproduces the penalized
# unified core of the original study (same penalty, same quadrature rule).
#
# Model (between-item structure): item j loads on dimension k(j),
#   p_j(theta) = c_j + (1 - c_j) plogis(d_j + alpha_j theta_k(j)),
#   alpha_j = exp(zeta_j), c_j = plogis(gamma_j) for free-guessing items,
#   theta ~ N_D(0, R), R = L L', L from canonical partial correlations
#   z = tanh(phi) (lower triangle, row order (2,1), (3,1), (3,2), ...).
# Parameter vector: eta = (d[J], zeta[J], gamma[free], phi[D(D-1)/2]).
#
# Criterion (minimized):
#   F(eta) = - sum_i omega_i loglik_i(eta) / normalizer + P(eta),
# with omega_i = w_i (full data / pilot) or w_i / pi_i (final IPW draw) and a
# fixed normalizer (W_N for full data, r * W_N for the final draw).

# ---------------------------------------------------------------------------
# Quadrature and latent correlation
# ---------------------------------------------------------------------------

gh_normal_1d <- function(nodes) {
  nodes <- as.integer(nodes)
  if (nodes < 3L) stop("At least three nodes are required.", call. = FALSE)
  jacobi <- matrix(0, nodes, nodes)
  off <- sqrt(seq_len(nodes - 1L))
  jacobi[cbind(seq_len(nodes - 1L), 2:nodes)] <- off
  jacobi[cbind(2:nodes, seq_len(nodes - 1L))] <- off
  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  list(node = eig$values[ord], weight = eig$vectors[1L, ord]^2)
}

# Product rule for N_D(0, I): Z is (nodes^D x D), weights sum to one.
gh_normal_product <- function(nodes, D) {
  one <- gh_normal_1d(nodes)
  grid <- as.matrix(expand.grid(rep(list(seq_len(nodes)), D)))
  Z <- matrix(one$node[as.vector(grid)], ncol = D)
  w <- apply(matrix(one$weight[as.vector(grid)], ncol = D), 1L, prod)
  list(Z = Z, weight = w / sum(w), log_weight = log(w / sum(w)),
       nodes_per_dim = nodes, D = D)
}

n_corr <- function(D) D * (D - 1L) / 2L

corr_pairs <- function(D) {
  if (D < 2L) return(matrix(integer(0), 0L, 2L))
  out <- NULL
  for (i in 2:D) for (j in seq_len(i - 1L)) out <- rbind(out, c(i, j))
  out
}

# Cholesky factor from canonical partial correlations. Written with
# arithmetic only so that it accepts complex input (complex-step derivatives).
cpc_cholesky <- function(phi, D) {
  L <- matrix(0 * phi[1L] + 0, D, D)
  if (D == 1L) { L[1L, 1L] <- 1; return(L) }
  L[1L, 1L] <- 1
  pairs <- corr_pairs(D)
  z <- tanh(phi)
  for (i in 2:D) {
    acc <- 0
    for (j in seq_len(i - 1L)) {
      m <- which(pairs[, 1L] == i & pairs[, 2L] == j)
      L[i, j] <- z[m] * sqrt(1 - acc)
      acc <- acc + L[i, j]^2
    }
    L[i, i] <- sqrt(1 - acc)
  }
  L
}

corr_from_phi <- function(phi, D) {
  L <- cpc_cholesky(phi, D)
  L %*% t(L)
}

# Complex-step derivative of a matrix-valued function of phi.
complex_step_jacobian <- function(fun, phi, h = 1e-20) {
  base <- Re(fun(phi))
  out <- vector("list", length(phi))
  for (m in seq_along(phi)) {
    step <- complex(real = phi, imaginary = 0)
    step[m] <- complex(real = phi[m], imaginary = h)
    out[[m]] <- Im(fun(step)) / h
  }
  list(value = base, derivative = out)
}

# ---------------------------------------------------------------------------
# Model specification and penalty
# ---------------------------------------------------------------------------

mirt_penalty <- function(guessing_strength = 0.005, guessing_alpha = 2.5,
                         guessing_beta = 8.5,
                         log_discrimination_strength = 0.001,
                         log_discrimination_sd = 0.5) {
  list(guessing_strength = guessing_strength, guessing_alpha = guessing_alpha,
       guessing_beta = guessing_beta,
       log_discrimination_strength = log_discrimination_strength,
       log_discrimination_sd = log_discrimination_sd)
}

# Latent means of G groups (G > 1): eta holds (G - 1) x D free values and the
# full G x D matrix is mu = mean_map %*% eta_mu (row-major by group).
#   "reference": group 1 has mean 0 (mean_map = [0; I]);
#   "sum_zero":  sum_g m_g mu_g = 0 with group masses m_g (default equal), so
#                the latent scale is centered at the mass-weighted average of
#                the groups and no single group anchors the item locations.
mirt_spec <- function(dimension_of_item, free_guessing,
                      fixed_guessing = rep(0, length(dimension_of_item)),
                      D = max(dimension_of_item), G = 1L,
                      mean_constraint = c("reference", "sum_zero"),
                      group_mass = NULL) {
  mean_constraint <- match.arg(mean_constraint)
  J <- length(dimension_of_item)
  stopifnot(length(free_guessing) == J, length(fixed_guessing) == J,
            all(dimension_of_item %in% seq_len(D)))
  free_index <- which(free_guessing)
  nF <- length(free_index)
  nR <- n_corr(D)
  nMu <- (G - 1L) * D                       # free group-mean parameters
  mean_map <- NULL
  if (G > 1L) {
    if (is.null(group_mass)) group_mass <- rep(1, G)
    group_mass <- as.numeric(group_mass) / sum(group_mass)
    mean_map <- rbind(matrix(0, D, nMu), diag(nMu))
    if (mean_constraint == "sum_zero") {
      for (g in 2:G) for (k in seq_len(D)) {
        mean_map[k, (g - 2L) * D + k] <- -group_mass[g] / group_mass[1]
      }
    }
  }
  list(J = J, D = D, G = as.integer(G), dim = as.integer(dimension_of_item),
       mean_constraint = mean_constraint, group_mass = group_mass, mean_map = mean_map,
       free_guessing = as.logical(free_guessing), free_index = free_index,
       nF = nF, nR = nR, nMu = nMu, fixed_guessing = fixed_guessing,
       p = 2L * J + nF + nR + nMu,
       idx_d = seq_len(J), idx_zeta = J + seq_len(J),
       idx_gamma = 2L * J + seq_len(nF),
       idx_phi = 2L * J + nF + seq_len(nR),
       idx_mu = 2L * J + nF + nR + seq_len(nMu))
}

parse_eta <- function(eta, spec) {
  d <- eta[spec$idx_d]
  zeta <- eta[spec$idx_zeta]
  alpha <- exp(zeta)
  guess <- spec$fixed_guessing
  gamma <- eta[spec$idx_gamma]
  if (spec$nF) guess[spec$free_index] <- plogis(gamma)
  phi <- eta[spec$idx_phi]
  G <- if (is.null(spec$G)) 1L else spec$G
  mu <- matrix(0, G, spec$D)
  if (G > 1L) {
    mm <- if (is.null(spec$mean_map)) rbind(matrix(0, spec$D, spec$nMu), diag(spec$nMu)) else spec$mean_map
    mu <- matrix(mm %*% eta[spec$idx_mu], G, spec$D, byrow = TRUE)
  }
  list(d = d, zeta = zeta, alpha = alpha, guessing = guess, gamma = gamma,
       phi = phi, mu = mu, difficulty = -d / alpha,
       correlation = if (spec$D > 1L) corr_from_phi(phi, spec$D) else matrix(1))
}

mean_map_of <- function(spec) {
  if (!is.null(spec$mean_map)) return(spec$mean_map)
  rbind(matrix(0, spec$D, spec$nMu), diag(spec$nMu))
}

penalty_value <- function(eta, spec, penalty) {
  zeta <- eta[spec$idx_zeta]
  la <- penalty$log_discrimination_strength
  sa <- penalty$log_discrimination_sd
  obj <- la / (2 * sa^2) * sum(zeta^2)
  grad <- numeric(spec$p)
  hdiag <- numeric(spec$p)
  grad[spec$idx_zeta] <- la * zeta / sa^2
  hdiag[spec$idx_zeta] <- la / sa^2
  if (spec$nF && penalty$guessing_strength > 0) {
    lc <- penalty$guessing_strength
    a <- penalty$guessing_alpha
    b <- penalty$guessing_beta
    gamma <- eta[spec$idx_gamma]
    cc <- plogis(gamma)
    obj <- obj - lc * sum((a - 1) * plogis(gamma, log.p = TRUE) +
                            (b - 1) * plogis(-gamma, log.p = TRUE))
    grad[spec$idx_gamma] <- -lc * ((a - 1) * (1 - cc) - (b - 1) * cc)
    hdiag[spec$idx_gamma] <- lc * (a + b - 2) * cc * (1 - cc)
  }
  list(objective = obj, gradient = grad, hessian_diag = hdiag)
}

# ---------------------------------------------------------------------------
# Item quantities at the quadrature nodes
# ---------------------------------------------------------------------------

# Group-invariant part (computed once per evaluation): parameters, the
# correlated node grid and its derivatives with respect to phi.
node_base <- function(eta, spec, quad) {
  par <- parse_eta(eta, spec)
  D <- spec$D
  L <- cpc_cholesky(par$phi, D)
  dTheta <- NULL
  if (D > 1L) {
    cs <- complex_step_jacobian(function(ph) cpc_cholesky(ph, D), par$phi)
    dTheta <- lapply(cs$derivative, function(dL) {
      t((quad$Z %*% t(dL))[, spec$dim, drop = FALSE])   # J x Qn
    })
  }
  list(par = par, Theta0 = quad$Z %*% t(L), dTheta = dTheta)   # Theta0: Qn x D
}

node_quantities <- function(eta, spec, quad, group = 1L, base = NULL) {
  if (is.null(base)) base <- node_base(eta, spec, quad)
  par <- base$par
  Theta <- sweep(base$Theta0, 2L, par$mu[group, ], "+")  # group-specific latent mean
  Tk <- t(Theta[, spec$dim, drop = FALSE])         # J x Qn: theta_{q,k(j)}
  x <- par$d + par$alpha * Tk                      # J x Qn
  # logistic quantities from one exp and one log1p (stable for all x):
  # s = plogis(x), log s = -(max(-x,0) + log1p(e^{-|x|})), log(1-s) likewise
  ax <- abs(x)
  e <- exp(-ax)
  l1pe <- log1p(e)
  s <- 1 / (1 + e)
  neg <- x < 0
  s[neg] <- e[neg] * s[neg]
  log_s <- -(0.5 * (ax - x) + l1pe)                # (|x| - x)/2 = max(-x, 0), exact
  log_1ms <- log_s - x
  cc <- par$guessing
  P <- cc + (1 - cc) * s
  logP <- log_s                                    # exact for 2PL rows
  g3 <- cc > 0
  if (any(g3)) logP[g3, ] <- log(P[g3, , drop = FALSE])
  log1mP <- log1p(-cc) + log_1ms
  Ds <- (1 - cc) * s * (1 - s)                     # dp/dx
  list(par = par, Tk = Tk, s = s, P = P, logP = logP, log1mP = log1mP,
       Ds = Ds, dTheta = base$dTheta)
}

# Restrict node quantities to the items observed in a chunk (sparse designs:
# with persons sorted by booklet, chunks touch only a few items).
nq_cols <- function(nq, cols) {
  list(logP = nq$logP[cols, , drop = FALSE], log1mP = nq$log1mP[cols, , drop = FALSE])
}

# Person-by-node log likelihood, log marginal likelihood and posterior.
person_posterior <- function(Yo, M, nq, log_weight) {
  posterior_from_logf(Yo %*% nq$logP + (M - Yo) %*% nq$log1mP, log_weight)
}

posterior_from_logf <- function(logf, log_weight) {   # logf: n x Qn
  n <- nrow(logf)
  logf <- logf + rep(log_weight, each = n)
  mx <- logf[cbind(seq_len(n), max.col(logf, ties.method = "first"))]
  ex <- exp(logf - mx)
  tot <- rowSums(ex)
  list(loglik = mx + log(tot), tau = ex / tot)
}

prepare_responses <- function(Y, group = NULL, block = NULL) {
  Y <- as.matrix(Y)
  M <- 1 * !is.na(Y)
  Yo <- Y
  Yo[is.na(Yo)] <- 0
  list(Yo = Yo, M = M, n = nrow(Y), J = ncol(Y),
       group = if (is.null(group)) NULL else as.integer(group),
       block = if (is.null(block)) NULL else as.integer(factor(block)),
       cache = new.env(parent = emptyenv()))
}

# Evaluation engines. "dense": rows are processed in chunks within each
# (group, block), using only the items observed in the block (efficient when
# blocks are large, e.g. the full data). "sparse": each group is processed in
# one sparse-dense product over its observed responses (efficient when blocks
# are small, e.g. subsamples, where per-chunk overhead dominates). "auto"
# chooses per group by the mean number of rows per block. Both engines give
# the same values up to rounding. Set with options(subsampleIRT.engine = ...).
engine_choice <- function(rows_per_block) {
  eng <- getOption("subsampleIRT.engine", "auto")
  if (eng == "auto") {
    if (rows_per_block < getOption("subsampleIRT.sparse_rows", 128)) "sparse" else "dense"
  } else eng
}

# Per-response-set cache (groups, block sizes, observed count, sparse
# matrices); stored in resp$cache when present, so it is built once per fit.
eval_cache <- function(resp, G) {
  key <- paste0("G", G)
  if (!is.null(resp$cache) && !is.null(resp$cache[[key]])) return(resp$cache[[key]])
  groups <- group_index(resp, G)
  info <- lapply(groups, function(rows) {
    nb <- if (is.null(resp$block) || !length(rows)) 1L else length(unique(resp$block[rows]))
    list(rows = rows, rows_per_block = length(rows) / nb)
  })
  ec <- new.env(parent = emptyenv())
  ec$groups <- info
  ec$n_obs <- sum(resp$M)
  ec$sparse <- vector("list", length(info))
  if (!is.null(resp$cache)) resp$cache[[key]] <- ec
  ec
}

group_sparse <- function(resp, ec, gi) {
  if (!is.null(ec$sparse[[gi]])) return(ec$sparse[[gi]])
  rows <- ec$groups[[gi]]$rows
  Yo <- resp$Yo[rows, , drop = FALSE]
  Fo <- resp$M[rows, , drop = FALSE] - Yo
  sp1 <- function(A) {
    ij <- which(A != 0, arr.ind = TRUE)
    Matrix::sparseMatrix(i = ij[, 1], j = ij[, 2], x = as.numeric(A[ij]),
                         dims = dim(A))
  }
  sp <- list(Y = sp1(Yo), F = sp1(Fo))
  ec$sparse[[gi]] <- sp
  sp
}

# Row chunks within a set of rows: split by resp$block (e.g. booklet) when
# present, so that each chunk touches only that block's items, then by size.
row_chunks <- function(resp, rows, chunk) {
  if (is.null(resp$block)) {
    return(lapply(chunk_index(length(rows), chunk), function(ch) rows[ch]))
  }
  out <- list()
  for (b in split(rows, resp$block[rows])) {
    if (!length(b)) next
    out <- c(out, lapply(chunk_index(length(b), chunk), function(ch) b[ch]))
  }
  out
}

# Person indices by latent group (a single group when resp$group is NULL).
group_index <- function(resp, G = 1L) {
  if (is.null(resp$group) || G == 1L) return(list(`1` = seq_len(resp$n)))
  out <- split(seq_len(resp$n), factor(resp$group, levels = seq_len(G)))
  out
}

chunk_index <- function(n, chunk) split(seq_len(n), ceiling(seq_len(n) / chunk))

# ---------------------------------------------------------------------------
# Objective, gradient and individual scores
# ---------------------------------------------------------------------------

new_counter <- function() {
  e <- new.env(parent = emptyenv())
  e$evaluations <- 0
  e$person_node_item <- 0
  e
}

evaluate_mirt <- function(eta, resp, omega, normalizer, spec, quad, penalty,
                          chunk = 20000L, counter = NULL) {
  J <- spec$J
  G <- if (is.null(spec$G)) 1L else spec$G
  g <- numeric(spec$p)
  ll <- 0
  ec <- eval_cache(resp, G)
  base <- node_base(eta, spec, quad)
  gmu <- numeric(G * spec$D)                     # gradient w.r.t. all G x D means
  for (gi in seq_len(G)) {
    rows <- ec$groups[[gi]]$rows
    if (!length(rows)) next
    nq <- node_quantities(eta, spec, quad, gi, base = base)
    if (engine_choice(ec$groups[[gi]]$rows_per_block) == "sparse") {
      sp <- group_sparse(resp, ec, gi)
      post <- posterior_from_logf(as.matrix(sp$Y %*% nq$logP) +
                                    as.matrix(sp$F %*% nq$log1mP), quad$log_weight)
      wt <- post$tau * omega[rows]
      R1 <- as.matrix(Matrix::crossprod(sp$Y, wt))
      N1 <- R1 + as.matrix(Matrix::crossprod(sp$F, wt))
      ll <- ll + sum(omega[rows] * post$loglik)
    } else {
      R1 <- matrix(0, J, length(quad$weight))
      N1 <- R1
      for (idx in row_chunks(resp, rows, chunk)) {
        cols <- which(colSums(resp$M[idx, , drop = FALSE]) > 0)
        Yo <- resp$Yo[idx, cols, drop = FALSE]
        M <- resp$M[idx, cols, drop = FALSE]
        post <- person_posterior(Yo, M, nq_cols(nq, cols), quad$log_weight)
        wt <- post$tau * omega[idx]
        R1[cols, ] <- R1[cols, ] + crossprod(Yo, wt)
        N1[cols, ] <- N1[cols, ] + crossprod(M, wt)
        ll <- ll + sum(omega[idx] * post$loglik)
      }
    }
    Gq <- R1 / nq$P - (N1 - R1) / (1 - nq$P)          # J x Qn
    GD <- Gq * nq$Ds
    par <- nq$par
    g[spec$idx_d] <- g[spec$idx_d] + rowSums(GD)
    g[spec$idx_zeta] <- g[spec$idx_zeta] + rowSums(GD * nq$Tk) * par$alpha
    if (spec$nF) {
      f <- spec$free_index
      cc <- par$guessing[f]
      g[spec$idx_gamma] <- g[spec$idx_gamma] +
        rowSums(Gq[f, , drop = FALSE] * (1 - nq$s[f, , drop = FALSE])) * cc * (1 - cc)
    }
    if (spec$nR) {
      for (m in seq_len(spec$nR)) {
        g[spec$idx_phi[m]] <- g[spec$idx_phi[m]] + sum(GD * par$alpha * nq$dTheta[[m]])
      }
    }
    if (G > 1L) {
      ga <- rowSums(GD) * par$alpha                       # per item
      for (k in seq_len(spec$D)) {
        gmu[(gi - 1L) * spec$D + k] <- gmu[(gi - 1L) * spec$D + k] + sum(ga[spec$dim == k])
      }
    }
  }
  if (G > 1L) g[spec$idx_mu] <- g[spec$idx_mu] + as.vector(crossprod(mean_map_of(spec), gmu))
  if (!is.null(counter)) {
    counter$evaluations <- counter$evaluations + 1
    counter$person_node_item <- counter$person_node_item +
      ec$n_obs * length(quad$weight)
  }
  pen <- penalty_value(eta, spec, penalty)
  list(objective = -ll / normalizer + pen$objective,
       gradient = -g / normalizer + pen$gradient,
       loglik_sum = ll, penalty = pen)
}

# Unweighted per-person scores u_i (n x p) at eta.
person_scores <- function(eta, resp, spec, quad, chunk = 20000L) {
  G <- if (is.null(spec$G)) 1L else spec$G
  U <- matrix(0, resp$n, spec$p)
  groups <- group_index(resp, G)
  base <- node_base(eta, spec, quad)
  for (gi in seq_len(G)) {
    rows <- groups[[gi]]
    if (!length(rows)) next
    nq <- node_quantities(eta, spec, quad, gi, base = base)
    par <- nq$par
    invP <- 1 / nq$P
    inv1mP <- 1 / (1 - nq$P)
    Dd1 <- t(nq$Ds * invP); Dd0 <- t(nq$Ds * inv1mP)                    # Qn x J
    Tz <- nq$Ds * nq$Tk * par$alpha
    Dz1 <- t(Tz * invP); Dz0 <- t(Tz * inv1mP)
    if (spec$nF) {
      f <- spec$free_index
      cc <- par$guessing[f]
      Tg <- (1 - nq$s[f, , drop = FALSE]) * (cc * (1 - cc))
      Dg1 <- t(Tg * invP[f, , drop = FALSE]); Dg0 <- t(Tg * inv1mP[f, , drop = FALSE])
    }
    if (spec$nR) {
      Dp <- lapply(nq$dTheta, function(dT) {
        Tp <- nq$Ds * par$alpha * dT
        list(one = t(Tp * invP), zero = t(Tp * inv1mP))
      })
    }
    for (idx in row_chunks(resp, rows, chunk)) {
      cols <- which(colSums(resp$M[idx, , drop = FALSE]) > 0)
      Yo <- resp$Yo[idx, cols, drop = FALSE]
      M <- resp$M[idx, cols, drop = FALSE]
      Mo <- M - Yo
      tau <- person_posterior(Yo, M, nq_cols(nq, cols), quad$log_weight)$tau
      Ud <- Yo * (tau %*% Dd1[, cols, drop = FALSE]) - Mo * (tau %*% Dd0[, cols, drop = FALSE])
      U[idx, spec$idx_d[cols]] <- Ud
      U[idx, spec$idx_zeta[cols]] <- Yo * (tau %*% Dz1[, cols, drop = FALSE]) -
        Mo * (tau %*% Dz0[, cols, drop = FALSE])
      if (spec$nF) {
        fc <- which(f %in% cols)                    # positions within free set
        if (length(fc)) {
          jc <- match(f[fc], cols)                  # their columns in the chunk
          U[idx, spec$idx_gamma[fc]] <- Yo[, jc, drop = FALSE] * (tau %*% Dg1[, fc, drop = FALSE]) -
            Mo[, jc, drop = FALSE] * (tau %*% Dg0[, fc, drop = FALSE])
        }
      }
      if (spec$nR) {
        for (m in seq_len(spec$nR)) {
          U[idx, spec$idx_phi[m]] <- rowSums(Yo * (tau %*% Dp[[m]]$one[, cols, drop = FALSE]) -
                                               Mo * (tau %*% Dp[[m]]$zero[, cols, drop = FALSE]))
        }
      }
      if (G > 1L) {
        Bg <- mean_map_of(spec)[(gi - 1L) * spec$D + seq_len(spec$D), , drop = FALSE]
        nz <- which(colSums(Bg != 0) > 0)
        if (length(nz)) {
          Ua <- sweep(Ud, 2L, par$alpha[cols], "*")
          dc <- spec$dim[cols]
          uk <- vapply(seq_len(spec$D), function(k) rowSums(Ua[, dc == k, drop = FALSE]),
                       numeric(length(idx)))
          uk <- matrix(uk, length(idx), spec$D)
          U[idx, spec$idx_mu[nz]] <- uk %*% Bg[, nz, drop = FALSE]
        }
      }
    }
  }
  U
}

# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

mirt_bounds <- function(spec, intercept = c(-6, 6), discrimination = c(0.2, 4),
                        guessing = c(0.005, 0.40), phi = c(-3, 3),
                        mean = c(-3, 3)) {
  nMu <- if (is.null(spec$nMu)) 0L else spec$nMu
  lower <- c(rep(intercept[1], spec$J), rep(log(discrimination[1]), spec$J),
             rep(qlogis(guessing[1]), spec$nF), rep(phi[1], spec$nR),
             rep(mean[1], nMu))
  upper <- c(rep(intercept[2], spec$J), rep(log(discrimination[2]), spec$J),
             rep(qlogis(guessing[2]), spec$nF), rep(phi[2], spec$nR),
             rep(mean[2], nMu))
  list(lower = lower, upper = upper)
}

default_start <- function(resp, omega, spec, penalty, guess0 = NULL,
                          alpha0 = 1) {
  prop <- colSums(resp$Yo * omega) / colSums(resp$M * omega)
  c0 <- spec$fixed_guessing
  if (is.null(guess0)) {
    guess0 <- (penalty$guessing_alpha - 1) /
      (penalty$guessing_alpha + penalty$guessing_beta - 2)
  }
  if (spec$nF) c0[spec$free_index] <- guess0
  adj <- pmin(pmax((prop - c0) / (1 - c0), 0.02), 0.98)
  c(qlogis(adj), rep(log(alpha0), spec$J),
    rep(qlogis(guess0), spec$nF), rep(0, spec$nR),
    rep(0, if (is.null(spec$nMu)) 0L else spec$nMu))
}

projected_gradient_inf <- function(par, grad, lower, upper) {
  max(abs(par - pmin(pmax(par - grad, lower), upper)))
}

fit_mirt <- function(resp, spec, quad, penalty, omega = NULL,
                     normalizer = NULL, start = NULL, bounds = NULL,
                     maxit = 1000L, factr = 1e7, pgtol = 0,
                     chunk = 20000L, counter = new_counter()) {
  if (is.null(omega)) omega <- rep(1, resp$n)
  if (is.null(normalizer)) normalizer <- sum(omega)
  if (is.null(start)) start <- default_start(resp, omega, spec, penalty)
  if (is.null(bounds)) bounds <- mirt_bounds(spec)
  start <- pmin(pmax(start, bounds$lower), bounds$upper)
  cache <- new.env(parent = emptyenv())
  ev <- function(par) {
    if (is.null(cache$par) || any(par != cache$par)) {
      cache$par <- par
      cache$res <- evaluate_mirt(par, resp, omega, normalizer, spec, quad,
                                 penalty, chunk, counter)
    }
    cache$res
  }
  t0 <- proc.time()[["elapsed"]]
  opt <- optim(start, function(p) ev(p)$objective, function(p) ev(p)$gradient,
               method = "L-BFGS-B", lower = bounds$lower, upper = bounds$upper,
               control = list(maxit = maxit, factr = factr, pgtol = pgtol))
  elapsed <- proc.time()[["elapsed"]] - t0
  res <- ev(opt$par)
  tol_contact <- 1e-5
  at_bound <- (opt$par - bounds$lower < tol_contact) |
    (bounds$upper - opt$par < tol_contact)
  list(eta = opt$par, objective = opt$value, gradient = res$gradient,
       projected_gradient_inf = projected_gradient_inf(opt$par, res$gradient,
                                                       bounds$lower, bounds$upper),
       gradient_inf_interior = if (any(!at_bound)) max(abs(res$gradient[!at_bound])) else 0,
       boundary_contacts = sum(at_bound), converged = opt$convergence == 0,
       convergence_code = opt$convergence, message = opt$message,
       counts = opt$counts, elapsed = elapsed, normalizer = normalizer,
       evaluations = counter$evaluations,
       person_node_item = counter$person_node_item,
       parsed = parse_eta(opt$par, spec), bounds = bounds)
}

# Multistart wrapper: keeps the converged fit with the smallest objective.
fit_mirt_multistart <- function(resp, spec, quad, penalty, starts, ...) {
  fits <- lapply(starts, function(s) fit_mirt(resp, spec, quad, penalty,
                                              start = s, ...))
  ok <- vapply(fits, function(f) f$converged, logical(1L))
  obj <- vapply(fits, function(f) f$objective, numeric(1L))
  pick <- if (any(ok)) which(ok)[which.min(obj[ok])] else which.min(obj)
  best <- fits[[pick]]
  best$start_objectives <- obj
  best$start_converged <- ok
  best$elapsed_total <- sum(vapply(fits, function(f) f$elapsed, numeric(1L)))
  best
}

# Three deterministic starts of the original 3PL study, extended with phi = 0.
deterministic_starts <- function(resp, omega, spec, penalty) {
  lapply(list(c(1 / 6, 1), c(0.05, 0.8), c(0.25, 1.2)), function(v)
    default_start(resp, omega, spec, penalty, guess0 = v[1], alpha0 = v[2]))
}

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

target_jacobian <- function(eta, spec,
                            target = c("difficulty", "discrimination",
                                       "guessing", "correlation", "joint",
                                       "native", "mean"),
                            scales = c(difficulty = 1, discrimination = 0.5,
                                       guessing = 0.1, correlation = 0.1)) {
  target <- match.arg(target)
  par <- parse_eta(eta, spec)
  J <- spec$J
  block_b <- function() {
    G <- matrix(0, J, spec$p)
    G[cbind(seq_len(J), spec$idx_d)] <- -1 / par$alpha
    G[cbind(seq_len(J), spec$idx_zeta)] <- -par$difficulty
    G
  }
  block_a <- function() {
    G <- matrix(0, J, spec$p)
    G[cbind(seq_len(J), spec$idx_zeta)] <- par$alpha
    G
  }
  block_c <- function() {
    G <- matrix(0, spec$nF, spec$p)
    if (spec$nF) {
      cc <- par$guessing[spec$free_index]
      G[cbind(seq_len(spec$nF), spec$idx_gamma)] <- cc * (1 - cc)
    }
    G
  }
  block_r <- function() {
    G <- matrix(0, spec$nR, spec$p)
    if (spec$nR) {
      pairs <- corr_pairs(spec$D)
      cs <- complex_step_jacobian(function(ph) corr_from_phi(ph, spec$D),
                                  par$phi)
      for (m in seq_len(spec$nR)) {
        G[, spec$idx_phi[m]] <- cs$derivative[[m]][pairs]
      }
    }
    G
  }
  switch(target,
         difficulty = block_b(),
         discrimination = block_a(),
         guessing = block_c(),
         correlation = block_r(),
         native = diag(spec$p),
         mean = {
           nMu <- if (is.null(spec$nMu)) 0L else spec$nMu
           if (!nMu) return(matrix(0, 0, spec$p))
           rows <- mean_rows(spec)
           G <- matrix(0, length(rows), spec$p)
           G[, spec$idx_mu] <- mean_map_of(spec)[rows, , drop = FALSE]
           G
         },
         joint = rbind(block_b() / scales[["difficulty"]],
                       block_a() / scales[["discrimination"]],
                       block_c() / scales[["guessing"]],
                       block_r() / scales[["correlation"]]))
}

# Reported group means: groups 2..G under "reference" (group 1 is 0), all G
# groups under "sum_zero"; row-major by group.
mean_rows <- function(spec) {
  all <- seq_len(spec$G * spec$D)
  if (identical(spec$mean_constraint, "sum_zero")) all else all[-seq_len(spec$D)]
}

target_value <- function(eta, spec, target,
                         scales = c(difficulty = 1, discrimination = 0.5,
                                    guessing = 0.1, correlation = 0.1)) {
  par <- parse_eta(eta, spec)
  rho <- if (spec$nR) par$correlation[corr_pairs(spec$D)] else numeric(0)
  cF <- par$guessing[spec$free_index]
  switch(target,
         difficulty = par$difficulty,
         discrimination = par$alpha,
         guessing = cF,
         correlation = rho,
         native = eta,
         mean = if (is.null(spec$nMu) || !spec$nMu) numeric(0) else
           as.vector(t(par$mu))[mean_rows(spec)],
         joint = c(par$difficulty / scales[["difficulty"]],
                   par$alpha / scales[["discrimination"]],
                   cF / scales[["guessing"]],
                   rho / scales[["correlation"]]))
}

# ---------------------------------------------------------------------------
# Sampling probabilities
# ---------------------------------------------------------------------------

lower_bounded_probabilities <- function(raw, total_mass = 1, lower_bound = 0) {
  raw <- as.numeric(raw)
  count <- length(raw)
  stopifnot(count > 0, all(is.finite(raw)), all(raw >= 0),
            count * lower_bound <= total_mass + 1e-12)
  if (total_mass - count * lower_bound <= 1e-14 || sum(raw) <= 0) {
    return(rep(total_mass / count, count))
  }
  diff_mass <- function(m) sum(pmax(lower_bound, m * raw)) - total_mass
  upper <- total_mass / sum(raw)
  while (diff_mass(upper) < 0) upper <- 2 * upper
  m <- uniroot(diff_mass, c(0, upper), tol = 1e-14)$root
  prob <- pmax(lower_bound, m * raw)
  k <- which.max(prob)
  prob[k] <- prob[k] + total_mass - sum(prob)
  prob
}

stratified_probabilities <- function(raw, strata = NULL, base_weights,
                                     lower_bound) {
  if (is.null(strata)) {
    return(lower_bounded_probabilities(raw, 1, lower_bound))
  }
  strata <- as.factor(strata)
  prob <- numeric(length(raw))
  W <- sum(base_weights)
  for (lev in levels(strata)) {
    i <- which(strata == lev)
    prob[i] <- lower_bounded_probabilities(raw[i], sum(base_weights[i]) / W,
                                           lower_bound)
  }
  prob
}

# Score pass at the pilot estimate: OPG + penalty curvature, influence
# vectors v_i = G A^{-1} u_i, probabilities and planning variances.
# Returns everything needed for several targets from one pass.
design_from_pilot <- function(pilot_eta, resp, spec, quad, penalty,
                              base_weights = rep(1, resp$n),
                              targets = c("difficulty", "guessing",
                                          "correlation", "joint"),
                              strata = NULL, lower_bound_factor = 0.10,
                              curvature = NULL, chunk = 20000L,
                              counter = NULL) {
  t0 <- proc.time()[["elapsed"]]
  U <- person_scores(pilot_eta, resp, spec, quad, chunk)
  if (!is.null(counter)) {
    counter$evaluations <- counter$evaluations + 1
    counter$person_node_item <- counter$person_node_item +
      sum(resp$M) * length(quad$weight)
  }
  t_score <- proc.time()[["elapsed"]] - t0
  W <- sum(base_weights)
  pen <- penalty_value(pilot_eta, spec, penalty)
  if (is.null(curvature)) {
    curvature <- crossprod(U * sqrt(base_weights)) / W + diag(pen$hessian_diag)
  }
  Ainv <- solve(curvature)
  N <- resp$n
  lb <- lower_bound_factor / N
  pi_unif <- if (is.null(strata)) rep(1 / N, N) else {
    s <- as.factor(strata)
    mass <- tapply(base_weights, s, sum) / W
    as.numeric(mass[as.character(s)] / table(s)[as.character(s)])
  }
  nu2_of <- function(V, prob) {
    (colSums(base_weights^2 * V^2 / prob) - colSums(base_weights * V)^2) / W^2
  }
  out <- list()
  t1 <- proc.time()[["elapsed"]]
  for (tg in targets) {
    G <- target_jacobian(pilot_eta, spec, tg)
    if (!nrow(G)) next
    V <- U %*% t(G %*% Ainv)                       # n x K, v_i = G A^{-1} u_i
    h <- base_weights * sqrt(rowSums(V^2))
    pi_opt <- stratified_probabilities(h, strata, base_weights, lb)
    out[[tg]] <- list(probability = pi_opt, h = h,
                      planning_nu2 = nu2_of(V, pi_opt),
                      planning_nu2_uniform = nu2_of(V, rep(1 / N, N)),
                      planning_nu2_stratified_uniform = nu2_of(V, pi_unif))
  }
  # Original-sample score covariance (per observation), for total intervals:
  # S0 = N / W^2 sum_i w_i^2 (u_i - grad P)(u_i - grad P)'.
  Uc <- sweep(U, 2L, pen$gradient, "-") * base_weights
  S0 <- N * crossprod(Uc) / W^2
  list(designs = out, uniform = rep(1 / N, N), stratified_uniform = pi_unif,
       curvature = curvature, S0 = S0, score_seconds = t_score,
       probability_seconds = proc.time()[["elapsed"]] - t1,
       base_weights = base_weights, W = W)
}

# ---------------------------------------------------------------------------
# Final IPW fit, sandwich, two-stage precision rule, path
# ---------------------------------------------------------------------------

draw_indices <- function(r, prob) sample.int(length(prob), r, replace = TRUE,
                                             prob = prob)

subset_resp <- function(resp, idx) {
  list(Yo = resp$Yo[idx, , drop = FALSE], M = resp$M[idx, , drop = FALSE],
       n = length(idx), J = resp$J,
       group = if (is.null(resp$group)) NULL else resp$group[idx],
       block = if (is.null(resp$block)) NULL else resp$block[idx],
       cache = new.env(parent = emptyenv()))
}

fit_ipw <- function(resp, idx, prob, base_weights, spec, quad, penalty,
                    start, ...) {
  W <- sum(base_weights)
  omega <- base_weights[idx] / prob[idx]
  fit <- fit_mirt(subset_resp(resp, idx), spec, quad, penalty, omega = omega,
                  normalizer = length(idx) * W, start = start, ...)
  fit$idx <- idx
  fit$omega <- omega
  fit$W <- W
  fit$ht_ratio <- sum(omega) / (length(idx) * W)
  fit
}

numerical_hessian <- function(gr, eta, step = 1e-5,
                              type = c("central", "forward")) {
  type <- match.arg(type)
  p <- length(eta)
  H <- matrix(0, p, p)
  g0 <- if (type == "forward") gr(eta) else NULL
  for (k in seq_len(p)) {
    h <- step * max(1, abs(eta[k]))
    e <- numeric(p); e[k] <- h
    H[, k] <- if (type == "central") (gr(eta + e) - gr(eta - e)) / (2 * h) else
      (gr(eta + e) - g0) / h
  }
  (H + t(H)) / 2
}

# `bread` may be supplied (for example the stage-one Hessian in the two-stage
# rule); otherwise the Hessian of the fixed-normalization criterion is
# computed by central differences of the analytic gradient (2p evaluations).
sandwich_ipw <- function(fit, resp, spec, quad, penalty, chunk = 20000L,
                         bread = NULL, hessian = c("central", "forward")) {
  hessian <- match.arg(hessian)
  sub <- subset_resp(resp, fit$idx)
  r <- length(fit$idx)
  gr <- function(e) evaluate_mirt(e, sub, fit$omega, fit$normalizer, spec,
                                  quad, penalty, chunk)$gradient
  A <- if (is.null(bread)) numerical_hessian(gr, fit$eta, type = hessian) else bread
  Z <- person_scores(fit$eta, sub, spec, quad, chunk) * fit$omega
  Zc <- sweep(Z, 2L, colMeans(Z), "-")
  meat <- crossprod(Zc) / (r * (r - 1) * fit$W^2)       # = B_hat / r
  Ainv <- solve(A)
  V <- Ainv %*% meat %*% t(Ainv)
  list(covariance = (V + t(V)) / 2, bread = A, meat = meat,
       min_eig_bread = min(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
}

target_covariance <- function(fit, sw, spec, target) {
  G <- target_jacobian(fit$eta, spec, target)
  G %*% sw$covariance %*% t(G)
}

# Two-stage precision rule (Theorem 1 of the extension note). The index
# stream is nested: stage-two draws continue the stage-one sequence.
two_stage_fit <- function(resp, prob, base_weights, spec, quad, penalty,
                          target, epsilon, r1, start, r_max = Inf,
                          coordinates = NULL,
                          criterion = c("max", "mean"),
                          final_bread = c("stage1", "recompute"), ...) {
  criterion <- match.arg(criterion)
  final_bread <- match.arg(final_bread)
  tm <- function() proc.time()[["elapsed"]]
  t0 <- tm()
  idx1 <- draw_indices(r1, prob)
  f1 <- fit_ipw(resp, idx1, prob, base_weights, spec, quad, penalty, start, ...)
  t1 <- tm()
  sw1 <- sandwich_ipw(f1, resp, spec, quad, penalty)
  t2 <- tm()
  v1 <- diag(target_covariance(f1, sw1, spec, target)) * r1     # nu_k^2
  if (!is.null(coordinates)) v1 <- v1[coordinates]
  need <- if (criterion == "max") max(v1) else mean(v1)
  r_hat <- min(r_max, max(r1, ceiling(need / epsilon^2)))
  idx <- c(idx1, if (r_hat > r1) draw_indices(r_hat - r1, prob) else integer(0))
  f2 <- if (r_hat > r1) {
    fit_ipw(resp, idx, prob, base_weights, spec, quad, penalty, f1$eta, ...)
  } else f1
  t3 <- tm()
  # The stage-one Hessian is a consistent bread (r1 -> infinity); reusing it
  # avoids 2p gradient evaluations on all r_hat draws.
  sw2 <- if (r_hat == r1) sw1 else sandwich_ipw(
    f2, resp, spec, quad, penalty,
    bread = if (final_bread == "stage1") sw1$bread else NULL)
  t4 <- tm()
  list(stage1 = f1, final = f2, sandwich = sw2, r_hat = r_hat,
       nu2_stage1 = v1, epsilon = epsilon, criterion = criterion,
       seconds = c(stage1_fit = t1 - t0, stage1_sandwich = t2 - t1,
                   stage2_fit = t3 - t2, final_variance = t4 - t3))
}

# Amortized fitting: one draw, many penalties/models, warm-started in order.
# `path` is a list of list(penalty = ..., spec = ...); when the spec changes
# dimension the start is rebuilt from the default start.
fit_path <- function(resp, idx, prob, base_weights, quad, path, start_eta,
                     start_spec, ...) {
  out <- vector("list", length(path))
  prev_eta <- start_eta
  prev_spec <- start_spec
  for (l in seq_along(path)) {
    sp <- path[[l]]$spec
    st <- if (identical(sp$p, prev_spec$p) &&
              identical(sp$free_index, prev_spec$free_index)) prev_eta else NULL
    if (is.null(idx)) {
      out[[l]] <- fit_mirt(resp, sp, quad, path[[l]]$penalty,
                           omega = base_weights, start = st, ...)
    } else {
      if (is.null(st)) {
        sub <- subset_resp(resp, idx)
        st <- default_start(sub, base_weights[idx] / prob[idx], sp,
                            path[[l]]$penalty)
      }
      out[[l]] <- fit_ipw(resp, idx, prob, base_weights, sp, quad,
                          path[[l]]$penalty, st, ...)
    }
    prev_eta <- out[[l]]$eta
    prev_spec <- sp
  }
  out
}

# Held-out negative log-likelihood (per person) for selection along a path.
heldout_nll <- function(eta, resp_holdout, spec, quad, chunk = 20000L) {
  G <- if (is.null(spec$G)) 1L else spec$G
  groups <- group_index(resp_holdout, G)
  tot <- 0
  for (gi in seq_len(G)) {
    rows <- groups[[gi]]
    if (!length(rows)) next
    nq <- node_quantities(eta, spec, quad, gi)
    for (idx in row_chunks(resp_holdout, rows, chunk)) {
      cols <- which(colSums(resp_holdout$M[idx, , drop = FALSE]) > 0)
      tot <- tot + sum(person_posterior(resp_holdout$Yo[idx, cols, drop = FALSE],
                                        resp_holdout$M[idx, cols, drop = FALSE],
                                        nq_cols(nq, cols), quad$log_weight)$loglik)
    }
  }
  -tot / resp_holdout$n
}

# ---------------------------------------------------------------------------
# Data generation
# ---------------------------------------------------------------------------

# items: data.frame with columns dim, difficulty, discrimination, guessing.
# booklets: list of item-index vectors; NULL for a complete design.
simulate_mirt <- function(N, items, R, booklets = NULL) {
  D <- nrow(R)
  theta <- matrix(rnorm(N * D), N, D) %*% chol(R)
  Tk <- theta[, items$dim, drop = FALSE]
  x <- sweep(sweep(Tk, 2L, items$discrimination, "*"), 2L,
             -items$discrimination * items$difficulty, "+")
  P <- sweep(plogis(x), 2L, 1 - items$guessing, "*")
  P <- sweep(P, 2L, items$guessing, "+")
  Y <- matrix(rbinom(length(P), 1L, P), N, nrow(items))
  form <- rep(1L, N)
  if (!is.null(booklets)) {
    form <- sample.int(length(booklets), N, replace = TRUE)
    for (b in seq_along(booklets)) {
      miss <- setdiff(seq_len(nrow(items)), booklets[[b]])
      Y[form == b, miss] <- NA
    }
  }
  list(Y = Y, form = form, theta = theta)
}

true_eta <- function(items, R, spec) {
  D <- nrow(R)
  phi <- numeric(spec$nR)
  if (D > 1L) {
    # invert the CPC map: z_ij = L_ij / sqrt(1 - sum_{k<j} L_ik^2)
    L <- t(chol(R))
    pairs <- corr_pairs(D)
    for (m in seq_len(nrow(pairs))) {
      i <- pairs[m, 1L]; j <- pairs[m, 2L]
      acc <- if (j > 1L) sum(L[i, seq_len(j - 1L)]^2) else 0
      phi[m] <- atanh(L[i, j] / sqrt(1 - acc))
    }
  }
  c(-items$discrimination * items$difficulty, log(items$discrimination),
    qlogis(items$guessing[spec$free_index]), phi,
    rep(0, if (is.null(spec$nMu)) 0L else spec$nMu))
}
