#' Simulate data from multiple Gaussian graphical models
#'
#' Generates K precision matrices with controlled shared structure and draws
#' multivariate normal data from each. Follows the simulation design in
#' Peterson et al. (2015, JASA) Section 5.1. Group 1 uses the base graph
#' directly; groups 2 through K perturb the base graph by randomly flipping
#' edges with probability \code{perturb_prob}.
#'
#' @param K Integer; number of sample groups. Default 2.
#' @param p Integer; number of variables (nodes). Default 20.
#' @param n Integer (scalar or length-K vector); sample size per group.
#'   If scalar, all groups have the same sample size. Default 100.
#' @param graph_type Character; type of base graph:
#' * `band`: AR(2)-like banded structure: edges between all
#'       nodes within distance 2 (i.e., (i, i+1) and (i, i+2)). Produces
#'       approximately \code{2p - 3} edges. Matches Peterson et al. Section
#'       5.1.
#' * `random`: Erdos-Renyi random graph where each edge is
#'       included independently with probability \code{edge_prob}.
#' * `hub`: Star/hub graph with \code{floor(p/5)} hub nodes,
#'       each connected to approximately 40\% of other nodes.
#' @param edge_prob Numeric; probability of each edge in the base graph
#'   (used only for \code{graph_type = "random"}). Default 0.1.
#' @param perturb_prob Numeric; for groups \code{k = 2, ..., K}, the
#'   probability that each edge is flipped (added if absent, removed if
#'   present) relative to the base graph. Controls how different the groups
#'   are. Set to 0 for identical graphs. Default 0.1.
#' @param signal Numeric vector of length 2; magnitude range
#'   \code{c(lo, hi)} for off-diagonal precision entries where edges exist.
#'   Signs are chosen randomly. Default \code{c(0.3, 0.6)}.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return A list with components:
#' * `Omega_list`: List of K true precision matrices (each
#'       p x p, symmetric positive definite). Off-diagonal entries are
#'       non-zero only where edges exist, after row-normalization to ensure
#'       positive definiteness.
#' * `adj_list`: List of K true binary adjacency matrices (each
#'       p x p, 0/1). \code{adj_list[[k]][i,j] = 1} if edge (i,j) exists
#'       in group k.
#' * `data_list`: List of K data matrices (each n_k x p),
#'       drawn from \eqn{N(0, \Omega_k^{-1})} and column-centered.
#' * `S_list`: List of K cross-product matrices (each p x p),
#'       where \code{S_list[[k]] = t(X_k) \%*\% X_k} after centering.
#' * `n_vec`: Integer vector of length K with sample sizes.
#' * `K`: Integer; number of groups.
#' * `p`: Integer; number of variables.
#'
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 10, n = 100, seed = 42)
#' str(sim, max.level = 1)
#'
#' # Check true edge counts
#' sapply(sim$adj_list, function(a) sum(a[upper.tri(a)]))
#'
#' @seealso [multiggm_mcmc()]
#' @export
simulate_multiggm <- function(K = 2, p = 20, n = 100,
                               graph_type = c("band", "random", "hub"),
                               edge_prob = 0.1,
                               perturb_prob = 0.1,
                               signal = c(0.3, 0.6),
                               seed = NULL) {

  graph_type <- match.arg(graph_type)
  if (!is.null(seed)) set.seed(seed)
  n_vec <- rep_len(as.integer(n), K)

  # --- Step 1: Generate base adjacency matrix ---
  adj_base <- matrix(0L, p, p)

 if (graph_type == "band") {
    # AR(2) structure as in Peterson et al. Section 5.1
    for (i in 1:(p - 1)) adj_base[i, i + 1] <- adj_base[i + 1, i] <- 1L
    if (p > 2) {
      for (i in 1:(p - 2)) adj_base[i, i + 2] <- adj_base[i + 2, i] <- 1L
    }
  } else if (graph_type == "random") {
    # Erdos-Renyi
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        if (stats::runif(1) < edge_prob) {
          adj_base[i, j] <- adj_base[j, i] <- 1L
        }
      }
    }
  } else if (graph_type == "hub") {
    n_hubs <- max(1L, floor(p / 5))
    hub_nodes <- seq(1, p, length.out = n_hubs + 1)[1:n_hubs]
    hub_nodes <- as.integer(round(hub_nodes))
    for (h in hub_nodes) {
      spokes <- setdiff(seq_len(p), h)
      # Connect hub to ~40% of other nodes
      connected <- spokes[stats::runif(length(spokes)) < 0.4]
      for (s in connected) {
        adj_base[h, s] <- adj_base[s, h] <- 1L
      }
    }
  }

  # --- Step 2: Generate group-specific adjacencies by perturbing base ---
  adj_list <- vector("list", K)
  adj_list[[1]] <- adj_base

  for (k in 2:K) {
    adj_k <- adj_base
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        if (stats::runif(1) < perturb_prob) {
          # Flip edge
          adj_k[i, j] <- 1L - adj_k[i, j]
          adj_k[j, i] <- adj_k[i, j]
        }
      }
    }
    adj_list[[k]] <- adj_k
  }

  # --- Step 3: Generate precision matrices ---
  Omega_list <- vector("list", K)

  for (k in seq_len(K)) {
    Omega_k <- matrix(0, p, p)
    diag(Omega_k) <- 1

    # Fill off-diagonal entries where edges exist
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        if (adj_list[[k]][i, j] == 1L) {
          val <- stats::runif(1, signal[1], signal[2]) * sample(c(-1, 1), 1)
          Omega_k[i, j] <- val
          Omega_k[j, i] <- val
        }
      }
    }

    # Make positive definite using the approach from Peterson et al. / Danaher et al.:
    # divide each off-diagonal by (1.5 * sum of abs off-diagonals in its row), then
    # average with transpose
    Omega_k <- .make_pd(Omega_k)
    Omega_list[[k]] <- Omega_k
  }

  # --- Step 4: Draw data from N(0, Omega_k^{-1}) ---
  data_list <- vector("list", K)
  S_list    <- vector("list", K)

  for (k in seq_len(K)) {
    Sigma_k <- solve(Omega_list[[k]])
    # Ensure symmetry for chol
    Sigma_k <- (Sigma_k + t(Sigma_k)) / 2

    L <- chol(Sigma_k)  # upper Cholesky: Sigma = L'L
    Xk <- matrix(stats::rnorm(n_vec[k] * p), nrow = n_vec[k], ncol = p) %*% L
    # Column-center
    Xk <- scale(Xk, center = TRUE, scale = FALSE)
    data_list[[k]] <- Xk
    S_list[[k]]    <- crossprod(Xk)
  }

  list(
    Omega_list = Omega_list,
    adj_list   = adj_list,
    data_list  = data_list,
    S_list     = S_list,
    n_vec      = n_vec,
    K          = K,
    p          = p
  )
}


#' Simulate multi-platform data from multiple Gaussian graphical models
#'
#' Generates S platforms of data, each with K groups sharing structure across
#' groups (controlled by \code{perturb_prob}) and across platforms (controlled
#' by \code{platform_perturb_prob}). Follows Shaddox et al. (2020) simulation
#' design. Platform 1 uses the base graph; platforms 2..S perturb it.
#'
#' @param K Integer; number of sample groups. Default 2.
#' @param S Integer; number of platforms (data types). Default 2.
#' @param p_vec Integer vector of length S; number of variables per platform.
#'   Default \code{c(20, 15)}.
#' @param n Integer (scalar or length-K vector); sample size per group
#'   (same across platforms). Default 100.
#' @param graph_type Character; type of base graph. See
#'   \code{\link{simulate_multiggm}}.
#' @param perturb_prob Numeric; probability of edge flip between groups within
#'   a platform. Default 0.1.
#' @param platform_perturb_prob Numeric; probability of edge flip between
#'   platforms. Default 0.2.
#' @param signal Numeric vector of length 2; magnitude range for off-diagonal
#'   precision entries. Default \code{c(0.3, 0.6)}.
#' @param seed Optional integer random seed.
#'
#' @return A list with components:
#' * \code{platform_data}: List of S elements, each suitable as input
#'       to \code{multiggm_mcmc(method = "ssvs_platform", platform_data = ...)}.
#'       Each element is a list with \code{data_list}, \code{S_list},
#'       \code{n_vec}, \code{Omega_list}, and \code{adj_list}.
#' * \code{K}, \code{S}, \code{p_vec}: Dimensions.
#'
#' @examples
#' sim <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(10, 8),
#'                                    n = 80, seed = 42)
#' str(sim, max.level = 2)
#'
#' @seealso [simulate_multiggm()], [multiggm_mcmc()]
#' @export
simulate_multiggm_platform <- function(K = 2, S = 2, p_vec = c(20, 15),
                                        n = 100,
                                        graph_type = c("band", "random", "hub"),
                                        perturb_prob = 0.1,
                                        platform_perturb_prob = 0.2,
                                        signal = c(0.3, 0.6),
                                        seed = NULL) {

  graph_type <- match.arg(graph_type)
  if (!is.null(seed)) set.seed(seed)
  p_vec <- as.integer(p_vec)
  if (length(p_vec) != S) stop("p_vec must have length S.")
  n_vec <- rep_len(as.integer(n), K)

  platform_data <- vector("list", S)

  for (s in seq_len(S)) {
    p_s <- p_vec[s]

    # Use simulate_multiggm for this platform, optionally perturbing from
    # a common base graph pattern
    sim_s <- simulate_multiggm(
      K = K, p = p_s, n = n_vec,
      graph_type = graph_type,
      perturb_prob = perturb_prob,
      signal = signal,
      seed = NULL  # seed already set globally
    )

    # For platforms 2..S, additionally perturb the adjacency
    if (s > 1 && platform_perturb_prob > 0) {
      for (k in seq_len(K)) {
        adj_k <- sim_s$adj_list[[k]]
        for (i in 1:(p_s - 1)) {
          for (j in (i + 1):p_s) {
            if (stats::runif(1) < platform_perturb_prob) {
              adj_k[i, j] <- 1L - adj_k[i, j]
              adj_k[j, i] <- adj_k[i, j]
            }
          }
        }
        sim_s$adj_list[[k]] <- adj_k

        # Regenerate precision from new adjacency
        Omega_k <- matrix(0, p_s, p_s)
        diag(Omega_k) <- 1
        for (i in 1:(p_s - 1)) {
          for (j in (i + 1):p_s) {
            if (adj_k[i, j] == 1L) {
              val <- stats::runif(1, signal[1], signal[2]) * sample(c(-1, 1), 1)
              Omega_k[i, j] <- val
              Omega_k[j, i] <- val
            }
          }
        }
        Omega_k <- .make_pd(Omega_k)
        sim_s$Omega_list[[k]] <- Omega_k

        # Regenerate data
        Sigma_k <- solve(Omega_k)
        Sigma_k <- (Sigma_k + t(Sigma_k)) / 2
        L <- chol(Sigma_k)
        Xk <- matrix(stats::rnorm(n_vec[k] * p_s), nrow = n_vec[k], ncol = p_s) %*% L
        Xk <- scale(Xk, center = TRUE, scale = FALSE)
        sim_s$data_list[[k]] <- Xk
        sim_s$S_list[[k]] <- crossprod(Xk)
      }
    }

    platform_data[[s]] <- list(
      data_list  = sim_s$data_list,
      S_list     = sim_s$S_list,
      n_vec      = n_vec,
      Omega_list = sim_s$Omega_list,
      adj_list   = sim_s$adj_list
    )
  }

  list(
    platform_data = platform_data,
    K             = K,
    S             = S,
    p_vec         = p_vec
  )
}


#' Make a symmetric matrix positive definite
#'
#' Rescales off-diagonal entries row-wise and adds diagonal loading if needed.
#'
#' @param A Symmetric matrix with unit diagonal.
#' @param denom_factor Scale factor (larger = more shrinkage). Default 1.5.
#' @return A symmetric positive definite matrix.
#' @keywords internal
.make_pd <- function(A, denom_factor = 1.5) {
  p <- nrow(A)

  # Row-normalize off-diagonals (similar to fix_matrix but keeps original diagonal)
  for (i in seq_len(p)) {
    off_diag_sum <- sum(abs(A[i, ])) - abs(A[i, i])
    if (off_diag_sum > 0) {
      scale_factor <- 1 / (denom_factor * off_diag_sum)
      A[i, ] <- A[i, ] * scale_factor
      A[i, i] <- 1  # restore diagonal
    }
  }

  # Symmetrize
  A <- (A + t(A)) / 2

  # If still not PD, add diagonal loading
  eig_min <- min(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
  if (eig_min <= 0) {
    A <- A + (abs(eig_min) + 0.1) * diag(p)
  }

  A
}
