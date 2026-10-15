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

  adj_base <- .sim_base_graph(p, graph_type, edge_prob)
  adj_list <- .sim_group_graphs(adj_base, K, perturb_prob)
  sim <- .sim_from_graphs(adj_list, n_vec, signal)

  list(
    Omega_list = sim$Omega_list,
    adj_list   = adj_list,
    data_list  = sim$data_list,
    S_list     = sim$S_list,
    n_vec      = n_vec,
    K          = K,
    p          = p
  )
}

# Base adjacency matrix (0/1, p x p).
.sim_base_graph <- function(p, graph_type, edge_prob) {
  adj_base <- matrix(0L, p, p)
  if (p < 2) return(adj_base)
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
        if (stats::runif(1) < edge_prob) adj_base[i, j] <- adj_base[j, i] <- 1L
      }
    }
  } else if (graph_type == "hub") {
    n_hubs <- max(1L, floor(p / 5))
    hub_nodes <- as.integer(round(seq(1, p, length.out = n_hubs + 1)[1:n_hubs]))
    for (h in hub_nodes) {
      spokes <- setdiff(seq_len(p), h)
      # Connect hub to ~40% of other nodes
      connected <- spokes[stats::runif(length(spokes)) < 0.4]
      for (s in connected) adj_base[h, s] <- adj_base[s, h] <- 1L
    }
  }
  adj_base
}

# Flip each pair (i < j) independently with probability prob.
.perturb_graph <- function(adj, prob) {
  p <- nrow(adj)
  if (p < 2) return(adj)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      if (stats::runif(1) < prob) {
        adj[i, j] <- 1L - adj[i, j]
        adj[j, i] <- adj[i, j]
      }
    }
  }
  adj
}

# Group 1 keeps adj_base; groups 2..K perturb it independently.
.sim_group_graphs <- function(adj_base, K, perturb_prob) {
  adj_list <- vector("list", K)
  adj_list[[1]] <- adj_base
  for (k in seq_len(K)[-1]) adj_list[[k]] <- .perturb_graph(adj_base, perturb_prob)
  adj_list
}

# Precision matrices and centred data for given graphs.
.sim_from_graphs <- function(adj_list, n_vec, signal) {
  K <- length(adj_list); p <- nrow(adj_list[[1]])
  Omega_list <- vector("list", K)
  for (k in seq_len(K)) {
    Omega_k <- diag(1, p)
    if (p > 1) {
      for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
          if (adj_list[[k]][i, j] == 1L) {
            val <- stats::runif(1, signal[1], signal[2]) * sample(c(-1, 1), 1)
            Omega_k[i, j] <- Omega_k[j, i] <- val
          }
        }
      }
    }
    # Row-normalize off-diagonals as in Peterson et al. / Danaher et al. so
    # that Omega_k is positive definite
    Omega_list[[k]] <- .make_pd(Omega_k)
  }
  data_list <- vector("list", K)
  S_list    <- vector("list", K)
  for (k in seq_len(K)) {
    Sigma_k <- solve(Omega_list[[k]])
    Sigma_k <- (Sigma_k + t(Sigma_k)) / 2
    L <- chol(Sigma_k)  # upper Cholesky: Sigma = L'L
    Xk <- matrix(stats::rnorm(n_vec[k] * p), nrow = n_vec[k], ncol = p) %*% L
    Xk <- scale(Xk, center = TRUE, scale = FALSE)
    data_list[[k]] <- Xk
    S_list[[k]]    <- crossprod(Xk)
  }
  list(Omega_list = Omega_list, data_list = data_list, S_list = S_list)
}


#' Simulate multi-platform data from multiple Gaussian graphical models
#'
#' Generates S platforms (data types) of data, each with the same K groups,
#' with graphs that share structure across platforms and across groups.
#'
#' The data-generating process is hierarchical:
#' 1. One base graph is drawn on \code{max(p_vec)} nodes (see
#'    \code{graph_type}).
#' 2. Platform \code{s} uses the leading \code{p_vec[s] x p_vec[s]} block of
#'    the base graph. Platform 1 uses it unchanged; platforms 2, ..., S flip
#'    each pair independently with probability \code{platform_perturb_prob}.
#'    Variable \code{j} therefore plays the same role in every platform that
#'    has it.
#' 3. Within each platform, group 1 uses the platform graph and groups
#'    2, ..., K flip each pair with probability \code{perturb_prob}, exactly as
#'    in [simulate_multiggm()]. Group perturbations are independent across
#'    platforms.
#' 4. Precision matrices and data are generated from each graph as in
#'    [simulate_multiggm()].
#'
#' @param K Integer; number of sample groups. Default 2.
#' @param S Integer; number of platforms (data types). Default 2.
#' @param p_vec Integer vector of length S; number of variables per platform.
#'   Default \code{c(20, 15)}.
#' @param n Integer (scalar or length-K vector); sample size per group
#'   (same across platforms). Default 100.
#' @param graph_type Character; type of base graph. See
#'   \code{\link{simulate_multiggm}}.
#' @param edge_prob Numeric; edge probability of the base graph for
#'   \code{graph_type = "random"}. Default 0.1.
#' @param perturb_prob Numeric; probability of edge flip between groups within
#'   a platform. Default 0.1.
#' @param platform_perturb_prob Numeric; probability of edge flip between
#'   platform 1 and each other platform. Default 0.2.
#' @param signal Numeric vector of length 2; magnitude range for off-diagonal
#'   precision entries. Default \code{c(0.3, 0.6)}.
#' @param seed Optional integer random seed.
#'
#' @return A list with components:
#' * \code{platform_data}: List of S elements, each suitable as input
#'       to \code{multiggm_mcmc(method = "ssvs_platform", platform_data = ...)}.
#'       Each element is a list with \code{data_list}, \code{S_list},
#'       \code{n_vec}, \code{Omega_list}, and \code{adj_list}.
#' * \code{adj_base}: The base graph on \code{max(p_vec)} nodes.
#' * \code{platform_adj}: List of S platform graphs (before group
#'       perturbation).
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
                                        edge_prob = 0.1,
                                        perturb_prob = 0.1,
                                        platform_perturb_prob = 0.2,
                                        signal = c(0.3, 0.6),
                                        seed = NULL) {

  graph_type <- match.arg(graph_type)
  if (!is.null(seed)) set.seed(seed)
  p_vec <- as.integer(p_vec)
  if (length(p_vec) != S) stop("p_vec must have length S.")
  n_vec <- rep_len(as.integer(n), K)

  adj_base <- .sim_base_graph(max(p_vec), graph_type, edge_prob)
  platform_adj <- vector("list", S)
  platform_data <- vector("list", S)

  for (s in seq_len(S)) {
    idx <- seq_len(p_vec[s])
    adj_s <- adj_base[idx, idx, drop = FALSE]
    if (s > 1) adj_s <- .perturb_graph(adj_s, platform_perturb_prob)
    platform_adj[[s]] <- adj_s

    adj_list <- .sim_group_graphs(adj_s, K, perturb_prob)
    sim_s <- .sim_from_graphs(adj_list, n_vec, signal)
    platform_data[[s]] <- list(
      data_list  = sim_s$data_list,
      S_list     = sim_s$S_list,
      n_vec      = n_vec,
      Omega_list = sim_s$Omega_list,
      adj_list   = adj_list
    )
  }

  list(
    platform_data = platform_data,
    adj_base      = adj_base,
    platform_adj  = platform_adj,
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
