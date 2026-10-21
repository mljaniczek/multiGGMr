#' Plot method for multiggm_fit objects
#'
#' Produces diagnostic and summary plots for a fitted multi-GGM model.
#'
#' @param x A \code{multiggm_fit} or \code{multiggm_fit_list} object returned
#'   by \code{\link{multiggm_mcmc}}. For several chains, trace plots overlay
#'   the chains and the PIP plots pool them.
#' @param type Character string specifying the plot type:
#'   \describe{
#'     \item{\code{"trace_theta"}}{Trace plots of the graph similarity
#'       parameter \eqn{\theta_{km}} across post-burn-in iterations. Values
#'       > 0 indicate the model is borrowing strength between groups.
#'       Requires K >= 2.}
#'     \item{\code{"trace_edges"}}{Trace of the number of included edges per
#'       group across iterations. Useful for assessing convergence and
#'       stationarity of model complexity.}
#'     \item{\code{"pip"}}{Heatmap of posterior inclusion probabilities (PIP)
#'       for each edge, faceted by group. Uses viridis color scale.}
#'     \item{\code{"network"}}{Network visualization of the estimated graph
#'       at the given PIP threshold, using igraph. One plot per group.}
#'   }
#' @param pip_threshold Numeric; PIP threshold for edge selection in
#'   network plots. Default 0.5. Ignored for other plot types.
#' @param ... Additional arguments (currently ignored).
#'
#' @return For \code{"trace_theta"}, \code{"trace_edges"}, and \code{"pip"}:
#'   a \code{ggplot} object (can be further customized). For
#'   \code{"network"}: invisible \code{NULL} (plots are drawn as a side
#'   effect using igraph).
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' plot(fit, type = "trace_edges")
#' plot(fit, type = "pip")
#'
#' @seealso [plot_trace()], [plot_pip_heatmap()], [plot_network()],
#'   [plot_roc()], [plot_recovery()]
#' @export
plot.multiggm_fit <- function(x, type = c("trace_theta", "trace_edges", "pip", "network"),
                               pip_threshold = 0.5, ...) {
  .plot_fit(x, match.arg(type), pip_threshold)
}

#' @rdname plot.multiggm_fit
#' @export
plot.multiggm_fit_list <- function(x, type = c("trace_theta", "trace_edges", "pip", "network"),
                                    pip_threshold = 0.5, ...) {
  .plot_fit(x, match.arg(type), pip_threshold)
}

.plot_fit <- function(x, type, pip_threshold) {
  switch(type,
    trace_theta = plot_trace(x, what = "theta"),
    trace_edges = plot_trace(x, what = "edges"),
    pip         = plot_pip_heatmap(pip_edges(x)),
    network     = {
      pip <- pip_edges(x)
      K <- x$K
      for (k in seq_len(K)) {
        adj_k <- (pip[, , k] >= pip_threshold) * 1L
        diag(adj_k) <- 0L
        plot_network(adj_k, main = paste("Group", k))
      }
      invisible(NULL)
    }
  )
}

#' Trace plots for MCMC diagnostics
#'
#' Creates trace plots to assess MCMC convergence. Supports plotting either
#' the graph similarity parameter (theta) or the number of edges per group
#' across saved iterations.
#'
#' @param fit A \code{multiggm_fit} or \code{multiggm_fit_list} object
#'   returned by \code{\link{multiggm_mcmc}}; chains are overlaid in colour.
#' @param what Character; what to plot:
#'   \describe{
#'     \item{\code{"theta"}}{Graph similarity parameters \eqn{\theta_{km}}
#'       for each pair of groups. Requires K >= 2.}
#'     \item{\code{"edges"}}{Number of edges per group at each saved
#'       iteration. Color-coded by group.}
#'   }
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' plot_trace(fit, "theta")
#' plot_trace(fit, "edges")
#'
#' @seealso [plot.multiggm_fit()], [edge_counts()]
#' @export
plot_trace <- function(fit, what = c("theta", "edges")) {
  what <- match.arg(what)
  chains <- .fit_chains(fit)
  multi <- length(chains) > 1L
  K <- chains[[1]]$K

  if (what == "theta") {
    if (K < 2) stop("Theta trace requires K >= 2.")
    is_mghs <- identical(chains[[1]]$method, "mghs")
    nm <- if (is_mghs) "R" else "theta"
    df <- do.call(rbind, lapply(seq_along(chains), function(ch) {
      S <- if (is_mghs) chains[[ch]]$R_save else chains[[ch]]$Theta_save
      do.call(rbind, lapply(seq_len(K - 1), function(k) do.call(rbind, lapply((k + 1):K, function(m)
        data.frame(iteration = seq_len(dim(S)[3]), value = S[k, m, ],
                   pair = sprintf("%s[%d,%d]", nm, k, m), chain = paste("Chain", ch),
                   stringsAsFactors = FALSE)))))
    }))
    aes <- if (multi) ggplot2::aes(x = .data$iteration, y = .data$value, color = .data$chain)
           else ggplot2::aes(x = .data$iteration, y = .data$value)
    ggplot2::ggplot(df, aes) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::facet_wrap(~ pair, scales = "free_y") +
      ggplot2::labs(title = if (is_mghs) "Trace: correlation of edge scales (R)"
                            else "Trace: graph similarity (theta)",
                    x = "Iteration (post-burnin)", y = nm, color = "Chain") +
      ggplot2::theme_minimal()
  } else {
    df <- do.call(rbind, lapply(seq_along(chains), function(ch) {
      ec <- edge_counts(chains[[ch]])
      do.call(rbind, lapply(seq_len(ncol(ec)), function(k)
        data.frame(iteration = seq_len(nrow(ec)), edges = ec[, k],
                   group = paste("Group", k), chain = paste("Chain", ch),
                   stringsAsFactors = FALSE)))
    }))
    if (multi) {
      ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$edges, color = .data$chain)) +
        ggplot2::geom_line(alpha = 0.6) +
        ggplot2::facet_wrap(~ group, scales = "free_y") +
        ggplot2::labs(title = "Trace: edge count per group",
                      x = "Iteration (post-burnin)", y = "Number of edges", color = "Chain") +
        ggplot2::theme_minimal()
    } else {
      ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$edges, color = .data$group)) +
        ggplot2::geom_line(alpha = 0.6) +
        ggplot2::labs(title = "Trace: edge count per group",
                      x = "Iteration (post-burnin)", y = "Number of edges", color = "Group") +
        ggplot2::theme_minimal()
    }
  }
}

#' Heatmap of posterior inclusion probabilities
#'
#' Creates a heatmap of posterior inclusion probabilities (PIP) using
#' ggplot2 with a viridis color scale. Only the upper triangle is shown.
#'
#' @param pip A numeric array of PIPs: either \code{[p, p]} (single group)
#'   or \code{[p, p, K]} (multiple groups, as returned by
#'   \code{\link{pip_edges}}).
#' @param node_labels Optional character vector of length p with node names.
#'   If \code{NULL}, nodes are labeled 1 through p.
#'
#' @return A \code{ggplot} object. For K > 1, faceted by group.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' pip <- pip_edges(fit)
#' plot_pip_heatmap(pip)
#'
#' @seealso [pip_edges()], [plot.multiggm_fit()]
#' @export
plot_pip_heatmap <- function(pip, node_labels = NULL) {
  dims <- dim(pip)

  if (length(dims) == 2) {
    # Single matrix
    pip <- array(pip, dim = c(dims, 1))
    dims <- dim(pip)
  }

  p <- dims[1]
  K <- dims[3]
  if (is.null(node_labels)) node_labels <- as.character(seq_len(p))

  df_list <- vector("list", K)
  for (k in seq_len(K)) {
    mat <- pip[, , k]
    # Only upper triangle
    mat[lower.tri(mat, diag = TRUE)] <- NA
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (!is.na(mat[i, j])) {
          df_list[[length(df_list) + 1]] <- data.frame(
            row = node_labels[i], col = node_labels[j],
            pip = mat[i, j],
            group = paste("Group", k),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  df <- do.call(rbind, df_list)

  # Factor levels to preserve order
  df$row <- factor(df$row, levels = rev(node_labels))
  df$col <- factor(df$col, levels = node_labels)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data$col, y = .data$row, fill = .data$pip)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), name = "PIP") +
    ggplot2::labs(title = "Posterior Inclusion Probabilities", x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (K > 1) g <- g + ggplot2::facet_wrap(~ group)
  g
}

#' Network plot of an adjacency matrix
#'
#' Visualizes a graph using igraph with Fruchterman-Reingold layout.
#'
#' @param adj A symmetric adjacency matrix (0/1 or weighted). Non-zero
#'   entries are treated as edges.
#' @param layout An igraph layout function or a two-column matrix of node
#'   coordinates. Default: Fruchterman-Reingold (\code{layout_with_fr}).
#' @param node_labels Optional character vector of node names. If
#'   \code{NULL}, nodes are labeled 1 through p.
#' @param main Character; plot title. Default \code{""}.
#' @param vertex_size Numeric; node size. Default 15.
#' @param ... Additional arguments passed to \code{igraph::plot.igraph()}.
#'
#' @return Invisible \code{NULL}. The plot is drawn as a side effect.
#'
#' @examples
#' adj <- matrix(0, 5, 5)
#' adj[1,2] <- adj[2,1] <- adj[2,3] <- adj[3,2] <- 1
#' plot_network(adj, main = "Simple graph")
#'
#' @seealso [plot.multiggm_fit()], [pip_edges()]
#' @export
plot_network <- function(adj, layout = NULL, node_labels = NULL,
                          main = "", vertex_size = 15, ...) {
  adj <- as.matrix(adj)
  diag(adj) <- 0
  # Ensure symmetric
  adj <- ((adj + t(adj)) > 0) * 1

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)

  if (is.null(node_labels)) {
    node_labels <- as.character(seq_len(nrow(adj)))
  }
  igraph::V(g)$label <- node_labels

  if (is.null(layout)) {
    layout <- igraph::layout_with_fr(g)
  }

  igraph::plot.igraph(g, layout = layout,
                       vertex.size = vertex_size,
                       vertex.color = "steelblue",
                       vertex.label.color = "white",
                       vertex.label.cex = 0.7,
                       edge.color = "gray40",
                       main = main, ...)
  invisible(NULL)
}

#' Compare estimated graph to ground truth
#'
#' Shows side-by-side heatmaps of the true adjacency matrix, PIP, and
#' thresholded estimate. Returns confusion metrics invisibly.
#'
#' @param fit A \code{multiggm_fit} object returned by
#'   \code{\link{multiggm_mcmc}}.
#' @param true_adj A list of K true adjacency matrices (0/1), or a single
#'   matrix (recycled for all groups).
#' @param pip_threshold Numeric; threshold for edge selection. Default 0.5.
#' @param groups Integer vector; which groups to plot. Default: all groups.
#'
#' @return A named list of confusion metric vectors per group (invisible).
#'   Each element is a named numeric vector as returned by
#'   \code{\link{confusion_at_threshold}}, with components \code{TP},
#'   \code{FP}, \code{TN}, \code{FN}, \code{TPR}, \code{FPR}.
#'
#' @examples
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 200, nsave = 100)
#' cm <- plot_recovery(fit, sim$adj_list)
#' cm$Group_1  # confusion metrics for group 1
#'
#' @seealso [confusion_at_threshold()], [roc_auc()], [plot_roc()]
#' @export
plot_recovery <- function(fit, true_adj, pip_threshold = 0.5, groups = NULL) {
  pip <- pip_edges(fit)
  K <- fit$K
  p <- fit$p

  if (!is.list(true_adj)) true_adj <- rep(list(true_adj), K)
  if (is.null(groups)) groups <- seq_len(K)

  node_labels <- as.character(seq_len(p))
  confusion_list <- list()

  df_all <- data.frame()
  for (k in groups) {
    truth_k <- as.matrix(true_adj[[k]])
    pip_k <- pip[, , k]
    est_k <- (pip_k >= pip_threshold) * 1L

    # Confusion metrics
    cm <- confusion_at_threshold(pip_k, truth_k, pip_threshold)
    confusion_list[[paste0("Group_", k)]] <- cm

    # Build data frames for three panels
    for (panel_name in c("True graph", "PIP", "Estimated")) {
      mat <- switch(panel_name,
        "True graph" = truth_k,
        "PIP"        = pip_k,
        "Estimated"  = est_k
      )
      mat[lower.tri(mat, diag = TRUE)] <- NA
      for (i in seq_len(p)) {
        for (j in seq_len(p)) {
          if (!is.na(mat[i, j])) {
            df_all <- rbind(df_all, data.frame(
              row = node_labels[i], col = node_labels[j],
              value = mat[i, j],
              panel = panel_name,
              group = paste("Group", k),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }

  df_all$row <- factor(df_all$row, levels = rev(node_labels))
  df_all$col <- factor(df_all$col, levels = node_labels)
  df_all$panel <- factor(df_all$panel, levels = c("True graph", "PIP", "Estimated"))

  g <- ggplot2::ggplot(df_all, ggplot2::aes(x = .data$col, y = .data$row, fill = .data$value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1), name = "Value") +
    ggplot2::labs(title = sprintf("Graph recovery (threshold = %.2f)", pip_threshold),
                  x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (length(groups) > 1) {
    g <- g + ggplot2::facet_grid(group ~ panel)
  } else {
    g <- g + ggplot2::facet_wrap(~ panel)
  }

  print(g)
  invisible(confusion_list)
}
