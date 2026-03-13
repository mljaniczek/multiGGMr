#' Plot method for multiggm_fit objects
#'
#' Produces diagnostic and summary plots for a fitted multi-GGM model.
#'
#' @param x A \code{multiggm_fit} object.
#' @param type Character string specifying the plot type:
#'   \describe{
#'     \item{\code{"trace_theta"}}{Trace plots of theta (graph similarity) parameters.}
#'     \item{\code{"trace_edges"}}{Trace of edge count per group across iterations.}
#'     \item{\code{"pip"}}{Heatmap of posterior inclusion probabilities.}
#'     \item{\code{"network"}}{Network plot of estimated graph at given PIP threshold.}
#'   }
#' @param pip_threshold PIP threshold for network plot. Default 0.5.
#' @param ... Additional arguments (currently ignored).
#' @return A ggplot object (for trace/pip/roc) or NULL (for network, which uses igraph).
#' @export
plot.multiggm_fit <- function(x, type = c("trace_theta", "trace_edges", "pip", "network"),
                               pip_threshold = 0.5, ...) {
  type <- match.arg(type)
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
#' @param fit A \code{multiggm_fit} object.
#' @param what Character: \code{"theta"} for graph similarity parameters,
#'   \code{"edges"} for edge counts per group.
#' @return A ggplot object.
#' @export
plot_trace <- function(fit, what = c("theta", "edges")) {
  what <- match.arg(what)

  if (what == "theta") {
    K <- fit$K
    nsave <- dim(fit$Theta_save)[3]
    if (K < 2) stop("Theta trace requires K >= 2.")

    pairs <- list()
    for (k in 1:(K - 1)) {
      for (m in (k + 1):K) {
        pairs[[length(pairs) + 1]] <- list(k = k, m = m)
      }
    }

    df_list <- vector("list", length(pairs))
    for (i in seq_along(pairs)) {
      k <- pairs[[i]]$k
      m <- pairs[[i]]$m
      df_list[[i]] <- data.frame(
        iteration = seq_len(nsave),
        value = fit$Theta_save[k, m, ],
        pair = sprintf("theta[%d,%d]", k, m),
        stringsAsFactors = FALSE
      )
    }
    df <- do.call(rbind, df_list)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value)) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::facet_wrap(~ pair, scales = "free_y") +
      ggplot2::labs(title = "Trace: graph similarity (theta)",
                    x = "Iteration (post-burnin)", y = "theta") +
      ggplot2::theme_minimal()
    return(p)

  } else {
    # Edge counts
    ec <- edge_counts(fit)
    K <- ncol(ec)
    nsave <- nrow(ec)

    df_list <- vector("list", K)
    for (k in seq_len(K)) {
      df_list[[k]] <- data.frame(
        iteration = seq_len(nsave),
        edges = ec[, k],
        group = paste("Group", k),
        stringsAsFactors = FALSE
      )
    }
    df <- do.call(rbind, df_list)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$edges, color = .data$group)) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::labs(title = "Trace: edge count per group",
                    x = "Iteration (post-burnin)", y = "Number of edges",
                    color = "Group") +
      ggplot2::theme_minimal()
    return(p)
  }
}

#' Heatmap of posterior inclusion probabilities
#'
#' @param pip An array of PIPs: either \code{[p, p]} (single group) or
#'   \code{[p, p, K]} (multiple groups from \code{pip_edges()}).
#' @param node_labels Optional character vector of node names.
#' @return A ggplot object.
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
#' Uses igraph to visualize a graph.
#'
#' @param adj A symmetric adjacency matrix (0/1 or weighted).
#' @param layout igraph layout function or matrix. Default: Fruchterman-Reingold.
#' @param node_labels Optional character vector of node names.
#' @param main Plot title.
#' @param vertex_size Node size. Default 15.
#' @param ... Additional arguments passed to \code{igraph::plot.igraph()}.
#' @return Invisible NULL. Plot is drawn as a side effect.
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

#' Plot ROC curve from roc_auc output
#'
#' @param roc_obj Output from \code{roc_auc()}: a list with FPR, TPR, auc.
#' @param main Optional title.
#' @return A ggplot object.
#' @export
plot_roc <- function(roc_obj, main = "ROC Curve") {
  df <- data.frame(FPR = roc_obj$FPR, TPR = roc_obj$TPR)
  auc_label <- sprintf("AUC = %.3f", roc_obj$auc)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$FPR, y = .data$TPR)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::annotate("text", x = 0.7, y = 0.2, label = auc_label, size = 5) +
    ggplot2::labs(title = main, x = "False Positive Rate", y = "True Positive Rate") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal()
}

#' Compare estimated graph to ground truth
#'
#' Shows side-by-side heatmaps of true adjacency, PIP, and thresholded
#' estimate. Returns confusion metrics invisibly.
#'
#' @param fit A \code{multiggm_fit} object.
#' @param true_adj A list of K true adjacency matrices, or a single matrix
#'   (recycled for all groups).
#' @param pip_threshold Threshold for edge selection. Default 0.5.
#' @param groups Which groups to plot. Default: all.
#' @return A list of confusion metrics per group (invisible).
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
