sim2 <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)

# --- Convergence statistics ---------------------------------------------------

test_that("R-hat and ESS behave as expected on known chains", {
  set.seed(1)
  iid <- matrix(rnorm(4000), 1000, 4)
  v <- multiGGM:::.convergence(iid)
  expect_lt(v[["rhat"]], 1.01)
  expect_gt(v[["ess_bulk"]], 3000); expect_lt(v[["ess_bulk"]], 5200)
  expect_gt(v[["ess_tail"]], 2500)
  # AR(1) with rho = 0.9: ESS about N (1 - rho) / (1 + rho)
  ar <- sapply(1:4, function(i) as.numeric(stats::arima.sim(list(ar = 0.9), 5000)))
  e <- multiGGM:::.ess_basic(multiGGM:::.split_chains(ar))
  expect_gt(e, 0.6 * 20000 * 0.1 / 1.9); expect_lt(e, 1.5 * 20000 * 0.1 / 1.9)
  # chains centred at different values
  shifted <- iid + rep(c(0, 0, 0, 1), each = 1000)
  expect_gt(multiGGM:::.convergence(shifted)[["rhat"]], 1.05)
  # a constant trace gives NA, not an error
  expect_true(is.na(multiGGM:::.convergence(matrix(3, 100, 2))[["rhat"]]))
})

test_that("batch-means MCSE matches a hand computation, single and pooled chains", {
  fit <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 30, nsave = 50,
                       nchains = 2, seed = 4, mcse_batches = 5)
  bm_nvar <- function(x, B) {
    L <- length(x) %/% B
    b <- pmin((seq_along(x) - 1) %/% L, B - 1)
    m_b <- tapply(x, b, mean); len <- tabulate(b + 1)
    length(x) * sum(len * (m_b - mean(x))^2) / (B - 1)
  }
  G <- lapply(fit$chains, function(ch) ch$adj_save[1, 2, 1, ] != ch$adj_save[1, 2, 2, ])
  ep1 <- edge_probs(fit, chain = 1)
  row <- ep1$type == "pair" & ep1$i == 1 & ep1$j == 2
  expect_equal(ep1$mcse_struct[row], sqrt(bm_nvar(G[[1]], 5)) / 50)
  ep <- edge_probs(fit)
  expect_equal(ep$mcse_struct[row], sqrt(bm_nvar(G[[1]], 5) + bm_nvar(G[[2]], 5)) / 100)
  # without batches the column is NA
  f0 <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 10, nsave = 10,
                      seed = 4, mcse_batches = 0)
  expect_true(all(is.na(edge_probs(f0)$mcse_struct)))
})

test_that("mcmc_diagnostics() works without draws and flags disagreeing chains", {
  fit <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 100, nsave = 200,
                       nchains = 2, seed = 1, save_draws = FALSE)
  d <- mcmc_diagnostics(fit)
  expect_s3_class(d, "multiggm_diagnostics")
  expect_true(all(c("n_edges[1]", "loglik[2]", "theta[1,2]", "n_diff_struct[1,2]") %in%
                    d$scalars$quantity))
  expect_true(all(c("rhat", "ess_bulk", "ess_tail", "mcse_mean", "flag") %in% names(d$scalars)))
  expect_true(all(is.finite(d$edges$max_chain_diff)))
  expect_output(print(d), "Scalar summaries")
  expect_s3_class(plot_chain_agreement(fit), "ggplot")
  expect_equal(dim(edge_counts(fit$chains[[1]])), c(200L, 2L))
  for (w in c("theta", "edges", "loglik", "diff"))
    expect_s3_class(plot_trace(fit, w), "ggplot")
  expect_s3_class(plot_trace(fit, "diff", rule = "union", c = 0.05), "ggplot")
  m <- multiggm_mcmc(data_list = sim2$data_list, method = "mghs", burnin = 50, nsave = 50,
                     seed = 1, save_draws = FALSE)
  dm <- mcmc_diagnostics(m)
  expect_true(all(c("tau[1]", "R[1,2]") %in% dm$scalars$quantity))
  expect_null(dm$acceptance)
  expect_s3_class(plot_trace(m, "edges"), "ggplot")
  expect_s3_class(plot_trace(m, "diff"), "ggplot")
})

# --- Directional, node and network summaries -------------------------------

test_that("directional, node-level and count summaries match the draws", {
  f <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 30, nsave = 60,
                     seed = 2, diff_deltas = c(0.05, 0.1))
  pc <- posterior_pcor(f)
  D <- pc[, , 1, ] - pc[, , 2, ]
  ep <- edge_probs(f, scale = "pcor")
  row <- ep$type == "pair" & ep$i == 2 & ep$j == 5
  expect_equal(ep$p_up_0.1[row], mean(D[2, 5, ] > 0.1))
  expect_equal(ep$p_down_0.1[row], mean(D[2, 5, ] < -0.1))
  # direction flips with the pair order
  up21 <- edge_prob_matrix(f, "p_up", pair = c(2, 1), c = 0.1, scale = "pcor")
  expect_equal(up21[2, 5], mean(D[2, 5, ] < -0.1))

  G <- f$adj_save[, , 1, ] != f$adj_save[, , 2, ]
  deg_u <- apply(G | abs(D) > 0.05, 3, function(m) { diag(m) <- FALSE; rowSums(m) })
  nd <- node_diff_summary(f, rule = "union", c = 0.05)
  expect_equal(nd$diff_degree, rowMeans(deg_u))
  expect_equal(nd$p_any, rowMeans(deg_u > 0))
  rw <- apply(abs(D), 3, function(m) { diag(m) <- 0; rowSums(m) })
  expect_equal(nd$rewiring, rowMeans(rw))
  expect_error(node_diff_summary(f, rule = "magnitude"), "needs a single threshold")

  dc <- diff_count_posterior(f, rule = "magnitude", c = 0.1)
  cnt <- apply(abs(D) > 0.1, 3, function(m) sum(m[upper.tri(m)]))
  expect_equal(dc$draws, cnt)
  expect_equal(dc$p_any, mean(cnt > 0))
  expect_output(print(dc), "P\\(at least one")
  # off-grid c from saved draws
  expect_equal(diff_count_posterior(f, rule = "magnitude", c = 0.07)$draws,
               apply(abs(D) > 0.07, 3, function(m) sum(m[upper.tri(m)])))
})

# --- q-values -----------------------------------------------------------------

test_that("bfdr_qvalues() reproduces bayes_fdr() at every alpha, with ties", {
  set.seed(3)
  p <- c(runif(40), rep(0.9, 5), rep(0.5, 3), 1, 1)
  q <- bfdr_qvalues(p)
  for (a in c(0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.4))
    expect_identical(q <= a, as.logical(bayes_fdr(p, alpha = a)), label = paste("alpha", a))
  expect_true(all(tapply(q, p, function(z) length(unique(z))) == 1))
  m <- matrix(0, 4, 4); m[upper.tri(m)] <- c(0.99, 0.2, 0.95, 0.5, 0.7, 0.1); m <- m + t(m)
  qm <- bfdr_qvalues(m)
  expect_true(isSymmetric(qm))
  expect_true(all(is.na(diag(qm))))
  f <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 30, nsave = 40, seed = 2)
  s <- select_diff_edges(f, rule = "union", c = 0.05, alpha = 0.2)
  expect_true("q_value" %in% names(s$edges))
  expect_true(all(s$edges$q_value <= 0.2))
})

# --- WAIC ---------------------------------------------------------------------

test_that("compute_waic() equals the WAIC from the pointwise log-likelihood draws", {
  f <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 30, nsave = 40,
                     nchains = 2, seed = 7, waic = TRUE)
  ll <- do.call(rbind, lapply(1:2, function(k) {
    X <- scale(sim2$data_list[[k]], scale = FALSE)
    draws <- do.call(c, lapply(f$chains, function(ch) lapply(seq_len(40), function(s) ch$C_save[, , k, s])))
    sapply(draws, function(O) {
      U <- chol(O)
      sum(log(diag(U))) - 0.5 * ncol(X) * log(2 * pi) - 0.5 * colSums((U %*% t(X))^2)
    })
  }))
  lppd <- apply(ll, 1, function(l) log(mean(exp(l - max(l)))) + max(l))
  pw <- apply(ll, 1, var)
  w <- suppressWarnings(compute_waic(f))
  expect_equal(w$estimates["waic", "estimate"], -2 * sum(lppd - pw), tolerance = 1e-8)
  expect_equal(w$estimates["p_waic", "estimate"], sum(pw), tolerance = 1e-8)
  expect_equal(nrow(w$pointwise), 120L)
  expect_error(compute_waic(multiggm_mcmc(data_list = sim2$data_list, method = "ssvs",
                                          burnin = 5, nsave = 5)), "waic = TRUE")
  expect_error(multiggm_mcmc(S_list = sim2$S_list, n_vec = sim2$n_vec, method = "ssvs",
                             burnin = 5, nsave = 5, waic = TRUE), "data_list")
})

# --- Prior-only fits and Bayes factors -------------------------------------

test_that("prior_only fits ignore the data and give Bayes factors", {
  other <- simulate_multiggm(K = 2, p = 6, n = 30, seed = 9)
  for (m in c("ssvs", "gwishart")) {
    a <- multiggm_mcmc(data_list = sim2$data_list, method = m, burnin = 20, nsave = 30,
                       seed = 1, prior_only = TRUE, save_draws = FALSE)
    b <- multiggm_mcmc(data_list = other$data_list, method = m, burnin = 20, nsave = 30,
                       seed = 1, prior_only = TRUE, save_draws = FALSE)
    expect_identical(a$edge_acc, b$edge_acc, label = m)
    expect_true(a$prior_only)
    expect_null(a$trace$loglik)
  }
  expect_error(multiggm_mcmc(data_list = sim2$data_list, method = "mghs", burnin = 5,
                             nsave = 5, prior_only = TRUE), "improper")
  fit <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 50, nsave = 100,
                       seed = 1, save_draws = FALSE)
  pri <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 50, nsave = 100,
                       seed = 2, prior_only = TRUE, save_draws = FALSE)
  bf <- edge_bayes_factors(fit, pri)
  expect_equal(nrow(bf), 3 * 15)
  expect_true(all(is.finite(bf$bayes_factor)))
  expect_equal(length(unique(bf$prior[bf$type == "group" & bf$group1 == 1])), 1L)
  expect_error(edge_bayes_factors(fit, fit), "prior_only")
})

# --- FDR calibration ---------------------------------------------------------

test_that("fdr_calibration() reports expected and achieved FDR", {
  sim <- simulate_multiggm(K = 2, p = 8, n = 200, seed = 3)
  f <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs", burnin = 100, nsave = 100,
                     seed = 3, save_draws = FALSE)
  cal <- fdr_calibration(f, sim, rule = "union", c = 0.1, alpha = c(0.05, 0.1, 0.2))
  expect_named(cal, c("alpha", "n_selected", "expected_fdr", "achieved_fdr", "power", "n_true"))
  expect_true(all(diff(cal$n_selected) >= 0))
  expect_true(all(cal$expected_fdr <= cal$alpha + 1e-12))
  cg <- fdr_calibration(f, sim$Omega_list, group = 1, rule = "structural")
  O <- sim$Omega_list[[1]]
  expect_equal(cg$n_true[1], sum(abs(O[upper.tri(O)]) > 1e-8))
})
