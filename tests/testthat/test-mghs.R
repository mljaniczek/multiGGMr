sim <- simulate_multiggm(K = 2, p = 7, n = 120, seed = 5)

test_that("mghs fit has the multiggm_fit layout and no indicators", {
  fit <- multiggm_mcmc(data_list = sim$data_list, method = "mghs",
                       burnin = 100, nsave = 60, thin = 2, seed = 1)
  expect_s3_class(fit, "multiggm_fit")
  expect_identical(fit$method, "mghs")
  expect_equal(dim(fit$C_save), c(7, 7, 2, 60))
  expect_equal(dim(fit$R_save), c(2, 2, 60))
  expect_equal(fit$n_post, 120L)
  expect_null(fit$adj_save)
  expect_true(all(abs(fit$C_save[, , 1, 5] - t(fit$C_save[, , 1, 5])) < 1e-12))

  expect_error(pip_edges(fit), "edge-inclusion indicators")
  expect_error(pip_diff_edge(fit), "edge-inclusion indicators")
  expect_error(edge_counts(fit), "edge-inclusion indicators")
  # since 0.5.0 the structural rule uses the shrinkage weight: P(S_1 != S_2)
  expect_s3_class(select_diff_edges(fit, rule = "structural"), "multiggm_diff_selection")
  expect_error(posterior_pcor(fit, use = "Sig_save"), "no Sig_save draws")

  expect_output(print(fit), "mGHS")
  expect_output(print(summary(fit)), "Posterior mean of R")
  expect_s3_class(plot_trace(fit, "theta"), "ggplot")
  expect_equal(dim(posterior_pcor(fit)), c(7, 7, 2, 60))
  expect_s3_class(select_diff_edges(fit, alpha = 0.1, delta = 0.1),
                  "multiggm_diff_selection")
})

test_that("running summaries equal draw-based summaries when thin = 1", {
  args <- list(data_list = sim$data_list, method = "mghs", burnin = 80,
               nsave = 100, seed = 2, diff_deltas = c(0.02, 0.1))
  full <- do.call(multiggm_mcmc, args)
  lean <- do.call(multiggm_mcmc, c(args, save_draws = FALSE))
  expect_null(lean$C_save)
  expect_equal(lean$C_mean, full$C_mean)

  for (sc in c("pcor", "precision")) {
    a <- diff_summary(full, scale = sc, deltas = c(0.02, 0.1), probs = NULL)
    b <- diff_summary(lean, scale = sc, deltas = c(0.02, 0.1), probs = NULL)
    expect_equal(b, a, tolerance = 1e-10)
  }
  expect_equal(coef(lean), coef(full), tolerance = 1e-10)
  expect_equal(fitted(lean), fitted(full), tolerance = 1e-10)

  rev <- diff_summary(lean, groups = c(2, 1), deltas = 0.1, probs = NULL)
  fwd <- diff_summary(lean, groups = c(1, 2), deltas = 0.1, probs = NULL)
  expect_equal(rev$mean, -fwd$mean)
  expect_equal(rev$p_exceed_0.1, fwd$p_exceed_0.1)

  expect_message(diff_summary(lean, deltas = 0.1), "need saved draws")
  expect_error(diff_summary(lean, deltas = 0.07), "grid set at fit time")
  expect_error(diff_draws(lean), "save_draws = FALSE")
  expect_error(posterior_precision(lean, summary_fun = median), "only the posterior mean")
})

test_that("summary-only chains are pooled exactly", {
  fit <- multiggm_mcmc(data_list = sim$data_list, method = "mghs", burnin = 60,
                       nsave = 40, nchains = 2, seed = 3, save_draws = FALSE)
  pooled <- diff_summary(fit, deltas = 0.05, probs = NULL)
  one <- lapply(1:2, function(ch) diff_summary(fit, deltas = 0.05, probs = NULL, chain = ch))
  expect_equal(pooled$mean, (one[[1]]$mean + one[[2]]$mean) / 2)
  expect_equal(pooled$p_exceed_0.05, (one[[1]]$p_exceed_0.05 + one[[2]]$p_exceed_0.05) / 2)
  expect_identical(attr(pooled, "n_draws"), 80)
})

test_that("non-finite input stops instead of hanging", {
  S <- list(diag(5), diag(5)); S[[1]][2, 3] <- S[[1]][3, 2] <- NaN
  r <- mghs_engine_cpp(10L, 10L, 1L, c(20, 20), S, 5L, c(2, 1.3), TRUE, 0.05, FALSE, 0L)
  expect_identical(r$stuck, 1L)
  expect_error(multiggm_mcmc(S_list = S, n_vec = c(20, 20), method = "mghs",
                             burnin = 5, nsave = 5), "non-finite")
  expect_error(multiggm_mcmc(data_list = sim$data_list, method = "mghs", burnin = 5,
                             nsave = 5, hyper = list(t_var = 0.5)), "t_var")
})
