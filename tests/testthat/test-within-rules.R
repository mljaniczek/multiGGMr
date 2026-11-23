sim <- simulate_multiggm(K = 2, p = 8, n = 300, seed = 2)

test_that("group sign probabilities and the interval rule match the draws", {
  f <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs", burnin = 50, nsave = 80, seed = 1)
  pc <- posterior_pcor(f)
  ep <- edge_probs(f, scale = "pcor")
  row <- ep$type == "group" & ep$group1 == 2 & ep$i == 3 & ep$j == 5
  x <- pc[3, 5, 2, ]
  expect_equal(ep$p_sign[row], max(mean(x > 0), mean(x < 0)))
  s <- select_edges(f, group = 2, rule = "interval", scale = "pcor")
  ps <- edge_prob_matrix(f, "p_sign", group = 2, scale = "pcor")
  expect_identical(s$selected[upper.tri(ps)], ps[upper.tri(ps)] > 0.75)
  expect_true(is.na(s$expected_fdr))
  s9 <- select_edges(f, group = 2, rule = "interval", level = 0.9, scale = "pcor")
  expect_identical(s9$selected[upper.tri(ps)], ps[upper.tri(ps)] > 0.95)
  expect_output(print(s), "credible interval excludes 0")
  expect_error(select_edges(f, rule = "interval", level = 1), "level")
  expect_error(select_edges(f, rule = "shrinkage"), "mGHS")
})

test_that("mGHS shrinkage 1 - kappa is in [0, 1] and ranks true edges higher", {
  skip_on_cran()
  f <- multiggm_mcmc(data_list = sim$data_list, method = "mghs", burnin = 300, nsave = 300,
                     seed = 1, hyper = list(R_eta = 10))
  ep <- edge_probs(f, scale = "pcor")
  g <- ep[ep$type == "group", ]
  expect_true(all(g$mean_shrinkage >= 0 & g$mean_shrinkage <= 1))
  expect_true(all(g$p_shrinkage >= 0 & g$p_shrinkage <= 1))
  truth <- unlist(lapply(sim$Omega_list, function(O) abs(O[upper.tri(O)]) > 1e-8))
  expect_gt(mean(g$mean_shrinkage[truth]), mean(g$mean_shrinkage[!truth]) + 0.2)
  s <- select_edges(f, group = 1, rule = "shrinkage")
  expect_identical(s$decision, "median")
  m <- edge_prob_matrix(f, "mean_shrinkage", group = 1)
  expect_identical(s$selected[upper.tri(m)], m[upper.tri(m)] > 0.5)
  expect_identical(select_edges(f, group = 1, rule = "shrinkage", decision = "bfdr")$decision,
                   "bfdr")
  # pairs have no shrinkage summary
  expect_true(all(is.na(ep$mean_shrinkage[ep$type == "pair"])))
})

test_that("the shrinkage summary survives thresholds computed from draws", {
  skip_on_cran()
  f <- multiggm_mcmc(data_list = sim$data_list, method = "mghs", burnin = 50, nsave = 50,
                     seed = 3, save_weights = TRUE, hyper = list(R_eta = 10))
  a <- edge_probs(f)
  b <- edge_probs(f, c = 0.07)
  expect_equal(b$mean_shrinkage, a$mean_shrinkage)
  expect_false(anyNA(b$mean_shrinkage[b$type == "group"]))
})

test_that("mask_pcor() zeroes unselected edges", {
  f <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs", burnin = 50, nsave = 80,
                     seed = 1, save_draws = FALSE, diff_deltas = c(0.05, 0.1, 0.2))
  m <- mask_pcor(f, c = c(0.1, 0.2))
  P <- f$pcor_mean
  for (k in 1:2) {
    s <- select_edges(f, group = k, rule = "magnitude", c = c(0.1, 0.2)[k], scale = "pcor")
    expect_equal(m[, , k][s$selected], P[, , k][s$selected])
    expect_true(all(m[, , k][!s$selected & upper.tri(P[, , k])] == 0))
    expect_equal(diag(m[, , k]), rep(1, 8))
    expect_true(isSymmetric(m[, , k]))
  }
  expect_identical(attr(m, "c"), c(0.1, 0.2))
  ms <- mask_pcor(f, rule = "structural", decision = "median")
  expect_identical(attr(ms, "n_selected")[1],
                   select_edges(f, group = 1, decision = "median")$n_selected)
  only2 <- mask_pcor(f, c = 0.1, groups = 2)
  expect_equal(only2[, , 1], P[, , 1])
})
