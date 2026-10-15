test_that("simulate_multiggm works with a single group", {
  sim <- simulate_multiggm(K = 1, p = 6, n = 30, seed = 1)
  expect_length(sim$adj_list, 1)
  expect_length(sim$Omega_list, 1)
  expect_length(sim$data_list, 1)
  expect_equal(dim(sim$data_list[[1]]), c(30, 6))
})

test_that("platforms share one base graph", {
  sim <- simulate_multiggm_platform(K = 2, S = 3, p_vec = c(12, 9, 12), n = 40,
                                    graph_type = "random", edge_prob = 0.3,
                                    perturb_prob = 0, platform_perturb_prob = 0,
                                    seed = 3)
  base <- sim$adj_base
  expect_equal(dim(base), c(12, 12))
  expect_true(sum(base) > 0)
  for (s in 1:3) {
    idx <- seq_len(sim$p_vec[s])
    expect_identical(sim$platform_adj[[s]], base[idx, idx])
    # no group perturbation: every group graph equals the platform graph
    for (k in 1:2) expect_identical(sim$platform_data[[s]]$adj_list[[k]], base[idx, idx])
  }
})

test_that("platform and group perturbations act on the shared base", {
  sim <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(30, 30), n = 40,
                                    graph_type = "hub", perturb_prob = 0.05,
                                    platform_perturb_prob = 0.1, seed = 5)
  ut <- upper.tri(sim$adj_base)
  expect_identical(sim$platform_adj[[1]], sim$adj_base)
  flips_platform <- mean(sim$platform_adj[[2]][ut] != sim$adj_base[ut])
  expect_gt(flips_platform, 0.04)
  expect_lt(flips_platform, 0.18)
  g <- sim$platform_data[[2]]$adj_list
  expect_identical(g[[1]], sim$platform_adj[[2]])
  expect_lt(mean(g[[2]][ut] != g[[1]][ut]), 0.12)
  # data dimensions follow p_vec
  expect_equal(ncol(sim$platform_data[[2]]$data_list[[1]]), 30)
})
