
test_that("ssvs_platform runs with 2 platforms and 2 groups", {
  sim <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(8, 6),
                                     n = 60, seed = 123)
  expect_equal(sim$K, 2)
  expect_equal(sim$S, 2)
  expect_equal(sim$p_vec, c(8L, 6L))

  fit <- multiggm_mcmc(
    method = "ssvs_platform",
    platform_data = sim$platform_data,
    burnin = 200, nsave = 100
  )

  expect_s3_class(fit, "multiggm_platform_fit")
  expect_equal(fit$K, 2)
  expect_equal(fit$S, 2)
  expect_equal(fit$p_vec, c(8L, 6L))
  expect_equal(fit$method, "ssvs_platform")

  # Check platform-level saves

  expect_equal(dim(fit$Phi_save), c(2, 2, 100))
  expect_equal(dim(fit$w_save), c(2, 2, 100))

  # Check per-platform saves
  expect_equal(length(fit$platforms), 2)

  # Platform 1: p = 8
  p1 <- fit$platforms[[1]]
  expect_equal(p1$p, 8)
  expect_equal(dim(p1$C_save), c(8, 8, 2, 100))
  expect_equal(dim(p1$Sig_save), c(8, 8, 2, 100))
  expect_equal(dim(p1$adj_save), c(8, 8, 2, 100))
  expect_equal(dim(p1$Theta_save), c(2, 2, 100))
  expect_equal(dim(p1$nu_save), c(8, 8, 100))

  # Platform 2: p = 6
  p2 <- fit$platforms[[2]]
  expect_equal(p2$p, 6)
  expect_equal(dim(p2$C_save), c(6, 6, 2, 100))
})

test_that("pip_edges works with platform fit", {
  sim <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(6, 5),
                                     n = 50, seed = 456)
  fit <- multiggm_mcmc(
    method = "ssvs_platform",
    platform_data = sim$platform_data,
    burnin = 100, nsave = 50
  )

  pip1 <- pip_edges(fit, platform = 1)
  expect_equal(dim(pip1), c(6, 6, 2))
  expect_true(all(pip1 >= 0 & pip1 <= 1))

  pip2 <- pip_edges(fit, platform = 2)
  expect_equal(dim(pip2), c(5, 5, 2))

  # Should error without platform arg
  expect_error(pip_edges(fit), "platform")
})

test_that("S3 methods work for platform fit", {
  sim <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(6, 5),
                                     n = 50, seed = 789)
  fit <- multiggm_mcmc(
    method = "ssvs_platform",
    platform_data = sim$platform_data,
    burnin = 100, nsave = 50
  )

  # print
  expect_output(print(fit), "multiggm_platform_fit")

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.multiggm_platform_fit")
  expect_output(print(s), "Platform-level")

  # coef
  cc <- coef(fit, platform = 1)
  expect_equal(length(cc), 2)
  expect_equal(dim(cc$Group_1), c(6, 6))

  cc_all <- coef(fit)
  expect_equal(length(cc_all), 2)
  expect_equal(dim(cc_all$Platform_1$Group_1), c(6, 6))
  expect_equal(dim(cc_all$Platform_2$Group_1), c(5, 5))
})
