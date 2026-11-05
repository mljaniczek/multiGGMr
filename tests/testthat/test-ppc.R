sim <- simulate_multiggm(K = 2, p = 6, n = 150, seed = 1)
fit <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs", burnin = 300, nsave = 300,
                     seed = 1)

test_that("observed statistics match a direct computation", {
  pc <- posterior_predictive_check(fit, sim$data_list, ndraws = 50, seed = 1)
  X1 <- scale(sim$data_list[[1]], scale = FALSE)
  P1 <- -cov2cor(solve(crossprod(X1))); diag(P1) <- 1
  X2 <- scale(sim$data_list[[2]], scale = FALSE)
  P2 <- -cov2cor(solve(crossprod(X2))); diag(P2) <- 1
  ut <- upper.tri(P1)
  e <- pc$edges
  expect_equal(e$observed[e$type == "group" & e$group1 == 1], P1[ut])
  expect_equal(e$observed[e$type == "pair"], P1[ut] - P2[ut])
  g <- pc$global
  expect_equal(g$observed[g$statistic == "sum_sq_diff"], sum((P1[ut] - P2[ut])^2))
  R1 <- cov2cor(crossprod(X1))
  expect_equal(g$observed[g$statistic == "sum_sq_cor" & g$group1 == 1], sum(R1[ut]^2))
  # realized discrepancies average T(y, Omega_s) over the draws used
  pm <- posterior_predictive_check(fit, sim$data_list, ndraws = 5, use = "mean", seed = 1)
  U <- chol(fit$C_mean[, , 1])
  m <- rowSums((X1 %*% t(U))^2)
  gm <- pm$global
  expect_equal(gm$observed[gm$statistic == "kurtosis" & gm$group1 == 1], mean(m^2))
  expect_equal(gm$observed[gm$statistic == "max_mahal" & gm$group1 == 1], max(m))
  expect_true(all(c("ppp", "p_two", "z") %in% names(e)))
  expect_equal(nrow(e), 3 * 15)
})

test_that("a correctly specified model is not flagged; heavy tails are", {
  pc <- posterior_predictive_check(fit, sim$data_list, ndraws = 150, seed = 2)
  expect_true(all(pc$global$ppp > 0.01 & pc$global$ppp < 0.99))
  expect_lt(mean(pc$edges$p_two < 0.05), 0.1)
  # multivariate t with 3 degrees of freedom: same covariance structure, heavy tails
  set.seed(3)
  tdat <- lapply(sim$Omega_list, function(O) {
    Z <- matrix(rnorm(150 * 6), 150) %*% chol(solve(O))
    Z / sqrt(rchisq(150, 3) / 3)
  })
  ft <- multiggm_mcmc(data_list = tdat, method = "ssvs", burnin = 300, nsave = 300, seed = 1)
  pt <- posterior_predictive_check(ft, tdat, ndraws = 150, seed = 2)
  k <- pt$global[pt$global$statistic == "kurtosis", "ppp"]
  expect_true(all(k < 0.05))
})

test_that("options, other samplers and errors", {
  a <- posterior_predictive_check(fit, sim$data_list, ndraws = 20, seed = 5)
  b <- posterior_predictive_check(fit, sim$data_list, ndraws = 20, seed = 5)
  expect_identical(a, b)
  expect_equal(a$n_draws, 20L)
  expect_identical(a$edge_stat, "pcor")
  expect_identical(posterior_predictive_check(fit, sim$data_list, ndraws = 5,
                                              edge_stat = "cor")$edge_stat, "cor")
  expect_output(print(a), "two-sided ppp")
  expect_s3_class(plot_ppc(a), "ggplot")
  expect_s3_class(plot_ppc(a, "edges"), "ggplot")

  nd <- multiggm_mcmc(data_list = sim$data_list, method = "ssvs", burnin = 20, nsave = 20,
                      seed = 1, save_draws = FALSE)
  expect_error(posterior_predictive_check(nd, sim$data_list), "save_draws = TRUE")
  expect_identical(posterior_predictive_check(nd, sim$data_list, use = "mean",
                                              ndraws = 5)$use, "mean")
  expect_error(posterior_predictive_check(fit, sim$data_list[1]), "K = 2")
  small <- lapply(sim$data_list, function(x) x[1:6, ])
  expect_error(posterior_predictive_check(fit, small, edge_stat = "pcor"), "p \\+ 1")
  expect_identical(posterior_predictive_check(fit, small, ndraws = 5)$edge_stat, "cor")

  sim3 <- simulate_multiggm(K = 3, p = 5, n = 60, seed = 2)
  m <- multiggm_mcmc(data_list = sim3$data_list, method = "mghs", burnin = 50, nsave = 50,
                     nchains = 2, seed = 1)
  pm <- posterior_predictive_check(m, sim3$data_list, ndraws = 30, seed = 1)
  expect_equal(sum(pm$global$statistic == "sum_sq_diff"), 3L)
  expect_equal(nrow(pm$edges), (3 + 3) * 10)
})
