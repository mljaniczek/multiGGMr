sim_c <- simulate_multiggm(K = 2, p = 5, n = 50, seed = 8)
fit_c <- multiggm_mcmc(data_list = sim_c$data_list, burnin = 30, nsave = 20,
                       nchains = 2, seed = 2)

test_that("posterior summaries pool chains by default and chain = selects one", {
  p1 <- pip_edges(fit_c, chain = 1); p2 <- pip_edges(fit_c, chain = 2)
  expect_equal(pip_edges(fit_c), (p1 + p2) / 2)
  expect_equal(p1, pip_edges(fit_c$chains[[1]]))
  expect_equal(pip_diff_edge(fit_c),
               (pip_diff_edge(fit_c, chain = 1) + pip_diff_edge(fit_c, chain = 2)) / 2)

  pc <- posterior_pcor(fit_c)
  expect_equal(dim(pc), c(5, 5, 2, 40))
  expect_equal(pc[, , , 21:40], posterior_pcor(fit_c, chain = 2))

  pooled <- posterior_precision(fit_c)
  manual <- (apply(fit_c$chains[[1]]$C_save[, , 1, ], c(1, 2), sum) +
             apply(fit_c$chains[[2]]$C_save[, , 1, ], c(1, 2), sum)) / 40
  expect_equal(pooled$Group_1, manual)
  expect_equal(posterior_precision(fit_c, chain = 1), posterior_precision(fit_c$chains[[1]]))
})

test_that("multiggm_fit_list has summary, coef, fitted and plot methods", {
  expect_equal(coef(fit_c), posterior_precision(fit_c))
  f <- fitted(fit_c)
  expect_equal(f$Group_2, apply(posterior_pcor(fit_c)[, , 2, ], c(1, 2), mean))
  s <- summary(fit_c)
  expect_s3_class(s, "summary.multiggm_fit_list")
  expect_length(s$chains, 2)
  expect_length(s$max_chain_diff, 2)
  expect_output(print(s), "Pooled selected edges")
  for (ty in c("trace_theta", "trace_edges", "pip"))
    expect_s3_class(plot(fit_c, type = ty), "ggplot")
})

test_that("summary-only mGHS chains pool their running means", {
  fit_m <- multiggm_mcmc(data_list = sim_c$data_list, method = "mghs", burnin = 40,
                         nsave = 30, nchains = 2, seed = 3, save_draws = FALSE)
  expect_equal(coef(fit_m)$Group_1,
               (fit_m$chains[[1]]$C_mean[, , 1] + fit_m$chains[[2]]$C_mean[, , 1]) / 2)
  expect_output(print(summary(fit_m)), "posterior mean partial correlation")
})

test_that("ssvs_platform refuses several chains and warns about parallel", {
  sp <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(4, 3), n = 30, seed = 1)
  expect_error(multiggm_mcmc(method = "ssvs_platform", platform_data = sp$platform_data,
                             burnin = 5, nsave = 5, nchains = 2), "single chain")
  expect_warning(multiggm_mcmc(method = "ssvs_platform", platform_data = sp$platform_data,
                               burnin = 5, nsave = 5, parallel = TRUE, seed = 1),
                 "no effect")
})
