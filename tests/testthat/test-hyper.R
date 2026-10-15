sim_h <- simulate_multiggm(K = 2, p = 5, n = 40, seed = 4)

test_that("a partial hyper list keeps the other SSVS defaults", {
  fit <- multiggm_mcmc(data_list = sim_h$data_list, method = "ssvs", burnin = 20,
                       nsave = 10, seed = 1, hyper = list(v0 = 0.01))
  expected <- modifyList(multiGGM:::.default_hyper("ssvs"), list(v0 = 0.01))
  expect_equal(fit$hyper[names(expected)], expected)
  expect_equal(fit$hyper$v0, 0.01)
})

test_that("a partial hyper list keeps the other G-Wishart defaults", {
  fit <- multiggm_mcmc(data_list = sim_h$data_list, burnin = 20, nsave = 10,
                       seed = 1, hyper = list(w = 0.5))
  expect_equal(fit$hyper$w, 0.5)
  expect_equal(fit$hyper$b_prior, 3)
  expect_equal(fit$hyper$D_prior, diag(5))
  expect_equal(fit$hyper$beta_prop, 5)
})

test_that("unknown or malformed hyperparameters are an error", {
  expect_error(multiggm_mcmc(data_list = sim_h$data_list, method = "ssvs", burnin = 5,
                             nsave = 5, hyper = list(vo = 0.01)),
               "Unknown hyperparameter for method = \"ssvs\": vo")
  expect_error(multiggm_mcmc(data_list = sim_h$data_list, burnin = 5, nsave = 5,
                             hyper = list(v0 = 0.01)), "Unknown hyperparameter.*gwishart")
  expect_error(multiggm_mcmc(data_list = sim_h$data_list, method = "mghs", burnin = 5,
                             nsave = 5, hyper = list(t_mu = 0.5)), "Unknown hyperparameter")
  expect_error(multiggm_mcmc(data_list = sim_h$data_list, burnin = 5, nsave = 5,
                             hyper = list(0.5)), "named list")
})
