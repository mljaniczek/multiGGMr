sim_r <- simulate_multiggm(K = 2, p = 4, n = 80, seed = 6)
fit_twice <- function(...) list(multiggm_mcmc(data_list = sim_r$data_list, ...),
                                multiggm_mcmc(data_list = sim_r$data_list, ...))

test_that("the same seed gives identical output for each method", {
  for (m in c("gwishart", "ssvs", "mghs")) {
    f <- fit_twice(method = m, burnin = 30, nsave = 20, seed = 11)
    expect_identical(f[[1]]$C_save, f[[2]]$C_save, label = m)
  }
  g <- fit_twice(method = "gwishart", burnin = 30, nsave = 20, seed = 11)
  expect_identical(g[[1]]$adj_save, g[[2]]$adj_save)
  expect_identical(g[[1]]$Theta_save, g[[2]]$Theta_save)
  expect_identical(g[[1]]$nu_save, g[[2]]$nu_save)

  sp <- simulate_multiggm_platform(K = 2, S = 2, p_vec = c(4, 3), n = 40, seed = 2)
  run_p <- function() multiggm_mcmc(method = "ssvs_platform", platform_data = sp$platform_data,
                                    burnin = 20, nsave = 15, seed = 11)
  a <- run_p(); b <- run_p()
  for (s in 1:2) expect_identical(a$platforms[[s]]$C_save, b$platforms[[s]]$C_save)
  expect_identical(a$Phi_save, b$Phi_save)

  d <- multiggm_mcmc(data_list = sim_r$data_list, burnin = 30, nsave = 20, seed = 12)
  expect_false(identical(d$C_save, g[[1]]$C_save))
})

test_that("parallel and serial runs of two chains are identical", {
  skip_on_os("windows")
  args <- list(data_list = sim_r$data_list, burnin = 30, nsave = 20, nchains = 2, seed = 5)
  ser <- do.call(multiggm_mcmc, c(args, parallel = FALSE))
  par <- do.call(multiggm_mcmc, c(args, parallel = TRUE, ncores = 2))
  for (ch in 1:2) {
    expect_identical(par$chains[[ch]]$C_save, ser$chains[[ch]]$C_save)
    expect_identical(par$chains[[ch]]$adj_save, ser$chains[[ch]]$adj_save)
  }
  expect_false(identical(ser$chains[[1]]$C_save, ser$chains[[2]]$C_save))
})

test_that("method = \"ssvs\" works directly", {
  fit <- multiggm_mcmc(data_list = sim_r$data_list, method = "ssvs", burnin = 30,
                       nsave = 25, seed = 1)
  expect_s3_class(fit, "multiggm_fit")
  expect_identical(fit$method, "ssvs")
  expect_equal(dim(fit$C_save), c(4, 4, 2, 25))
  expect_equal(dim(fit$adj_save), c(4, 4, 2, 25))
  expect_equal(fit$hyper$v0, 0.02^2)
  pip <- pip_edges(fit)
  expect_true(all(pip >= 0 & pip <= 1))
  expect_true(all(is.finite(unlist(coef(fit)))))
})

test_that("R and C++ G-Wishart engines agree", {
  # Both engines call the same compiled helpers and draw random numbers in
  # the same order, so with one seed they give the same chain.
  args <- list(data_list = sim_r$data_list, burnin = 50, nsave = 100, seed = 7)
  fc <- do.call(multiggm_mcmc, c(args, engine = "cpp"))
  fr <- do.call(multiggm_mcmc, c(args, engine = "R"))
  expect_equal(pip_edges(fr), pip_edges(fc), tolerance = 1e-6)
  expect_equal(coef(fr), coef(fc), tolerance = 1e-6)
  expect_equal(fr$Theta_save, fc$Theta_save, tolerance = 1e-6)
})
