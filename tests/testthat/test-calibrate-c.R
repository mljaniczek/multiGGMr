# Null data: two groups of different sizes drawn from one network
null_data <- function(p = 8, n = c(60, 40), seed = 1) {
  sim <- simulate_multiggm(K = 2, p = p, n = max(n), seed = seed)
  set.seed(seed + 1)
  L <- chol(solve(sim$Omega_list[[1]]))
  list(matrix(rnorm(n[1] * p), n[1]) %*% L + 3,     # a mean shift, removed by centring
       matrix(rnorm(n[2] * p), n[2]) %*% L)
}

test_that("permutations keep group sizes, centre within original groups, and are reproducible", {
  X <- null_data()
  P <- multiGGM:::.make_perms(60, 40, 5, seed = 3)
  expect_identical(P, multiGGM:::.make_perms(60, 40, 5, seed = 3))
  expect_false(identical(P, multiGGM:::.make_perms(60, 40, 5, seed = 4)))
  Xp <- multiGGM:::.perm_data(X, c(1, 2), P[[1]])
  expect_equal(vapply(Xp, nrow, 1), c(60, 40))
  # the pooled rows are the original rows, each centred on its own group means
  Xc <- lapply(X, function(x) scale(x, scale = FALSE))
  pool <- rbind(Xc[[1]], Xc[[2]])
  expect_equal(unname(Xp[[1]]), unname(pool[P[[1]]$idx_a, ]))
  expect_equal(unname(rbind(Xp[[1]], Xp[[2]])[order(c(P[[1]]$idx_a, setdiff(1:100, P[[1]]$idx_a))), ]),
               unname(pool))
  # the RNG state of the caller is left untouched
  set.seed(10); a <- runif(1); set.seed(10); invisible(multiGGM:::.make_perms(60, 40, 5, 3))
  expect_identical(runif(1), a)
})

test_that("the Weibull tail estimate recovers a known quantile", {
  set.seed(1)
  grid <- seq(0.01, 0.6, by = 0.01)
  T <- rweibull(3000, shape = 2, scale = 0.15)
  Tg <- vapply(T, function(t) if (t < grid[1]) 0 else max(grid[grid <= t]), 1)
  f <- multiGGM:::.fit_last_call(Tg, grid, alpha = 0.05)
  expect_equal(f$c, qweibull(0.95, 2, 0.15), tolerance = 0.03)
  expect_lt(f$lower, qweibull(0.95, 2, 0.15)); expect_gt(f$upper, qweibull(0.95, 2, 0.15))
})

test_that("calibrate_c() runs, is reproducible, and matches the selection rule", {
  skip_on_cran()
  X <- null_data()
  fit <- multiggm_mcmc(data_list = X, method = "ssvs", burnin = 150, nsave = 150, seed = 1)
  a <- calibrate_c(X, fit = fit, n_perm = 6, min_perm = 6, verbose = FALSE)
  b <- calibrate_c(X, fit = fit, n_perm = 6, min_perm = 6, verbose = FALSE)
  expect_s3_class(a, "multiggm_c_calibration")
  expect_identical(a$calls, b$calls)
  expect_identical(a$suggested_c, b$suggested_c)
  expect_equal(nrow(a$perms), 6L)
  expect_identical(a$settings$burnin, 150)
  expect_equal(dim(a$calls), c(6L, length(a$summary$c)))
  expect_named(as.data.frame(a), c("c", "mean_calls", "se_calls", "share_any", "share_any_lower",
                                   "share_any_upper", "p_any_fitted", "meets_target"))
  expect_output(print(a), "Suggested c")
  expect_s3_class(plot(a), "ggplot")
  expect_s3_class(plot(a, "mean"), "ggplot")
  if (.Platform$OS.type == "unix") {
    p <- calibrate_c(X, fit = fit, n_perm = 6, min_perm = 6, verbose = FALSE,
                     parallel = TRUE, ncores = 2)
    expect_identical(p$calls, a$calls)
  }
  # rebuilding permutation 2 and refitting reproduces its call counts
  perm <- multiGGM:::.make_perms(60, 40, 6, a$seed)[[2]]
  expect_identical(perm$seed, as.integer(a$perms$seed[2]))
  g <- c(0.02, 0.05, 0.1)
  f2 <- multiggm_mcmc(data_list = multiGGM:::.perm_data(X, c(1, 2), perm), method = "ssvs",
                      burnin = 150, nsave = 150, seed = perm$seed, save_draws = FALSE,
                      diff_deltas = a$summary$c, mcse_batches = 0)
  for (cc in g)
    expect_equal(select_diff_edges(f2, c = cc, scale = "pcor")$n_selected,
                 unname(a$calls[2, match(cc, a$summary$c)]), label = paste("c =", cc))
  # mean_false_calls: first grid value from which the mean stays <= 1
  m <- calibrate_c(X, fit = fit, n_perm = 6, target = "mean_false_calls", verbose = FALSE)
  s <- m$summary
  expect_true(all(s$mean_calls[s$c >= m$suggested_c] <= 1))
  expect_true(m$suggested_c == s$c[1] || s$mean_calls[match(m$suggested_c, s$c) - 1] > 1)
})

test_that("settings come from the fit, with fallbacks for fits without fit$settings", {
  skip_on_cran()
  X <- null_data()
  fl <- multiggm_mcmc(data_list = X, method = "ssvs", burnin = 120, nsave = 50, thin = 2,
                      nchains = 2, seed = 1, hyper = list(v0 = 0.03^2))
  st <- multiGGM:::.settings_from_fit(fl)
  expect_identical(st$burnin, 120); expect_identical(st$nchains, 2)
  # an older fit: no settings recorded; burnin from the call, the rest from the chains
  old <- fl; old$settings <- NULL
  so <- multiGGM:::.settings_from_fit(old)
  expect_equal(so$burnin, 120); expect_equal(so$nsave, 50); expect_equal(so$thin, 2)
  expect_equal(so$nchains, 2L); expect_equal(so$hyper$v0, 0.03^2)
  expect_null(so$hyper$a_prop)
  # a single-chain fit without settings does not know its burnin
  one <- multiggm_mcmc(data_list = X, method = "ssvs", burnin = 30, nsave = 30, seed = 1)
  one$settings <- NULL
  expect_error(calibrate_c(X, fit = one, n_perm = 2, verbose = FALSE), "burnin")
  expect_error(calibrate_c(X, method = "ssvs", n_perm = 2, verbose = FALSE), "burnin")
  r <- suppressWarnings(calibrate_c(X, fit = one, mcmc_args = list(burnin = 30), n_perm = 2,
                                    verbose = FALSE))
  expect_identical(r$settings$burnin, 30)
  expect_warning(calibrate_c(X, fit = one, mcmc_args = list(burnin = 30, seed = 4),
                             n_perm = 2, verbose = FALSE), "Ignoring") |>
    suppressWarnings()
})

test_that("failed permutation fits are recorded and skipped", {
  skip_on_cran()
  X <- null_data()
  real <- multiGGM::multiggm_mcmc
  n_call <- 0
  local_mocked_bindings(multiggm_mcmc = function(...) {
    n_call <<- n_call + 1
    if (n_call %% 3 == 0) stop("chain stuck (test)")
    real(...)
  })
  w <- character(0)
  r <- withCallingHandlers(
    calibrate_c(X, method = "ssvs", mcmc_args = list(burnin = 50, nsave = 50), n_perm = 6,
                verbose = FALSE),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("Only 4 of 6 permutation fits succeeded", w)))
  expect_equal(sum(!r$perms$ok), 2L)
  expect_true(all(grepl("stuck", r$perms$error[!r$perms$ok])))
  expect_equal(r$n_ok, 4L)
  expect_true(all(is.na(r$calls[!r$perms$ok, ])))
})

test_that("the suggested c meets its target on fresh permutations", {
  skip_on_cran()
  X <- null_data(p = 10, n = c(100, 70), seed = 5)
  args <- list(burnin = 300, nsave = 300)
  cal <- calibrate_c(X, method = "ssvs", mcmc_args = args, n_perm = 30, c_tol = 0,
                     seed = 1, verbose = FALSE)
  expect_true(is.finite(cal$suggested_c))
  fresh <- calibrate_c(X, method = "ssvs", mcmc_args = args, n_perm = 40, c_tol = 0,
                       seed = 999, c_grid = cal$summary$c, verbose = FALSE)
  share <- fresh$summary$share_any[match(cal$suggested_c, fresh$summary$c)]
  # target 0.05; 40 fresh permutations, so allow for Monte Carlo error
  expect_lte(share, 0.15)
  # and the calibration is not wildly conservative: a much smaller c fails it
  expect_gt(fresh$summary$share_any[1], 0.3)
})

test_that("select_diff_edges() accepts a calibration, and input checks", {
  skip_on_cran()
  X <- null_data()
  fit <- multiggm_mcmc(data_list = X, method = "ssvs", burnin = 50, nsave = 50, seed = 1,
                       diff_deltas = seq(0.01, 0.8, by = 0.01))
  cal <- calibrate_c(X, fit = fit, n_perm = 5, verbose = FALSE)
  if (is.finite(cal$suggested_c)) {
    s <- select_diff_edges(fit, c = cal)
    expect_identical(s$delta, cal$suggested_c)
    expect_error(select_diff_edges(fit, c = cal, scale = "precision"), "pcor")
  }
  expect_error(calibrate_c(X, method = "ssvs", mcmc_args = list(burnin = 5, nsave = 5),
                           c_grid = c(0.1, 1.5), verbose = FALSE), "c_grid")
  expect_error(calibrate_c(X[1], method = "ssvs"), "at least two")
  w <- character(0)
  withCallingHandlers(
    try(calibrate_c(X, method = "mghs", mcmc_args = list(burnin = 5, nsave = 5), n_perm = 2,
                    verbose = FALSE), silent = TRUE),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("R_eta", w)))
})
