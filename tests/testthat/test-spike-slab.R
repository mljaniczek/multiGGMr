test_that("spike-and-slab log ratios match the explicit densities", {
  sh <- 2; sc <- 0.2; shp <- 1; scp <- 1 / 3; w <- 0.6
  lp <- function(x) if (x == 0) log(1 - w) else log(w) + dgamma(x, sh, scale = sc, log = TRUE)
  lq <- function(x) dgamma(x, shp, scale = scp, log = TRUE)
  # spike -> slab
  expect_equal(spike_slab_log_ratio_cpp(0, 0.7, sh, sc, shp, scp, w),
               lp(0.7) - lp(0) - lq(0.7))
  # slab -> spike
  expect_equal(spike_slab_log_ratio_cpp(0.7, 0, sh, sc, shp, scp, w),
               lp(0) - lp(0.7) + lq(0.7))
  # within slab
  expect_equal(spike_slab_log_ratio_cpp(0.7, 0.3, sh, sc, shp, scp, w),
               lp(0.3) - lp(0.7) + lq(0.7) - lq(0.3))
})

test_that("with no data the chain samples the Bernoulli-Gamma prior", {
  # alpha != alpha_prop and beta_prop != 1 exercise both bugs the old SSVS and
  # platform code had (scale/rate mix-up, sign of the log(theta) term).
  set.seed(1)
  shape <- 2; rate <- 5; shape_p <- 1; rate_p <- 3; w <- 0.6
  x <- spike_slab_prior_chain_cpp(2e5, shape, 1 / rate, shape_p, 1 / rate_p, w)
  slab <- x[x > 0]
  expect_equal(mean(x > 0), w, tolerance = 0.03)
  expect_equal(mean(slab), shape / rate, tolerance = 0.03)
  expect_equal(var(slab), shape / rate^2, tolerance = 0.06)
  expect_gt(suppressWarnings(ks.test(slab[seq(1, length(slab), by = 20)],
                                     "pgamma", shape, rate)$p.value), 0.001)
})

test_that("prior is recovered when w is small (catches the old sign error)", {
  # SSVS default shapes with w = 0.2: the old spike -> slab sign gave
  # P(slab) of about 0.35 instead of 0.20.
  set.seed(2)
  x <- spike_slab_prior_chain_cpp(1e5, 2, 1 / 5, 1, 1, 0.2)
  expect_equal(mean(x > 0), 0.2, tolerance = 0.1)
  expect_equal(mean(x[x > 0]), 0.4, tolerance = 0.05)
})
