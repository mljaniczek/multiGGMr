test_that("bayes_fdr selects the largest set with mean(1 - prob) <= alpha", {
  sel <- bayes_fdr(c(0.6, 0.99, 0.9, 0.97), alpha = 0.05)
  expect_equal(as.logical(sel), c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(attr(sel, "n_selected"), 3L)
  expect_equal(attr(sel, "expected_fdr"), mean(c(0.01, 0.03, 0.1)))
  expect_equal(attr(sel, "cutoff"), 0.9)

  none <- bayes_fdr(c(0.5, 0.4), alpha = 0.05)
  expect_false(any(none))
  expect_identical(attr(none, "expected_fdr"), 0)
  expect_error(bayes_fdr(c(0.5, 1.2)), "prob must lie")
})

test_that("bayes_fdr on a matrix uses the upper triangle and returns a symmetric matrix", {
  P <- matrix(0, 4, 4)
  P[upper.tri(P)] <- c(0.99, 0.2, 0.98, 0.1, 0.3, 0.97)
  P <- P + t(P)
  sel <- bayes_fdr(P, alpha = 0.05)
  expect_true(isSymmetric(unclass(sel) * 1))
  expect_false(any(diag(sel)))
  expect_equal(sum(sel[upper.tri(sel)]), 3)
})

sim <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 11)
fit2 <- multiggm_mcmc(data_list = sim$data_list, burnin = 100, nsave = 40,
                      nchains = 2, seed = 3)

test_that("diff_draws pools chains and matches C_save", {
  D <- diff_draws(fit2, scale = "precision")
  expect_equal(dim(D), c(6, 6, 80))
  ch2 <- fit2$chains[[2]]
  expect_equal(D[1, 2, 41], ch2$C_save[1, 2, 1, 1] - ch2$C_save[1, 2, 2, 1])
  expect_true(all(D[cbind(1:6, 1:6, 1)] == 0))
  expect_equal(dim(diff_draws(fit2, chain = 1)), c(6, 6, 40))

  Dp <- diff_draws(fit2, groups = c(2, 1), scale = "pcor", chain = 1)
  pc <- posterior_pcor(fit2, chain = 1)
  expect_equal(Dp[2, 5, 7], pc[2, 5, 2, 7] - pc[2, 5, 1, 7], tolerance = 1e-12)
  expect_error(diff_draws(fit2, groups = c(1, 1)), "distinct")
})

test_that("diff_summary agrees with the pooled draws", {
  ds <- diff_summary(fit2, deltas = c(0.05, 0.1))
  D <- diff_draws(fit2)
  row <- which(ds$i == 2 & ds$j == 4)
  d <- D[2, 4, ]
  expect_equal(nrow(ds), choose(6, 2))
  expect_equal(ds$mean[row], mean(d))
  expect_equal(ds$sd[row], sd(d))
  expect_equal(ds$p_exceed_0.1[row], mean(abs(d) > 0.1))
  expect_equal(ds$p_sign[row], max(mean(d > 0), mean(d < 0)))
  expect_equal(ds$lower[row], unname(quantile(d, 0.025)))
  expect_identical(attr(ds, "n_draws"), 80L)
  expect_false("lower" %in% names(diff_summary(fit2, probs = NULL)))
})

test_that("select_diff_edges supports magnitude and structural rules", {
  s1 <- select_diff_edges(fit2, alpha = 0.2, delta = 0.05)
  expect_s3_class(s1, "multiggm_diff_selection")
  ds <- diff_summary(fit2, deltas = 0.05, probs = NULL)
  expect_equal(s1$prob[upper.tri(s1$prob)], ds$p_exceed_0.05)
  expect_equal(s1$n_selected, nrow(s1$edges))
  expect_lte(s1$expected_fdr, 0.2)

  s2 <- select_diff_edges(fit2, rule = "structural", alpha = 0.5)
  pd <- (pip_diff_edge(fit2, chain = 1) + pip_diff_edge(fit2, chain = 2)) / 2
  expect_equal(s2$prob, pd, ignore_attr = TRUE)
  expect_output(print(s1), "multiggm_diff_selection")
})
