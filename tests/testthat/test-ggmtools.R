test_that("shared functions are the ggmtools objects", {
  expect_identical(roc_auc, ggmtools::roc_auc)
  expect_identical(precision_to_pcor, ggmtools::precision_to_pcor)
  expect_identical(confusion_at_threshold, ggmtools::confusion_at_threshold)
  expect_identical(plot_roc, ggmtools::plot_roc)
})

test_that("pip_table() and pcor_estimate() methods for multiggm fits", {
  sim <- simulate_multiggm(K = 2, p = 5, n = 60, seed = 1)
  fit <- multiggm_mcmc(data_list = sim$data_list, burnin = 50, nsave = 30)
  pt <- pip_table(fit)
  expect_named(pt, c("i", "j", "node1", "node2", "group", "pip"))
  expect_equal(nrow(pt), 2 * 10)
  pip <- pip_edges(fit)
  expect_equal(pt$pip[pt$group == "Group_2"], pip[, , 2][upper.tri(pip[, , 2])])
  pc <- pcor_estimate(fit, varnames = letters[1:5])
  expect_named(pc, c("Group_1", "Group_2"))
  expect_equal(unname(pc$Group_1), unname(fitted(fit)$Group_1))
  expect_equal(rownames(pc$Group_1), letters[1:5])
  expect_equal(pip_table(fit, varnames = letters[1:5])$node2[1:3], c("b", "c", "c"))
})

test_that("roc_auc from ggmtools handles tied PIPs", {
  truth <- matrix(0, 6, 6); truth[1, 2] <- truth[3, 4] <- 1; truth <- truth + t(truth)
  score <- matrix(0, 6, 6)                       # every edge tied at 0
  expect_equal(roc_auc(score, truth)$auc, 0.5)
})
