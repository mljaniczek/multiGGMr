sym <- function(v, p) { m <- matrix(0, p, p); m[upper.tri(m)] <- v; m + t(m) }

test_that("roc_auc is the Mann-Whitney AUC with ties counted as 1/2", {
  p <- 6; n <- choose(p, 2)
  truth <- sym(rep(c(1, 0, 0), length.out = n), p)
  score <- sym(c(1, 1, 0, 0.5, 1, 0, 0, 0.5, 1, 0, 1, 1, 0, 0.5, 0)[seq_len(n)], p)
  y <- truth[upper.tri(truth)] > 0; s <- score[upper.tri(score)]
  pairs <- outer(s[y], s[!y], "-")
  expected <- mean((pairs > 0) + 0.5 * (pairs == 0))
  r <- roc_auc(score, truth)
  expect_equal(r$auc, expected)
  # the trapezoidal area under the returned curve agrees
  trap <- sum(diff(r$FPR) * (head(r$TPR, -1) + tail(r$TPR, -1)) / 2)
  expect_equal(trap, expected)
  expect_equal(c(head(r$FPR, 1), tail(r$FPR, 1), tail(r$TPR, 1)), c(0, 1, 1))
})

test_that("roc_auc does not depend on the order of tied scores", {
  p <- 8; n <- choose(p, 2)
  set.seed(1)
  truth_v <- rbinom(n, 1, 0.3)
  score_v <- sample(c(0, 1), n, replace = TRUE)   # heavily tied, like PIPs at 0 / 1
  r1 <- roc_auc(sym(score_v, p), sym(truth_v, p))
  perm <- sample(n)
  r2 <- roc_auc(sym(score_v[perm], p), sym(truth_v[perm], p))
  expect_equal(r1$auc, r2$auc)
  expect_equal(r1$FPR, r2$FPR)
  expect_equal(r1$TPR, r2$TPR)
  # all scores tied: chance level
  expect_equal(roc_auc(sym(rep(0.5, n), p), sym(truth_v, p))$auc, 0.5)
})
