
test_that("calc_mrf_logC matches direct enumeration", {
  set.seed(1)
  K <- 4
  Theta <- matrix(rnorm(K*K, sd = 0.1), K, K)
  Theta <- (Theta + t(Theta))/2
  diag(Theta) <- 0
  nu <- c(-2, -0.5, 0, 1)

  # direct R enumeration
  B <- as.matrix(expand.grid(rep(list(0:1), K)))
  s <- rowSums(B)
  q <- rowSums((B %*% Theta) * B)

  r_logC <- vapply(nu, function(nui) {
    x <- nui * s + q
    m <- max(x)
    m + log(sum(exp(x - m)))
  }, numeric(1))

  cpp_logC <- calc_mrf_logC(Theta, nu)

  expect_equal(cpp_logC, r_logC, tolerance = 1e-10)
})
