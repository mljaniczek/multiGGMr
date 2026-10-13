# Validation against upstream mGHS (github.com/cbusatto/mGHS@1c9fd70).
# Skipped unless upstream is installed; set MULTIGGM_MGHS_LIB to the library
# that holds it if that is not on .libPaths().

mghs_lib <- Sys.getenv("MULTIGGM_MGHS_LIB")
if (nzchar(mghs_lib)) .libPaths(c(.libPaths(), mghs_lib))

run_both <- function(K, p, n, B, bn, seed) {
  set.seed(100 + seed)
  L <- chol(0.5^abs(outer(1:p, 1:p, "-")))
  S <- lapply(1:K, function(k) crossprod(scale(matrix(rnorm(n * p), n, p) %*% L, scale = FALSE)))
  nv <- rep(n, K)
  hyp <- multiGGM:::.mghs_hyper(NULL)$hyp_ta
  set.seed(seed)
  utils::capture.output(up <- mGHS::mGHS(B, bn, nv, S, p, hyp, chain = 1L, ping = B + 1L))
  set.seed(seed)
  me <- mghs_engine_cpp(bn, B - bn, 1L, nv, S, p, hyp, TRUE, 0.05, FALSE, 0L)
  list(up = up, me = me, K = K, p = p, ns = B - bn)
}

test_that("packaged sampler reproduces upstream draws exactly", {
  skip_on_cran()
  skip_if_not_installed("mGHS")
  cases <- list(c(K = 2, p = 8, n = 50, B = 300, bn = 100, seed = 2),
                c(K = 3, p = 10, n = 40, B = 300, bn = 100, seed = 3),
                c(K = 2, p = 15, n = 100, B = 200, bn = 50, seed = 4))
  for (cs in cases) {
    r <- do.call(run_both, as.list(cs))
    expect_identical(r$up$stuck, 0L)
    expect_identical(r$me$stuck, 0L)
    pairs <- do.call(rbind, lapply(seq_len(r$p - 1), function(i) cbind(i, (i + 1):r$p)))
    for (h in seq_len(r$K)) {
      cols <- (seq_len(r$ns) - 1) * r$K + h
      ours_ij <- vapply(seq_len(r$ns), function(s) r$me$C_save[, , h, s][pairs], numeric(nrow(pairs)))
      ours_jj <- vapply(seq_len(r$ns), function(s) diag(r$me$C_save[, , h, s]), numeric(r$p))
      expect_identical(ours_ij, unname(r$up$Omega_ij[, cols]))
      expect_identical(ours_jj, unname(r$up$Omega_jj[, cols]))
      expect_identical(r$me$tau_save[h, ], unname(r$up$tau[h, ]))
      ours_R <- vapply(seq_len(r$ns), function(s) r$me$R_save[h, , s], numeric(r$K))
      expect_equal(ours_R, unname(r$up$R[, cols]), tolerance = 1e-14)
    }
  }
})

test_that("upstream row-wise pair order maps onto upper.tri as documented", {
  skip_on_cran()
  skip_if_not_installed("mGHS")
  r <- run_both(K = 2, p = 6, n = 80, B = 150, bn = 50, seed = 7)
  m <- matrix(0, 6, 6)
  m[lower.tri(m)] <- rowMeans(r$up$Omega_ij[, seq(1, 2 * r$ns, by = 2)])
  m <- m + t(m)
  ours <- apply(r$me$C_save[, , 1, ], c(1, 2), mean)
  expect_equal(m[upper.tri(m)], ours[upper.tri(ours)], tolerance = 1e-12)
})
