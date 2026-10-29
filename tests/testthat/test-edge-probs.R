# --- helpers ----------------------------------------------------------------

# Exact comparison in a normal (optimized) build; debug builds such as
# devtools::load_all() (-O0) round differently, so allow rounding there.
# pkgload marks namespaces loaded with load_all() with .__DEVTOOLS__.
dev_build <- function() exists(".__DEVTOOLS__", envir = asNamespace("multiGGM"), inherits = FALSE)
expect_same <- function(object, expected, platform, ...) {
  if (!dev_build() && identical(platform, R.version$platform)) {
    expect_identical(object, expected, ...)
  } else {
    expect_equal(object, expected, tolerance = 1e-8, ...)
  }
}

sim2 <- simulate_multiggm(K = 2, p = 6, n = 60, seed = 1)
sim3 <- simulate_multiggm(K = 3, p = 5, n = 50, seed = 2)

# --- 1. Reproducibility -------------------------------------------------------

test_that("mGHS outputs are unchanged from 0.4.0 for the same seed", {
  ref <- readRDS(test_path("fixtures", "ref_mghs_0.4.0.rds"))
  expect_identical(ref$version, "0.4.0")
  f2 <- multiggm_mcmc(data_list = sim2$data_list, method = "mghs", burnin = 60,
                      nsave = 40, seed = 1, diff_deltas = c(0.02, 0.1))
  f3 <- multiggm_mcmc(data_list = sim3$data_list, method = "mghs", burnin = 50,
                      nsave = 30, thin = 2, seed = 3, diff_deltas = c(0.05))
  for (nm in names(ref$mghs_k2)) expect_same(f2[[nm]], ref$mghs_k2[[nm]], ref$platform, label = nm)
  for (nm in names(ref$mghs_k3)) expect_same(f3[[nm]], ref$mghs_k3[[nm]], ref$platform, label = nm)
})

test_that("G-Wishart and SSVS draws are unchanged by the accumulators", {
  ref <- readRDS(test_path("fixtures", "ref_mrf_main.rds"))
  fits <- list(
    ssvs = multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 40, nsave = 30, seed = 1),
    gwishart = multiggm_mcmc(data_list = sim2$data_list, burnin = 40, nsave = 30, seed = 1),
    ssvs_k3 = multiggm_mcmc(data_list = sim3$data_list, method = "ssvs", burnin = 30, nsave = 20,
                            thin = 2, seed = 4))
  for (f in names(fits)) for (nm in names(ref[[f]]))
    expect_same(fits[[f]][[nm]], ref[[f]][[nm]], ref$platform, label = paste(f, nm))
})

test_that("save_draws = FALSE gives the same running summaries and Theta trace", {
  for (m in c("ssvs", "gwishart", "mghs")) {
    a <- multiggm_mcmc(data_list = sim2$data_list, method = m, burnin = 30, nsave = 20, seed = 9)
    b <- multiggm_mcmc(data_list = sim2$data_list, method = m, burnin = 30, nsave = 20, seed = 9,
                       save_draws = FALSE)
    expect_identical(b$edge_acc, a$edge_acc, label = m)
    expect_null(b$C_save)
    if (m != "mghs") {
      expect_identical(b$Theta_save, a$Theta_save)
      expect_null(b$adj_save); expect_null(b$nu_save)
      expect_equal(pip_edges(b), apply(a$adj_save, 1:3, mean), tolerance = 1e-12)
    }
  }
})

# --- 2. Running sums equal the draw-based computation ------------------------

test_that("running sums equal the draw-based values (G-Wishart, SSVS)", {
  for (m in c("ssvs", "gwishart")) for (sim in list(sim2, sim3)) {
    f <- multiggm_mcmc(data_list = sim$data_list, method = m, burnin = 30, nsave = 40,
                       seed = 5, diff_deltas = c(0.02, 0.1, 0.25))
    ref <- multiGGM:::.edge_acc_from_draws(f$C_save, f$adj_save, f$edge_acc$deltas)
    expect_equal(f$edge_acc, ref, tolerance = 1e-12, label = m)
    expect_equal(f$C_mean, apply(f$C_save, 1:3, mean), tolerance = 1e-12)
    expect_equal(f$pcor_mean, apply(posterior_pcor(f), 1:3, mean), tolerance = 1e-12)
  }
})

test_that("running sums equal the draw-based values (mGHS, incl. weights and cut model)", {
  for (sim in list(sim2, sim3)) {
    f <- multiggm_mcmc(data_list = sim$data_list, method = "mghs", burnin = 40, nsave = 40,
                       seed = 6, diff_deltas = c(0.02, 0.1), save_weights = TRUE)
    ref <- multiGGM:::.edge_acc_from_draws(f$C_save, NULL, f$edge_acc$deltas,
                                           f$weight_save, f$cut_save)
    expect_equal(f$edge_acc, ref, tolerance = 1e-12)
    # cut-model indicator mapping agrees with the legacy running mean
    expect_equal(apply(f$cut_save, 1:3, mean)[upper.tri(diag(f$p))],
                 f$cut_z_mean[, , ][upper.tri(diag(f$p))], tolerance = 1e-12)
    # pair magnitude sums agree with the legacy diff_acc
    expect_equal(f$edge_acc$pair$precision$sum, f$diff_acc$precision$sum, tolerance = 1e-12)
    expect_equal(f$edge_acc$pair$pcor$n_exceed, f$diff_acc$pcor$n_exceed)
  }
})

test_that("edge_probs() matches hand computation from draws", {
  f <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 30, nsave = 50, seed = 2)
  ep <- edge_probs(f, scale = "pcor")
  pc <- posterior_pcor(f)
  row <- ep$type == "pair" & ep$i == 2 & ep$j == 4
  D <- pc[2, 4, 1, ] - pc[2, 4, 2, ]
  G <- f$adj_save[2, 4, 1, ] != f$adj_save[2, 4, 2, ]
  expect_equal(ep$mean[row], mean(D))
  expect_equal(ep$sd[row], sd(D))
  expect_equal(ep$p_struct[row], mean(G))
  expect_equal(ep$p_exceed_0.05[row], mean(abs(D) > 0.05))
  expect_equal(ep$p_union_0.05[row], mean(G | abs(D) > 0.05))
  expect_equal(ep$p_sign[row], max(mean(D > 0), mean(D < 0)))
  grow <- ep$type == "group" & ep$group1 == 2 & ep$i == 2 & ep$j == 4
  expect_equal(ep$p_struct[grow], mean(f$adj_save[2, 4, 2, ]))
  expect_equal(ep$p_exceed_0.1[grow], mean(abs(pc[2, 4, 2, ]) > 0.1))
  # thresholds outside the grid come from the saved draws
  ep2 <- edge_probs(f, scale = "pcor", c = 0.07)
  expect_equal(ep2$p_union_0.07[row], mean(G | abs(D) > 0.07))
  # ... and are an error without draws
  g <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 10, nsave = 10,
                     seed = 2, save_draws = FALSE)
  expect_error(edge_probs(g, c = 0.07), "not in the grid")
})

# --- 3. Structural vs magnitude vs union ------------------------------------

diff_sim <- function(n = 600, seed = 11) {
  p <- 20
  path <- function(nodes) cbind(nodes[-length(nodes)], nodes[-1])
  A <- path(1:5); B <- path(11:15); C <- path(6:10)   # C: shared, unchanged
  O1 <- diag(p); O2 <- diag(p)
  for (e in seq_len(nrow(A))) O1[A[e, 1], A[e, 2]] <- 0.3
  for (e in seq_len(nrow(B))) { O1[B[e, 1], B[e, 2]] <- 0.3; O2[B[e, 1], B[e, 2]] <- 0.45 }
  for (e in seq_len(nrow(C))) { O1[C[e, 1], C[e, 2]] <- 0.3; O2[C[e, 1], C[e, 2]] <- 0.3 }
  O1 <- O1 + t(O1) - diag(p); O2 <- O2 + t(O2) - diag(p)
  set.seed(seed)
  X <- lapply(list(O1, O2), function(O) matrix(rnorm(n * p), n, p) %*% chol(solve(O)))
  list(X = X, A = A, B = B, C = C)
}

flagged <- function(sel, edges) apply(edges, 1, function(e) sel$selected[e[1], e[2]])

test_that("structural rules see deletions, magnitude sees scaling, union sees both", {
  skip_on_cran()
  # True |D| = 0.30 (block A, deleted in group 2) and 0.15 (block B, scaled by
  # 1.5); block C is identical in both groups. n is large enough that sampling
  # noise in D (~0.02) is well below c = 0.1.
  d <- diff_sim(n = 4000)
  f <- multiggm_mcmc(data_list = d$X, method = "ssvs", burnin = 300, nsave = 600, seed = 1,
                     save_draws = FALSE, diff_deltas = c(0.05, 0.1))
  st <- select_diff_edges(f, rule = "structural", decision = "median")
  mg <- select_diff_edges(f, rule = "magnitude", c = 0.1, scale = "precision", decision = "median")
  un <- select_diff_edges(f, rule = "union", c = 0.1, scale = "precision", decision = "median")
  expect_true(all(flagged(st, d$A)))
  expect_false(any(flagged(st, d$B)))
  expect_true(all(flagged(mg, d$B)))
  expect_true(all(flagged(un, d$A)) && all(flagged(un, d$B)))
  expect_false(any(flagged(un, d$C)))
  expect_lte(un$n_selected, nrow(d$A) + nrow(d$B) + 1)

  # mGHS: the magnitude rule separates the blocks. Its structural indicator
  # S = 1(w > 0.5) does not flag the scaled block, but it is also weak for
  # deletions: w = lambda^2 / (1 + lambda^2) excludes tau, so for a null edge
  # w stays near its prior and P(S = 1) is about 0.5 (see ?edge_probs).
  g <- multiggm_mcmc(data_list = d$X, method = "mghs", burnin = 300, nsave = 600, seed = 1,
                     save_draws = FALSE, diff_deltas = c(0.05, 0.1))
  mg <- select_diff_edges(g, rule = "magnitude", c = 0.1, scale = "precision", decision = "median")
  st <- select_diff_edges(g, rule = "structural", decision = "median")
  expect_true(all(flagged(mg, d$A)) && all(flagged(mg, d$B)))
  expect_false(any(flagged(mg, d$C)))
  expect_false(any(flagged(st, d$B)))
})

# --- 4. Multi-chain pooling -------------------------------------------------

test_that("pooled probabilities are the n_post-weighted average of the chains", {
  f1 <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 20, nsave = 30,
                      seed = 1, save_draws = FALSE)
  f2 <- multiggm_mcmc(data_list = sim2$data_list, method = "ssvs", burnin = 20, nsave = 70,
                      seed = 2, save_draws = FALSE)
  fl <- structure(list(chains = list(f1, f2), K = 2, nchains = 2), class = "multiggm_fit_list")
  e1 <- edge_probs(f1); e2 <- edge_probs(f2); ep <- edge_probs(fl)
  w <- c(f1$n_post, f2$n_post) / (f1$n_post + f2$n_post)
  # p_sign is a maximum of two probabilities, so it is pooled from the counts,
  # not averaged; the other columns are linear in the counts.
  for (col in c("mean", "p_struct", "p_exceed_0.05", "p_union_0.1"))
    expect_equal(ep[[col]], w[1] * e1[[col]] + w[2] * e2[[col]], label = col)
  expect_equal(attr(ep, "n_draws"), 100)
  expect_equal(edge_probs(fl, chain = 2), e2, ignore_attr = TRUE)
  expect_equal(pip_edges(fl), w[1] * pip_edges(f1) + w[2] * pip_edges(f2))
})

# --- interface ----------------------------------------------------------------

test_that("edge_prob_matrix() and selection helpers behave as documented", {
  f <- multiggm_mcmc(data_list = sim3$data_list, method = "ssvs", burnin = 30, nsave = 40,
                     seed = 3, save_draws = FALSE)
  ep <- edge_probs(f)
  expect_setequal(unique(paste(ep$group1, ep$group2)[ep$type == "pair"]), c("1 2", "1 3", "2 3"))
  m12 <- edge_prob_matrix(f, "mean", pair = c(1, 2))
  expect_equal(edge_prob_matrix(f, "mean", pair = c(2, 1)), -m12)
  expect_true(isSymmetric(m12))
  expect_equal(edge_prob_matrix(f, "p_struct", group = 2), pip_edges(f)[, , 2])
  expect_error(edge_prob_matrix(f, "p_cut", pair = c(1, 2)), "not available")
  expect_error(edge_prob_matrix(f, "p_exceed", group = 1), "single threshold")

  s1 <- select_edges(f, group = 1, decision = "median")
  expect_s3_class(s1, "multiggm_edge_selection")
  expect_equal(s1$selected, pip_edges(f)[, , 1] > 0.5)
  expect_error(select_edges(f, rule = "cut"), "mGHS")
  expect_error(select_edges(f, rule = "union"), "select_diff_edges")
  expect_error(select_diff_edges(f, rule = "cut"), "mGHS")
  expect_output(print(select_diff_edges(f, groups = c(1, 3), rule = "union", c = 0.05)),
                "union")

  g <- multiggm_mcmc(data_list = sim2$data_list, method = "mghs", burnin = 30, nsave = 30,
                     seed = 3, save_draws = FALSE)
  expect_error(select_edges(g, rule = "median_shrinkage", decision = "bfdr"), "not a probability")
  ms <- select_edges(g, rule = "median_shrinkage")
  expect_equal(ms$decision, "median")
  expect_true(is.na(ms$expected_fdr))
  expect_length(attr(edge_probs(g), "t_hat"), 2)
  expect_error(pip_edges(g), "edge-inclusion indicators")
})

test_that("bayes_fdr keeps tied probabilities together", {
  # 0.99, 0.95 alone: mean(1 - p) = 0.03; adding the two tied 0.9s gives 0.065
  sel <- bayes_fdr(c(0.99, 0.95, 0.9, 0.9, 0.2), alpha = 0.05)
  expect_equal(as.logical(sel), c(TRUE, TRUE, FALSE, FALSE, FALSE))
  sel2 <- bayes_fdr(c(0.99, 0.95, 0.9, 0.9, 0.2), alpha = 0.07)
  expect_equal(as.logical(sel2), c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_equal(attr(sel2, "expected_fdr"), mean(c(0.01, 0.05, 0.1, 0.1)))
})
