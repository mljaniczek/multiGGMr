# Validation of the packaged mGHS sampler

Reference: upstream mGHS (github.com/cbusatto/mGHS, commit 1c9fd70),
installed from a cleaned copy, called as
`mGHS::mGHS(B, bn, n, S, p, hyp_ta, chain = 1L)`.

Packaged: `multiGGM:::mghs_engine_cpp(bn, B - bn, thin = 1, n, S, p, hyp_ta, ...)`,
started from the same `set.seed()`.

## Result: exact reproduction

The local changes do not alter the random-number stream (see UPSTREAM.md),
so the packaged sampler reproduces upstream draw for draw.

| Setting | Seeds | Precision off-diagonals | Diagonals | tau | R |
|---|---|---|---|---|---|
| K=2, p=8, n=50, B=300 | 2 | identical | identical | identical | max diff < 1e-15 |
| K=3, p=10, n=40, B=300 | 3 | identical | identical | identical | max diff < 1e-15 |
| K=2, p=15, n=100, B=200 | 4 | identical | identical | identical | max diff < 1e-15 |
| K=2, p=50, n=200, B=1500 (`simulate_multiggm`) | 1, 2, 3 | max abs diff 0 | | | |

The R differences come only from how each version writes the saved matrix:
upstream copies the upper triangle into the lower triangle, while the
package saves the sampled matrix as it is.

Expected difference: the cut-model selection probabilities (`z`) differ,
by up to 0.31 at p = 50, because of the indexing fix (UPSTREAM.md, item 4).
They are not used by any package function.

Seed 1 at K=2, p=8, n=50 got stuck in **both** versions. The packaged
version restarts a stuck chain automatically.

Run times were the same order of magnitude (13 to 21 s per fit at p = 50,
B = 1500 on an Apple-silicon laptop, both versions).

## Reproducing

`tests/testthat/test-mghs-upstream.R` runs the first three rows with
`expect_identical()`. It is skipped unless upstream mGHS is installed;
point `MULTIGGM_MGHS_LIB` at its library if needed:

```
NOT_CRAN=true MULTIGGM_MGHS_LIB=~/Rlib-dina-new Rscript -e 'devtools::test(filter = "mghs-upstream")'
```

`tests/testthat/test-mghs.R` checks that the running summaries kept with
`save_draws = FALSE` equal the draw-based summaries exactly, that chains pool
correctly, and that non-finite input stops at once instead of hanging.
