# mGHS sampler: provenance and local changes

The `method = "mghs"` sampler is adapted from the R package **mGHS** by
Claudio Busatto, <https://github.com/cbusatto/mGHS>, commit `1c9fd70`,
licensed GPL (>= 2) and redistributed here under GPL-3.

Method: Busatto, C. & Stingo, F.C. (2025). Inference of multiple
high-dimensional networks with the Graphical Horseshoe prior. *Journal of
Computational and Graphical Statistics*. doi:10.1080/10618600.2024.2435480

## Files

Imported unchanged in commit "Import mGHS sampler sources verbatim" (run
`git log -- src/mGHS.cpp` to see every later change):

`mGHS.cpp`, `rgamma3p.cpp`, `utils_mGHS.cpp`, `utils_mGHS.h`,
`utils_rgamma3p.h`, `utils_rIW.h`, `rIW.cpp`, `pbdv.cpp`, `pbdv.h`,
`DVLA.cpp/.h`, `DVSA.cpp/.h`, `VVLA.cpp/.h`, `Gamma2.cpp/.h`,
`factorial.cpp/.h`.

Not imported: `rgamma3p_tabulated.cpp` (`rg3p()`, only referenced in
commented-out code), the `R/` analysis scripts, and the compiled `.o`/`.dll`
files.

## Changes

1. **Build integration.** Removed `Rcpp::export` attributes from upstream's
   helpers; removed `lambda_pstep()` and its RcppParallel worker (never
   called, and it would call R's RNG from worker threads); removed
   declarations of functions that were not imported; added an attribution
   header to each file.
2. **Rejection samplers can no longer hang** (`rgamma3p.cpp`). Each loop
   stops after `MGHS_MAX_TRIES` (1e6) tries and returns `NaN`; `NaN`/`Inf`
   arguments return `NaN` at once. Previously a `NaN` argument sent
   `rg3p_c1()` to an all-zero constant table that never accepts, the
   `while (i < 1e+16)` loop had an `int` counter that overflows first, and the
   truncated-normal `while (out <= 0)` loops never ended when the mean was
   far below zero.
3. **Stuck checks catch non-finite values** (`mGHS.cpp`). Upstream tested
   `x > 1e+06`, which is false for `NaN`. The engine now also checks every
   local and global scale draw. A stuck chain returns early with
   `stuck = 1`, and the R wrapper restarts it (`max_attempts`, default 3).
4. **Cut-model indexing bug** (`mGHS.cpp`). Upstream read
   `Lambda[js][j][hs]` after `hs` had been left equal to `K` by an earlier
   loop. With bounds checks off this reads `Lambda[js][j+1][0]`: every
   group's cut-model step used group 1's local scales shifted by one column,
   and pairs `(j, j+1)` read the constant diagonal (1.0). Fixed to
   `Lambda[js][j][h]`. The cut model does not feed back into the other
   parameters, and `rbinom(1, q)` uses one uniform for any `0 < q < 1`, so
   the random-number stream and the precision draws are unchanged.
5. **Engine.** `mGHS_chain()` became `mghs_engine_cpp()`; `mGHS_mean()` and
   the `mGHS()` dispatcher were removed. The sampling steps are unchanged.
   New: thinning; draws stored as an R array `[p, p, K, nsave]` (fixes the
   row-wise pair order at the boundary: upstream stores pairs as (1,2),
   (1,3), ..., (1,p), (2,3), ..., which is R's `lower.tri()` order, not
   `upper.tri()`); running sums over every post-burn-in iteration of the
   precision and partial-correlation means and of group differences
   (`save_draws = FALSE` keeps only these); `checkUserInterrupt()` each
   iteration; progress output only when `verbose`; a stuck chain reports
   which check fired (`stuck_reason`); the number of extra rejection-sampler
   tries is counted (`reject_tries`, a diagnostic for slow fits).

Not changed: the cut model's use of `Lambda / (1 + Lambda)` without the
global scale `tau`. The package does not expose cut-model edge selection
(`cut_z_mean` is stored for reference only).
