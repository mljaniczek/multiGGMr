# Choosing the minimum-difference threshold c for the magnitude rule from a
# permutation (label-shuffled) null. Each permutation is a full refit with the
# analysis settings; only its call counts are kept.

# MCMC settings of a fit: from fit$settings (multiGGM >= 0.9.0), otherwise
# reconstructed from the chains. burnin is NA when it cannot be recovered.
.settings_from_fit <- function(fit) {
  if (!is.null(fit$settings)) return(fit$settings)
  chains <- .fit_chains(fit)
  ch <- chains[[1]]
  method <- ch$method
  tr <- if (identical(method, "mghs")) ch$R_save else ch$Theta_save
  nsave <- if (!is.null(tr)) dim(tr)[3] else NA_integer_
  thin <- if (is.finite(nsave) && !is.null(ch$n_post)) ch$n_post / nsave else NA_real_
  hyper <- ch$hyper
  if (!is.null(hyper)) {
    keep <- names(.default_hyper(method, ch$p))
    hyper <- hyper[intersect(names(hyper), keep)]
  }
  burnin <- NA_real_
  if (inherits(fit, "multiggm_fit_list") && !is.null(fit$call)) {
    b <- fit$call$burnin
    if (is.null(b)) burnin <- eval(formals(multiggm_mcmc)$burnin)
    else if (is.numeric(b)) burnin <- b
  }
  list(method = method, burnin = burnin, nsave = nsave, thin = thin,
       nchains = length(chains), hyper = hyper, engine = "cpp", dots = list())
}

# Permutation i: indices of the pooled rows (group a then group b) assigned to
# group a, the rest to group b; and the seed for its fit. All drawn up front
# from `seed`, so results do not depend on batching or parallel workers.
.make_perms <- function(n_a, n_b, n_perm, seed) {
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old, globalenv()))
  set.seed(seed)
  lapply(seq_len(n_perm), function(i)
    list(idx_a = sort(sample.int(n_a + n_b, n_a)),
         seed = sample.int(.Machine$integer.max, 1L)))
}

# Data of permutation `perm`: each original group is centred on its own means
# (as multiggm_mcmc() does), then groups a and b are pooled and relabelled.
.perm_data <- function(data_list, groups, perm) {
  X <- lapply(data_list, function(x) {
    x <- as.matrix(x); x <- scale(x, center = TRUE, scale = FALSE); matrix(x, nrow(x))
  })
  pool <- rbind(X[[groups[1]]], X[[groups[2]]])
  X[[groups[1]]] <- pool[perm$idx_a, , drop = FALSE]
  X[[groups[2]]] <- pool[-perm$idx_a, , drop = FALSE]
  X
}

# Fit one permutation and return its call counts on the grid. Errors are
# returned, not thrown.
.calibrate_one <- function(i, perm, data_list, groups, settings, c_grid, alpha) {
  t0 <- proc.time()[["elapsed"]]
  warn <- character(0)
  res <- tryCatch(withCallingHandlers({
    X <- .perm_data(data_list, groups, perm)
    args <- c(list(data_list = X, method = settings$method, burnin = settings$burnin,
                   nsave = settings$nsave, thin = settings$thin,
                   nchains = settings$nchains, hyper = settings$hyper,
                   engine = settings$engine, seed = perm$seed, parallel = FALSE),
              settings$dots)
    args$save_draws <- FALSE
    args$diff_deltas <- c_grid
    args$mcse_batches <- 0L
    args$waic <- FALSE
    fit <- do.call(multiggm_mcmc, args)
    acc <- .pooled_edge_acc(.fit_chains(fit))
    g <- which(acc$group_pairs[, 1] == min(groups) & acc$group_pairs[, 2] == max(groups))
    n <- acc$n
    pr <- acc$pair$pcor$n_exceed[, , g, drop = FALSE] / n
    calls <- vapply(seq_along(c_grid), function(t)
      attr(bayes_fdr(pr[, t, 1], alpha = alpha), "n_selected"), numeric(1))
    list(ok = TRUE, calls = calls, abs_mean = abs(acc$pair$pcor$sum[, g] / n), error = NA_character_)
  }, warning = function(w) {
    warn <<- c(warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e)
    list(ok = FALSE, calls = rep(NA_real_, length(c_grid)), abs_mean = NULL,
         error = conditionMessage(e)))
  res$perm <- i
  res$seed <- perm$seed
  res$warning <- if (length(warn)) paste(unique(warn), collapse = " | ") else NA_character_
  res$seconds <- proc.time()[["elapsed"]] - t0
  res
}

# Smallest grid c at and above which `ok` holds at every grid value.
.first_stable <- function(c_grid, ok) {
  ok[is.na(ok)] <- FALSE
  bad <- which(!ok)
  if (!length(bad)) return(c_grid[1])
  if (max(bad) == length(c_grid)) return(NA_real_)
  c_grid[max(bad) + 1]
}

# Censored Weibull fit to the per-permutation "last call" thresholds T_i
# (largest grid c with any call). T_i = 0 (no call anywhere) is left-censored
# at the smallest grid value, T_i = max(grid) is right-censored there, and
# other values are known only to within one grid step (interval-censored on
# [T, T + step)). Returns the (1 - alpha) quantile of T, i.e. the c with
# P(any call) = alpha, and a delta-method interval. In permutation studies at
# p = 10-30 the Weibull gave an achieved null rate close to alpha, where the
# lognormal was far too conservative and the normal and the raw empirical
# rule too liberal.
.fit_last_call <- function(T, c_grid, alpha, level = 0.95) {
  lo <- c_grid[1]; hi <- c_grid[length(c_grid)]
  left <- T < lo; right <- T >= hi; obs <- !left & !right
  steps <- c(diff(c_grid), diff(c_grid)[length(c_grid) - 1L])
  step <- steps[match(round(T, 12), round(c_grid, 12))]
  na <- list(c = NA_real_, lower = NA_real_, upper = NA_real_, loc = NA_real_,
             scale = NA_real_, converged = FALSE)
  if (sum(obs) + sum(right) < 3L) {
    if (all(left)) na$c <- lo
    return(na)
  }
  # Weibull on the log scale: log T = loc + scale * W, W standard minimum
  # Gumbel; F(x) = 1 - exp(-exp((log x - loc) / scale)).
  lF <- function(x, m, s) log(-expm1(-exp((log(x) - m) / s)))
  lS <- function(x, m, s) -exp((log(x) - m) / s)
  nll <- function(th) {
    m <- th[1]; s <- exp(th[2])
    ll <- 0
    if (any(left)) ll <- ll + sum(left) * lF(lo, m, s)
    if (any(right)) ll <- ll + sum(right) * lS(hi, m, s)
    if (any(obs)) {
      a <- lS(T[obs], m, s); b <- lS(T[obs] + step[obs], m, s)
      ll <- ll + sum(a + log(-expm1(b - a)))
    }
    if (!is.finite(ll)) return(1e10)
    -ll
  }
  x <- log(pmax(T, lo))
  op <- stats::optim(c(mean(x), log(max(stats::sd(x), 0.05))), nll, method = "BFGS",
                     hessian = TRUE)
  w <- log(-log(alpha))
  est <- op$par[1] + exp(op$par[2]) * w
  V <- tryCatch(solve(op$hessian), error = function(e) matrix(NA_real_, 2, 2))
  grad <- c(1, exp(op$par[2]) * w)
  se <- sqrt(max(drop(t(grad) %*% V %*% grad), 0))
  zz <- stats::qnorm(1 - (1 - level) / 2)
  list(c = exp(est), lower = exp(est - zz * se), upper = exp(est + zz * se),
       loc = op$par[1], scale = exp(op$par[2]), converged = op$convergence == 0)
}

# P(T >= c) under the fitted Weibull.
.p_any_fitted <- function(fitT, c_grid) {
  if (is.null(fitT) || !is.finite(fitT$loc)) return(rep(NA_real_, length(c_grid)))
  exp(-exp((log(c_grid) - fitT$loc) / fitT$scale))
}

.default_c_grid <- function(n_a, n_b) {
  cmax <- max(0.40, ceiling(5 * sqrt(1 / n_a + 1 / n_b) * 20) / 20)
  round(seq(0.01, min(cmax, 0.95), by = 0.01), 2)
}

#' Choose the minimum-difference threshold c from a permutation null
#'
#' Estimates the smallest threshold \eqn{c} for the magnitude rule of
#' [select_diff_edges()] (\eqn{P(|D| > c)} with Bayesian FDR control) at
#' which these data, analysed with these MCMC settings, make no false calls
#' when the two groups do not differ. The group labels are shuffled
#' \code{n_perm} times; each shuffled dataset is refitted with the settings of
#' the real analysis, and every call it makes is false.
#'
#' @param data_list The list of data matrices of the real analysis.
#' @param fit The fit of the real analysis (recommended). Its method, prior
#'   (\code{hyper}), \code{burnin}, \code{nsave}, \code{thin} and number of
#'   chains are reused, so the null is fitted exactly as the real data. Fits
#'   made before multiGGM 0.9.0 do not record \code{burnin} reliably; give it
#'   in \code{mcmc_args} if asked.
#' @param method Sampler when \code{fit} is not given: \code{"ssvs"},
#'   \code{"mghs"} or \code{"gwishart"}.
#' @param mcmc_args Named list of [multiggm_mcmc()] arguments. Without
#'   \code{fit}, it must give at least \code{burnin} and \code{nsave}; with
#'   \code{fit}, its entries override the fit's settings.
#'   \code{save_draws}, \code{diff_deltas}, \code{mcse_batches}, \code{waic},
#'   \code{seed} and \code{parallel} are set by \code{calibrate_c()}.
#' @param groups The two groups to calibrate. Only these two groups' samples
#'   are shuffled; any other groups are left as they are.
#' @param alpha Target FDR of the selection rule, as in
#'   [select_diff_edges()]. Also the default target for \code{"any_call"}.
#' @param target \code{"any_call"} (default): the share of null datasets
#'   with at least one call is at most \code{target_value} (default
#'   \code{alpha}). Under a null of no differences every call is false, so this
#'   is the FDR of the rule. \code{"mean_false_calls"}: the expected number
#'   of false calls per dataset is at most \code{target_value} (default 1),
#'   a weaker per-family error rate.
#' @param target_value See \code{target}.
#' @param c_grid Thresholds on the partial-correlation scale. Default
#'   \code{0.01} to \code{max(0.40, 5 * sqrt(1/n_1 + 1/n_2))} by 0.01, where
#'   \eqn{\sqrt{1/n_1 + 1/n_2}} is about the sampling SD of a difference in
#'   partial correlations.
#' @param n_perm Maximum number of permutations (default 50).
#' @param min_perm Permutations before the first stopping check (default 20).
#' @param c_tol Stop early once the 95% confidence interval for the
#'   estimated \eqn{c} is within \eqn{\pm} \code{c_tol} (default 0.02). Only for
#'   \code{target = "any_call"}; \code{c_tol = 0} always runs all
#'   \code{n_perm}.
#' @param seed Seed for the permutations and the fits' seeds.
#' @param parallel,ncores Fit permutations in parallel (forked processes on
#'   Unix-alikes; sequential on Windows). Each permutation's chains run
#'   sequentially inside its worker.
#' @param verbose Print progress.
#'
#' @details
#' \strong{What the null is.} Each original group is first centred on its
#' own means, exactly as [multiggm_mcmc()] centres the data it fits, so
#' differences in mean expression between groups do not leak into the
#' shuffled groups' covariances. The pooled samples of the two groups are
#' then reassigned at random with the original group sizes. Nothing is
#' rescaled. The null is therefore "the two groups share one covariance
#' structure", and \eqn{c_{null}} is the smallest difference in partial
#' correlation these data can tell from noise.
#'
#' \strong{Same settings as the analysis.} The posterior spread of \eqn{D},
#' and the resolution \eqn{1 / n_{save}} of its probabilities, depend on the
#' prior and the chain length, so a null fitted with other settings gives a
#' \eqn{c} that does not transfer. Use the fit's settings, and for mGHS the
#' same \code{hyper$R_eta}. Under a permutation null the groups share one
#' network, which is exactly the case where mGHS with the uniform prior on
#' \eqn{R} runs \eqn{R} to \eqn{\pm 1} at small \eqn{n} (see
#' [multiggm_mcmc()]): use \code{R_eta > 1} there, in the analysis too.
#'
#' \strong{Choosing c.} For each permutation, \eqn{T_i} is the largest grid
#' value of \eqn{c} at which it makes any call. A dataset makes a call at
#' \eqn{c} only if \eqn{T_i \ge c}, so the share of null datasets with a call
#' at \eqn{c} is \eqn{P(T \ge c)}. With \code{target = "any_call"}, a
#' Weibull distribution is fitted to the \eqn{T_i} (censored at the ends of
#' the grid and to the grid spacing), and \eqn{c_{null}} is its
#' \eqn{1 - } \code{target_value} quantile, with a delta-method confidence
#' interval. This borrows strength across permutations, so 20 to 50 suffice
#' where the raw share at one \eqn{c} would need hundreds. In permutation
#' studies (SSVS, p = 10 and 30, n = 70 to 400 per group, references from
#' 200 permutations), the \eqn{c} chosen from 20 permutations gave an
#' achieved share of null datasets with a call of 0.03 to 0.05 on average
#' (target 0.05), and above 0.10 in at most 9% of runs; the raw empirical
#' rule exceeded 0.10 in up to 37%. With 40 permutations the estimate of
#' \eqn{c} varied by about \eqn{\pm 0.02} to \eqn{0.035} (10th to 90th
#' percentile). \code{suggested_c} is the smallest grid
#' value at or above \eqn{c_{null}}, so it can be used directly with a fit
#' whose \code{diff_deltas} contains the grid. The raw share and the mean
#' number of calls at every \eqn{c} are reported as well. With
#' \code{target = "mean_false_calls"}, \code{suggested_c} is the smallest grid
#' value from which the mean number of calls stays at or below
#' \code{target_value}.
#' When fewer than three permutations make a call, the tail cannot be fitted
#' and \code{suggested_c} is the smallest grid value from which the raw share
#' stays at or below the target (\code{c_method = "empirical"}).
#'
#' \strong{How to use it.} \eqn{c_{null}} is the smallest difference these
#' data can distinguish from noise, not the smallest difference that matters
#' scientifically. Use \eqn{c = \max(c_{scientific}, c_{null})}. Fit the real
#' analysis with \code{diff_deltas} containing the grid (or at least the
#' chosen \eqn{c}), so that [select_diff_edges()] can use it without saved
#' draws.
#'
#' \strong{Cost.} Up to \code{n_perm} fits of the full analysis. Fits
#' are independent, so wall time divides by the number of cores. Early
#' stopping ends after 20 to 30 permutations when the null curve is steep (larger
#' \eqn{n}); at small \eqn{n} expect all \code{n_perm}. As a rough guide, on
#' an Apple M-series laptop (plugged in) SSVS took about 25 to 35 ms per
#' iteration at p = 100, so a fit of 5,000 + 1,000 iterations takes about
#' 3 min. mGHS speed depends strongly on the data, through its rejection
#' samplers: at p of about 30 we measured from 5 to over 70 ms per iteration
#' on different datasets. Time one fit of the analysis first; with
#' \code{verbose = TRUE} the time per permutation and an upper bound on the
#' remaining time are printed after the first batch. Power-saving modes can
#' slow fits by an order of magnitude.
#'
#' Failed fits (for example an mGHS chain stopped at the boundary of
#' \eqn{R}) are recorded in \code{perms} and skipped. Warnings from fits
#' (such as a restarted mGHS chain) are recorded there too.
#'
#' @return An object of class \code{"multiggm_c_calibration"}: a list with
#'   \describe{
#'     \item{\code{summary}}{Data frame, one row per grid value:
#'       \code{c}, \code{mean_calls}, \code{se_calls}, \code{share_any},
#'       \code{share_any_lower}, \code{share_any_upper} (Clopper-Pearson
#'       95%), \code{p_any_fitted} (the fitted Weibull \eqn{P(T \ge c)}) and
#'       \code{meets_target}. Also returned by \code{as.data.frame()}.}
#'     \item{\code{calls}}{Matrix of call counts, permutations x grid.}
#'     \item{\code{suggested_c}}{The chosen grid value (\code{NA} if no grid
#'       value meets the target).}
#'     \item{\code{c_null}, \code{c_null_lower}, \code{c_null_upper}}{For
#'       \code{"any_call"}: the continuous estimate and its 95% interval.}
#'     \item{\code{last_call}}{The \eqn{T_i} of the successful permutations.}
#'     \item{\code{noise_ref}}{0.99 quantile of the posterior mean
#'       \eqn{|D|} over edges and permutations, a quick noise scale.}
#'     \item{\code{perms}}{Data frame, one row per permutation run:
#'       \code{perm}, \code{seed}, \code{ok}, \code{error}, \code{warning},
#'       \code{seconds}. With \code{seed}, \code{groups} and the original
#'       data, any permutation can be rebuilt (see Examples).}
#'     \item{\code{n_perm_run}, \code{n_ok}, \code{stopped}}{How many
#'       permutations were run and succeeded, and why it stopped.}
#'     \item{\code{target}, \code{target_value}, \code{alpha}, \code{groups},
#'       \code{settings}, \code{seed}}{Settings.}
#'   }
#'
#' @examples
#' \donttest{
#' # Two groups of different sizes that share one network
#' sim <- simulate_multiggm(K = 2, p = 8, n = 80, seed = 1)
#' set.seed(2)
#' X <- list(sim$data_list[[1]],
#'           matrix(rnorm(60 * 8), 60) %*% chol(solve(sim$Omega_list[[1]])))
#' fit <- multiggm_mcmc(data_list = X, method = "ssvs", burnin = 200, nsave = 200,
#'                      seed = 1)
#' cal <- calibrate_c(X, fit = fit, n_perm = 10)
#' cal
#' plot(cal)
#' head(as.data.frame(cal))
#' }
#'
#' @seealso [select_diff_edges()], [fdr_calibration()]
#' @export
calibrate_c <- function(data_list, fit = NULL, method = c("ssvs", "mghs", "gwishart"),
                        mcmc_args = list(), groups = c(1, 2), alpha = 0.05,
                        target = c("any_call", "mean_false_calls"), target_value = NULL,
                        c_grid = NULL, n_perm = 50, min_perm = 20, c_tol = 0.02,
                        seed = 1, parallel = FALSE,
                        ncores = max(1L, parallel::detectCores() - 1L), verbose = TRUE) {
  target <- match.arg(target)
  if (is.null(target_value)) target_value <- if (target == "any_call") alpha else 1
  if (!is.list(data_list) || length(data_list) < 2L)
    stop("data_list must be a list of at least two data matrices.", call. = FALSE)
  K <- length(data_list)
  groups <- .check_groups(groups, K)
  n_a <- nrow(as.matrix(data_list[[groups[1]]])); n_b <- nrow(as.matrix(data_list[[groups[2]]]))
  n_perm <- as.integer(n_perm); min_perm <- min(as.integer(min_perm), n_perm)
  if (n_perm < 2L) stop("n_perm must be at least 2.", call. = FALSE)

  # Settings: from the fit, overridden by mcmc_args
  if (!is.null(fit)) {
    if (inherits(fit, "multiggm_platform_fit"))
      stop("calibrate_c() does not support method = \"ssvs_platform\".", call. = FALSE)
    settings <- .settings_from_fit(fit)
    if (.fit_chains(fit)[[1]]$K != K)
      stop("data_list has ", K, " groups but fit has ", .fit_chains(fit)[[1]]$K, ".",
           call. = FALSE)
  } else {
    settings <- list(method = match.arg(method), burnin = NA_real_, nsave = NA_real_, thin = 1,
                     nchains = 1L, hyper = NULL, engine = "cpp", dots = list())
  }
  core <- c("method", "burnin", "nsave", "thin", "nchains", "hyper", "engine")
  for (nm in intersect(names(mcmc_args), core)) settings[[nm]] <- mcmc_args[[nm]]
  extra <- mcmc_args[setdiff(names(mcmc_args), core)]
  ignored <- intersect(names(extra), c("save_draws", "diff_deltas", "mcse_batches", "waic",
                                       "seed", "parallel", "ncores", "data_list", "S_list",
                                       "n_vec", "prior_only"))
  if (length(ignored))
    warning("Ignoring mcmc_args set by calibrate_c(): ", paste(ignored, collapse = ", "),
            call. = FALSE)
  settings$dots[names(extra)] <- extra
  settings$dots[ignored] <- NULL
  settings$dots[c("save_draws", "diff_deltas", "waic", "prior_only")] <- NULL
  if (!settings$method %in% c("ssvs", "mghs", "gwishart"))
    stop("calibrate_c() supports method = \"ssvs\", \"mghs\" or \"gwishart\".", call. = FALSE)
  for (nm in c("burnin", "nsave"))
    if (length(settings[[nm]]) != 1L || is.na(settings[[nm]]))
      stop(if (is.null(fit)) paste0("Give ", nm, " in mcmc_args (the value used for the analysis).")
           else paste0("The fit does not record ", nm, " (fits made before multiGGM 0.9.0 ",
                       "with a variable or a single chain); give it in mcmc_args, e.g. ",
                       "mcmc_args = list(", nm, " = 5000)."), call. = FALSE)
  if (identical(settings$method, "mghs") && min(n_a, n_b) <= 200 &&
      (is.null(settings$hyper$R_eta) || settings$hyper$R_eta == 1))
    warning("mGHS with the uniform prior on R (R_eta = 1) and ", min(n_a, n_b),
            " samples per group: under the permutation null the groups share one network, ",
            "where the chain often stops at the boundary of R. Consider hyper$R_eta = 10, in ",
            "the analysis as well.", call. = FALSE)

  if (is.null(c_grid)) c_grid <- .default_c_grid(n_a, n_b)
  c_grid <- sort(unique(round(as.numeric(c_grid), 12)))
  if (length(c_grid) < 2L || any(c_grid <= 0) || any(c_grid >= 1))
    stop("c_grid must contain at least two values in (0, 1) (partial-correlation scale).",
         call. = FALSE)

  if (verbose) {
    cat(sprintf("calibrate_c: %s, burnin %g, nsave %g, thin %g, %d chain(s); groups %d (n = %d) vs %d (n = %d)\n",
                settings$method, settings$burnin, settings$nsave, settings$thin,
                as.integer(settings$nchains), groups[1], n_a, groups[2], n_b))
    if (length(settings$hyper))
      cat("  hyper:", paste(names(settings$hyper), vapply(settings$hyper, function(v)
        if (length(v) == 1L) format(v) else "<matrix>", ""), sep = " = ", collapse = ", "), "\n")
    cat(sprintf("  up to %d permutations; c grid %g to %g (%d values)\n", n_perm,
                min(c_grid), max(c_grid), length(c_grid)))
  }

  perms <- .make_perms(n_a, n_b, n_perm, seed)
  use_par <- isTRUE(parallel) && .Platform$OS.type == "unix" && ncores > 1L
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) ncores <- min(ncores, 2L)
  batch <- if (use_par) as.integer(ncores) else 5L
  results <- list()
  t_start <- proc.time()[["elapsed"]]
  stopped <- "n_perm reached"
  est <- NULL
  done <- 0L
  while (done < n_perm) {
    ids <- seq(done + 1L, min(done + batch, n_perm))
    run <- function(i) .calibrate_one(i, perms[[i]], data_list, groups, settings, c_grid, alpha)
    out <- if (use_par) parallel::mclapply(ids, run, mc.cores = ncores, mc.preschedule = FALSE)
           else lapply(ids, run)
    out <- lapply(seq_along(ids), function(k) {
      o <- out[[k]]
      if (inherits(o, "try-error") || !is.list(o))
        list(ok = FALSE, calls = rep(NA_real_, length(c_grid)), abs_mean = NULL,
             error = as.character(o), perm = ids[k], seed = perms[[ids[k]]]$seed,
             warning = NA_character_, seconds = NA_real_)
      else o
    })
    results <- c(results, out)
    done <- max(ids)
    ok <- vapply(results, `[[`, logical(1), "ok")
    if (verbose)
      cat(sprintf("  %d/%d permutations run, %d succeeded%s\n", done, n_perm, sum(ok),
                  if (length(results) == length(ids)) {
                    el <- proc.time()[["elapsed"]] - t_start
                    sprintf(" (%.0f s per batch of %d; at most %.0f min remaining)", el,
                            length(ids), el * ceiling((n_perm - done) / batch) / 60)
                  } else ""))
    if (target == "any_call" && done >= min_perm && sum(ok) >= 5L && done < n_perm) {
      Tl <- .last_calls(results[ok], c_grid)
      est <- .fit_last_call(Tl, c_grid, target_value)
      if (is.finite(est$lower) && is.finite(est$upper) &&
          max(est$upper - est$c, est$c - est$lower) <= c_tol) {
        stopped <- sprintf("confidence interval within +-%g", c_tol)
        break
      }
    }
  }

  ok <- vapply(results, `[[`, logical(1), "ok")
  perms_df <- data.frame(
    perm = vapply(results, `[[`, numeric(1), "perm"),
    seed = vapply(results, function(r) as.numeric(r$seed), numeric(1)),
    ok = ok,
    error = vapply(results, function(r) as.character(r$error), character(1)),
    warning = vapply(results, function(r) as.character(r$warning), character(1)),
    seconds = vapply(results, function(r) as.numeric(r$seconds), numeric(1)),
    stringsAsFactors = FALSE)
  if (sum(ok) == 0L)
    stop("All ", length(ok), " permutation fits failed. First error: ", perms_df$error[1],
         call. = FALSE)
  if (sum(ok) < 0.8 * length(ok))
    warning(sprintf("Only %d of %d permutation fits succeeded; see $perms.", sum(ok), length(ok)),
            call. = FALSE)

  calls <- do.call(rbind, lapply(results, `[[`, "calls"))
  dimnames(calls) <- list(perm = perms_df$perm, c = format(c_grid))
  cok <- calls[ok, , drop = FALSE]
  n_ok <- nrow(cok)
  any_call <- cok > 0
  k_any <- colSums(any_call)
  ci <- t(vapply(k_any, function(k) {
    if (n_ok == 0) return(c(NA_real_, NA_real_))
    b <- stats::binom.test(k, n_ok)$conf.int
    c(b[1], b[2])
  }, numeric(2)))
  Tl <- .last_calls(results[ok], c_grid)
  fitT <- if (n_ok >= 5L) .fit_last_call(Tl, c_grid, target_value) else NULL
  p_fit <- .p_any_fitted(fitT, c_grid)
  summ <- data.frame(c = c_grid, mean_calls = colMeans(cok),
                     se_calls = apply(cok, 2, stats::sd) / sqrt(n_ok),
                     share_any = k_any / n_ok, share_any_lower = ci[, 1],
                     share_any_upper = ci[, 2], p_any_fitted = p_fit, row.names = NULL)

  c_null <- c_lo <- c_hi <- NA_real_
  c_method <- if (target == "any_call") NA_character_ else "mean_calls"
  if (target == "any_call") {
    if (n_ok < 5L) {
      warning("Fewer than 5 successful permutations: no c suggested.", call. = FALSE)
      suggested <- NA_real_
    } else if (!is.finite(fitT$c) || !is.finite(fitT$loc)) {
      # Too few permutations with a call to fit the tail: use the raw share.
      suggested <- .first_stable(c_grid, summ$share_any <= target_value)
      c_method <- "empirical"
    } else {
      c_null <- fitT$c; c_lo <- fitT$lower; c_hi <- fitT$upper
      suggested <- if (is.finite(c_null) && c_null <= max(c_grid))
        c_grid[which(c_grid >= c_null - 1e-12)[1]] else NA_real_
      c_method <- "weibull"
    }
    summ$meets_target <- if (is.finite(suggested)) summ$c >= suggested else FALSE
  } else {
    summ$meets_target <- summ$mean_calls <= target_value
    suggested <- .first_stable(c_grid, summ$meets_target)
  }
  if (!is.finite(suggested) && n_ok >= 5L)
    warning("No c in the grid meets the target; extend c_grid (largest value ",
            max(c_grid), ").", call. = FALSE)
  am <- unlist(lapply(results[ok], `[[`, "abs_mean"))

  structure(list(summary = summ, calls = calls, suggested_c = suggested,
                 c_null = c_null, c_null_lower = c_lo, c_null_upper = c_hi,
                 last_call = Tl, weibull = fitT[c("loc", "scale")], c_method = c_method,
                 noise_ref = if (length(am)) unname(stats::quantile(am, 0.99)) else NA_real_,
                 perms = perms_df, n_perm_run = nrow(perms_df), n_ok = n_ok, stopped = stopped,
                 target = target, target_value = target_value, alpha = alpha,
                 groups = groups, n = c(n_a, n_b), settings = settings, seed = seed,
                 call = match.call()),
            class = "multiggm_c_calibration")
}

# Largest grid c with any call per permutation (0 if none).
.last_calls <- function(res, c_grid) {
  vapply(res, function(r) {
    w <- which(r$calls > 0)
    if (length(w)) c_grid[max(w)] else 0
  }, numeric(1))
}

#' @export
print.multiggm_c_calibration <- function(x, ...) {
  cat("<multiggm_c_calibration>", x$settings$method, "| groups", x$groups[1], "vs",
      x$groups[2], sprintf("(n = %d, %d)", x$n[1], x$n[2]), "\n")
  cat(sprintf("  %d permutations run, %d succeeded (%s)\n", x$n_perm_run, x$n_ok, x$stopped))
  tgt <- if (x$target == "any_call")
    sprintf("share of null datasets with any call <= %g", x$target_value)
  else sprintf("mean false calls per null dataset <= %g", x$target_value)
  cat("  Target:", tgt, sprintf("(selection: P(|D| > c), Bayesian FDR %g, pcor scale)\n", x$alpha))
  if (x$target == "any_call" && is.finite(x$c_null))
    cat(sprintf("  c_null = %.3f (95%% CI %.3f to %.3f)\n", x$c_null, x$c_null_lower, x$c_null_upper))
  cat(sprintf("  Suggested c: %s%s\n",
              if (is.finite(x$suggested_c)) format(x$suggested_c) else "none in grid",
              if (identical(x$c_method, "empirical"))
                " (raw share: too few null calls to fit the tail)" else ""))
  cat(sprintf("  Noise reference (0.99 quantile of |E[D]| under the null): %.3f\n", x$noise_ref))
  if (is.finite(x$suggested_c)) {
    r <- x$summary[match(x$suggested_c, x$summary$c), ]
    cat(sprintf("  At c = %g: share with any call %.2f (raw), mean calls %.2f\n",
                x$suggested_c, r$share_any, r$mean_calls))
  }
  if (any(!x$perms$ok)) cat(sprintf("  %d permutation fit(s) failed; see $perms.\n", sum(!x$perms$ok)))
  cat("  Use c = max(c_scientific, suggested_c). as.data.frame() gives the curve.\n")
  invisible(x)
}

#' @export
as.data.frame.multiggm_c_calibration <- function(x, ...) x$summary

#' @rdname calibrate_c
#' @param x A \code{"multiggm_c_calibration"} object.
#' @param type \code{"any"} (default): share of null datasets with any call
#'   against \eqn{c}, with Clopper-Pearson intervals, the fitted curve, the
#'   target and the chosen \eqn{c}. \code{"mean"}: mean number of false
#'   calls per dataset, with its SE.
#' @param ... Unused.
#' @export
plot.multiggm_c_calibration <- function(x, type = c("any", "mean"), ...) {
  type <- match.arg(type)
  s <- x$summary
  if (type == "any") {
    g <- ggplot2::ggplot(s, ggplot2::aes(x = .data$c)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$share_any_lower,
                                          ymax = .data$share_any_upper),
                             width = 0, color = "grey70") +
      ggplot2::geom_point(ggplot2::aes(y = .data$share_any), size = 1.2)
    if (!all(is.na(s$p_any_fitted)))
      g <- g + ggplot2::geom_line(ggplot2::aes(y = .data$p_any_fitted), color = "steelblue")
    if (x$target == "any_call")
      g <- g + ggplot2::geom_hline(yintercept = x$target_value, linetype = "dashed")
    ylab <- "Share of null datasets with any call"
  } else {
    g <- ggplot2::ggplot(s, ggplot2::aes(x = .data$c, y = .data$mean_calls)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(.data$mean_calls - 2 * .data$se_calls, 0),
                                        ymax = .data$mean_calls + 2 * .data$se_calls),
                           fill = "grey85") +
      ggplot2::geom_line()
    if (x$target == "mean_false_calls")
      g <- g + ggplot2::geom_hline(yintercept = x$target_value, linetype = "dashed")
    ylab <- "Mean false calls per null dataset"
  }
  if (is.finite(x$suggested_c))
    g <- g + ggplot2::geom_vline(xintercept = x$suggested_c, color = "firebrick")
  if (x$target == "any_call" && is.finite(x$c_null_lower) && type == "any")
    g <- g + ggplot2::annotate("rect", xmin = x$c_null_lower, xmax = x$c_null_upper,
                               ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "firebrick")
  g + ggplot2::labs(title = sprintf("Permutation null (%d datasets): suggested c = %s",
                                    x$n_ok, format(x$suggested_c)),
                    x = "c (partial correlation)", y = ylab) +
    ggplot2::theme_minimal()
}
