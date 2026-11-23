# Package index

## All functions

- [`bayes_fdr()`](https://mljaniczek.github.io/multiGGMr/reference/bayes_fdr.md)
  : Bayesian false discovery rate selection
- [`bfdr_qvalues()`](https://mljaniczek.github.io/multiGGMr/reference/bfdr_qvalues.md)
  : Bayesian FDR q-values
- [`calc_mrf_logC()`](https://mljaniczek.github.io/multiGGMr/reference/calc_mrf_logC.md)
  : Log normalizing constant for the MRF prior (per edge)
- [`calc_mrf_logC_cpp()`](https://mljaniczek.github.io/multiGGMr/reference/calc_mrf_logC_cpp.md)
  : Compute log normalizing constant for the edge-sharing MRF prior
- [`calibrate_c()`](https://mljaniczek.github.io/multiGGMr/reference/calibrate_c.md)
  [`plot(`*`<multiggm_c_calibration>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/calibrate_c.md)
  : Choose the minimum-difference threshold c from a permutation null
- [`coef(`*`<multiggm_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/coef.multiggm_fit.md)
  [`coef(`*`<multiggm_fit_list>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/coef.multiggm_fit.md)
  : Extract posterior mean precision matrices
- [`coef(`*`<multiggm_platform_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/coef.multiggm_platform_fit.md)
  : Extract posterior mean precision matrices from a multi-platform fit
- [`compute_waic()`](https://mljaniczek.github.io/multiGGMr/reference/compute_waic.md)
  : WAIC from running sums
- [`diff_count_posterior()`](https://mljaniczek.github.io/multiGGMr/reference/diff_count_posterior.md)
  : Posterior distribution of the number of differential edges
- [`diff_draws()`](https://mljaniczek.github.io/multiGGMr/reference/diff_draws.md)
  : Posterior draws of the difference between two groups
- [`diff_prob_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/diff_prob_pcor.md)
  : Posterior probability of differential partial correlation (K=2)
- [`diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/diff_summary.md)
  : Per-edge posterior summary of a group difference
- [`edge_bayes_factors()`](https://mljaniczek.github.io/multiGGMr/reference/edge_bayes_factors.md)
  : Prior versus posterior edge probabilities and Bayes factors
- [`edge_counts()`](https://mljaniczek.github.io/multiGGMr/reference/edge_counts.md)
  : Per-iteration edge counts for convergence checking
- [`edge_prob_matrix()`](https://mljaniczek.github.io/multiGGMr/reference/edge_prob_matrix.md)
  : One posterior edge probability as a p x p matrix
- [`edge_probs()`](https://mljaniczek.github.io/multiGGMr/reference/edge_probs.md)
  : Posterior edge probabilities for every group and pair of groups
- [`fdr_calibration()`](https://mljaniczek.github.io/multiGGMr/reference/fdr_calibration.md)
  : Calibration of Bayesian FDR selection against a known truth
- [`fitted(`*`<multiggm_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_fit.md)
  [`fitted(`*`<multiggm_fit_list>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_fit.md)
  : Extract posterior mean partial correlation matrices
- [`fitted(`*`<multiggm_platform_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/fitted.multiggm_platform_fit.md)
  : Extract posterior mean partial correlation matrices from a
  multi-platform fit
- [`fix_matrix()`](https://mljaniczek.github.io/multiGGMr/reference/fix_matrix.md)
  : Rescale a symmetric matrix toward positive definiteness
- [`mask_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/mask_pcor.md)
  : Mask small partial correlations to zero
- [`mcmc_diagnostics()`](https://mljaniczek.github.io/multiGGMr/reference/mcmc_diagnostics.md)
  : MCMC convergence diagnostics
- [`pip_table(`*`<multiggm_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/multiggm-accessors.md)
  [`pip_table(`*`<multiggm_fit_list>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/multiggm-accessors.md)
  [`pip_table(`*`<multiggm_platform_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/multiggm-accessors.md)
  [`pcor_estimate(`*`<multiggm_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/multiggm-accessors.md)
  [`pcor_estimate(`*`<multiggm_fit_list>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/multiggm-accessors.md)
  [`pcor_estimate(`*`<multiggm_platform_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/multiggm-accessors.md)
  : Shared-format accessors for multiGGMr fits
- [`multiggm_mcmc()`](https://mljaniczek.github.io/multiGGMr/reference/multiggm_mcmc.md)
  : Fit multiple Gaussian graphical models via MCMC
- [`node_diff_summary()`](https://mljaniczek.github.io/multiGGMr/reference/node_diff_summary.md)
  : Node-level differential network summaries
- [`pip_diff_edge()`](https://mljaniczek.github.io/multiGGMr/reference/pip_diff_edge.md)
  : Posterior probability of differential edge inclusion (K=2)
- [`pip_edges()`](https://mljaniczek.github.io/multiGGMr/reference/pip_edges.md)
  : Posterior inclusion probabilities (PIP) per group
- [`plot(`*`<multiggm_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/plot.multiggm_fit.md)
  [`plot(`*`<multiggm_fit_list>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/plot.multiggm_fit.md)
  : Plot method for multiggm_fit objects
- [`plot_chain_agreement()`](https://mljaniczek.github.io/multiGGMr/reference/plot_chain_agreement.md)
  : Plot edge probabilities of one chain against another
- [`plot_network()`](https://mljaniczek.github.io/multiGGMr/reference/plot_network.md)
  : Network plot of an adjacency matrix
- [`plot_pip_heatmap()`](https://mljaniczek.github.io/multiGGMr/reference/plot_pip_heatmap.md)
  : Heatmap of posterior inclusion probabilities
- [`plot_ppc()`](https://mljaniczek.github.io/multiGGMr/reference/plot_ppc.md)
  : Plot posterior predictive checks
- [`plot_recovery()`](https://mljaniczek.github.io/multiGGMr/reference/plot_recovery.md)
  : Compare estimated graph to ground truth
- [`plot_trace()`](https://mljaniczek.github.io/multiGGMr/reference/plot_trace.md)
  : Trace plots for MCMC diagnostics
- [`posterior_ci()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_ci.md)
  : Posterior credible intervals for partial correlations
- [`posterior_pcor()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_pcor.md)
  : Compute posterior partial correlations from a multiggm_fit object
- [`posterior_precision()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_precision.md)
  : Extract posterior summary of precision matrices
- [`posterior_predictive_check()`](https://mljaniczek.github.io/multiGGMr/reference/posterior_predictive_check.md)
  : Posterior predictive checks
- [`print(`*`<multiggm_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/print.multiggm_fit.md)
  : Print a multiggm_fit object
- [`print(`*`<multiggm_fit_list>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/print.multiggm_fit_list.md)
  : Print a multiggm_fit_list object
- [`print(`*`<multiggm_platform_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/print.multiggm_platform_fit.md)
  : Print a multiggm_platform_fit object
- [`select_diff_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_diff_edges.md)
  : Select differential edges with Bayesian FDR control
- [`select_edges()`](https://mljaniczek.github.io/multiGGMr/reference/select_edges.md)
  : Select edges within one group by posterior probability
- [`simulate_multiggm()`](https://mljaniczek.github.io/multiGGMr/reference/simulate_multiggm.md)
  : Simulate data from multiple Gaussian graphical models
- [`simulate_multiggm_platform()`](https://mljaniczek.github.io/multiGGMr/reference/simulate_multiggm_platform.md)
  : Simulate multi-platform data from multiple Gaussian graphical models
- [`summary(`*`<multiggm_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/summary.multiggm_fit.md)
  : Summarize a multiggm_fit object
- [`summary(`*`<multiggm_fit_list>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/summary.multiggm_fit_list.md)
  : Summarize a multi-chain multiggm fit
- [`summary(`*`<multiggm_platform_fit>`*`)`](https://mljaniczek.github.io/multiGGMr/reference/summary.multiggm_platform_fit.md)
  : Summarize a multiggm_platform_fit object
