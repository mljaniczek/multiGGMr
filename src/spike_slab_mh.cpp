#include <Rcpp.h>
#include "spike_slab_mh.h"

// Test hook: runs the between- and within-model moves with no likelihood, so
// the chain should sample the spike-and-slab prior itself.
// [[Rcpp::export]]
Rcpp::NumericVector spike_slab_prior_chain_cpp(int n_iter, double shape, double scale,
                                               double shape_prop, double scale_prop,
                                               double w) {
  GammaSpikeSlab g = {shape, scale, shape_prop, scale_prop, w};
  Rcpp::NumericVector out(n_iter);
  double x = 0.0;
  for (int it = 0; it < n_iter; it++) {
    double x_prop = (x == 0.0) ? spike_slab_propose(g) : 0.0;
    double log_ar = spike_slab_log_ratio_between(x, x_prop, g);
    if (std::isfinite(log_ar) && log_ar > std::log(R::runif(0.0, 1.0))) x = x_prop;
    if (x != 0.0) {
      double x_prop2 = spike_slab_propose(g);
      double log_ar2 = spike_slab_log_ratio_within(x, x_prop2, g);
      if (std::isfinite(log_ar2) && log_ar2 > std::log(R::runif(0.0, 1.0))) x = x_prop2;
    }
    out[it] = x;
  }
  return out;
}

// [[Rcpp::export]]
double spike_slab_log_ratio_cpp(double x_curr, double x_prop, double shape, double scale,
                                double shape_prop, double scale_prop, double w) {
  GammaSpikeSlab g = {shape, scale, shape_prop, scale_prop, w};
  if (x_curr != 0.0 && x_prop != 0.0) return spike_slab_log_ratio_within(x_curr, x_prop, g);
  return spike_slab_log_ratio_between(x_curr, x_prop, g);
}
