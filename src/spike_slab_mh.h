// Metropolis-Hastings pieces for a similarity parameter x with prior
//   x = 0 with probability 1 - w,  x ~ Gamma(shape, scale) with probability w,
// updated by a between-model move (spike <-> slab) and a within-slab move,
// both proposing from an independent Gamma(shape_prop, scale_prop).
// All gamma parameters are SCALES, matching R::rgamma and R::dgamma; callers
// convert rates with 1 / rate. The log ratios exclude the likelihood terms.

#ifndef MULTIGGM_SPIKE_SLAB_MH_H
#define MULTIGGM_SPIKE_SLAB_MH_H

#include <Rcpp.h>
#include <cmath>

struct GammaSpikeSlab {
  double shape, scale, shape_prop, scale_prop, w;
};

inline double spike_slab_propose(const GammaSpikeSlab& g) {
  return R::rgamma(g.shape_prop, g.scale_prop);
}

inline double spike_slab_log_prior(double x, const GammaSpikeSlab& g) {
  if (x == 0.0) return std::log(1.0 - g.w);
  return std::log(g.w) + R::dgamma(x, g.shape, g.scale, 1);
}

// Spike -> slab (x_curr == 0) or slab -> spike (x_prop == 0).
inline double spike_slab_log_ratio_between(double x_curr, double x_prop,
                                           const GammaSpikeSlab& g) {
  double log_q_fwd = (x_curr == 0.0) ? R::dgamma(x_prop, g.shape_prop, g.scale_prop, 1) : 0.0;
  double log_q_rev = (x_prop == 0.0) ? R::dgamma(x_curr, g.shape_prop, g.scale_prop, 1) : 0.0;
  return spike_slab_log_prior(x_prop, g) - spike_slab_log_prior(x_curr, g)
       + log_q_rev - log_q_fwd;
}

// Slab -> slab, both values > 0.
inline double spike_slab_log_ratio_within(double x_curr, double x_prop,
                                          const GammaSpikeSlab& g) {
  return R::dgamma(x_prop, g.shape, g.scale, 1) - R::dgamma(x_curr, g.shape, g.scale, 1)
       + R::dgamma(x_curr, g.shape_prop, g.scale_prop, 1)
       - R::dgamma(x_prop, g.shape_prop, g.scale_prop, 1);
}

#endif
