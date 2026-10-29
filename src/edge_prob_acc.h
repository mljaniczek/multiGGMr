// Running sums for per-edge posterior probabilities, updated once per
// post-burn-in iteration by every sampler. Nothing here draws random numbers,
// so adding the accumulator leaves each sampler's chain unchanged.
//
// Edges are indexed e = 0, ..., p(p-1)/2 - 1 in R's upper.tri() order
// (column j = 1..p-1 outer, row i < j inner). Group pairs (h1 < h2) are in
// the order (0,1), (0,2), ..., (1,2), ... Counts and sums (not means) are
// stored so chains pool exactly by adding them.
//
// Per group k, for scale s (0 = precision entry, 1 = partial correlation):
//   sum, sum of squares, and count of |x| > c for each c in deltas;
//   count of the structural indicator (edge indicator G, or for the
//   horseshoe S = 1(w > 0.5) with w = lambda^2 / (1 + lambda^2));
//   horseshoe only: sum of w and count of the cut-model indicator Z.
// Per pair (k, m), for D = x_k - x_m on each scale:
//   sum, sum of squares, counts of D > 0 and D < 0, of |D| > c, and of
//   (structural difference OR |D| > c), the "union" event; counts of
//   structural differences (and, for the horseshoe, cut-model differences,
//   which are also the structural part of the union).

#ifndef MULTIGGM_EDGE_PROB_ACC_H
#define MULTIGGM_EDGE_PROB_ACC_H

#include <Rcpp.h>
#include <vector>
#include <cmath>

class EdgeProbAcc {
public:
  int p, K, npe, ngp, nd;
  bool shrink;
  double n;
  std::vector<double> deltas;
  std::vector<int> gp1, gp2;

  std::vector<double> g_sum[2], g_sumsq[2], g_exceed[2];
  std::vector<double> g_struct, g_wsum, g_cut, g_diag;
  std::vector<double> d_sum[2], d_sumsq[2], d_pos[2], d_neg[2], d_exceed[2], d_union[2];
  std::vector<double> d_struct, d_cut;

  EdgeProbAcc(int p_, int K_, const std::vector<double>& deltas_, bool shrink_)
    : p(p_), K(K_), npe(p_ * (p_ - 1) / 2), ngp(K_ * (K_ - 1) / 2),
      nd((int)deltas_.size()), shrink(shrink_), n(0.0), deltas(deltas_) {
    for (int h = 0; h < K; h++)
      for (int m = h + 1; m < K; m++) { gp1.push_back(h); gp2.push_back(m); }
    for (int s = 0; s < 2; s++) {
      g_sum[s].assign(npe * K, 0.0);
      g_sumsq[s].assign(npe * K, 0.0);
      g_exceed[s].assign(npe * nd * K, 0.0);
      d_sum[s].assign(npe * ngp, 0.0);
      d_sumsq[s].assign(npe * ngp, 0.0);
      d_pos[s].assign(npe * ngp, 0.0);
      d_neg[s].assign(npe * ngp, 0.0);
      d_exceed[s].assign(npe * nd * ngp, 0.0);
      d_union[s].assign(npe * nd * ngp, 0.0);
      xs[s].assign(npe * K, 0.0);
    }
    g_struct.assign(npe * K, 0.0);
    g_diag.assign(p * K, 0.0);
    d_struct.assign(npe * ngp, 0.0);
    if (shrink) {
      g_wsum.assign(npe * K, 0.0);
      g_cut.assign(npe * K, 0.0);
      d_cut.assign(npe * ngp, 0.0);
    }
    st.assign(npe * K, 0);
    ct.assign(npe * K, 0);
    dg.assign(p * K, 0.0);
  }

  // omega(i, j, k): precision entry. structural(i, j, k): 0/1 edge indicator
  // (ignored when shrink). weight(i, j, k): w = lambda^2 / (1 + lambda^2) and
  // cut(i, j, k): 0/1 cut-model indicator (used only when shrink).
  template <class FO, class FS, class FW, class FC>
  void add(FO omega, FS structural, FW weight, FC cut) {
    n += 1.0;
    for (int k = 0; k < K; k++)
      for (int j = 0; j < p; j++) {
        double v = omega(j, j, k);
        g_diag[j + p * k] += v;
        dg[j + p * k] = std::sqrt(v);
      }

    for (int k = 0; k < K; k++) {
      int e = 0;
      for (int j = 1; j < p; j++) {
        for (int i = 0; i < j; i++) {
          int idx = e + npe * k;
          double w = omega(i, j, k);
          xs[0][idx] = w;
          xs[1][idx] = -w / (dg[i + p * k] * dg[j + p * k]);
          unsigned char s;
          if (shrink) {
            double ww = weight(i, j, k);
            g_wsum[idx] += ww;
            s = (ww > 0.5) ? 1 : 0;
            unsigned char c = cut(i, j, k) ? 1 : 0;
            ct[idx] = c;
            g_cut[idx] += c;
          } else {
            s = structural(i, j, k) ? 1 : 0;
          }
          st[idx] = s;
          g_struct[idx] += s;
          for (int sc = 0; sc < 2; sc++) {
            double x = xs[sc][idx];
            g_sum[sc][idx] += x;
            g_sumsq[sc][idx] += x * x;
            double ax = std::abs(x);
            for (int t = 0; t < nd; t++)
              if (ax > deltas[t]) g_exceed[sc][e + npe * (t + nd * k)] += 1.0;
          }
          e++;
        }
      }
    }

    for (int g = 0; g < ngp; g++) {
      int h1 = gp1[g], h2 = gp2[g];
      for (int e = 0; e < npe; e++) {
        int i1 = e + npe * h1, i2 = e + npe * h2, ig = e + npe * g;
        bool sdiff = st[i1] != st[i2];
        d_struct[ig] += sdiff;
        bool udiff = sdiff;
        if (shrink) {
          bool cdiff = ct[i1] != ct[i2];
          d_cut[ig] += cdiff;
          udiff = cdiff;
        }
        for (int sc = 0; sc < 2; sc++) {
          double D = xs[sc][i1] - xs[sc][i2];
          d_sum[sc][ig] += D;
          d_sumsq[sc][ig] += D * D;
          if (D > 0.0) d_pos[sc][ig] += 1.0;
          else if (D < 0.0) d_neg[sc][ig] += 1.0;
          double aD = std::abs(D);
          for (int t = 0; t < nd; t++) {
            int it = e + npe * (t + nd * g);
            bool ex = aD > deltas[t];
            if (ex) d_exceed[sc][it] += 1.0;
            if (ex || udiff) d_union[sc][it] += 1.0;
          }
        }
      }
    }
  }

  Rcpp::List to_list() const {
    Rcpp::IntegerMatrix pairs(ngp, 2);
    for (int g = 0; g < ngp; g++) { pairs(g, 0) = gp1[g] + 1; pairs(g, 1) = gp2[g] + 1; }
    const char* scale_names[2] = {"precision", "pcor"};
    Rcpp::List grp, pr;
    for (int sc = 0; sc < 2; sc++) {
      grp[scale_names[sc]] = Rcpp::List::create(
        Rcpp::_["sum"] = mat(g_sum[sc], npe, K),
        Rcpp::_["sumsq"] = mat(g_sumsq[sc], npe, K),
        Rcpp::_["n_exceed"] = arr3(g_exceed[sc], npe, nd, K));
      pr[scale_names[sc]] = Rcpp::List::create(
        Rcpp::_["sum"] = mat(d_sum[sc], npe, ngp),
        Rcpp::_["sumsq"] = mat(d_sumsq[sc], npe, ngp),
        Rcpp::_["n_pos"] = mat(d_pos[sc], npe, ngp),
        Rcpp::_["n_neg"] = mat(d_neg[sc], npe, ngp),
        Rcpp::_["n_exceed"] = arr3(d_exceed[sc], npe, nd, ngp),
        Rcpp::_["n_union"] = arr3(d_union[sc], npe, nd, ngp));
    }
    grp["n_struct"] = mat(g_struct, npe, K);
    grp["diag_sum"] = mat(g_diag, p, K);
    pr["n_struct"] = mat(d_struct, npe, ngp);
    if (shrink) {
      grp["w_sum"] = mat(g_wsum, npe, K);
      grp["n_cut"] = mat(g_cut, npe, K);
      pr["n_cut"] = mat(d_cut, npe, ngp);
    }
    return Rcpp::List::create(
      Rcpp::_["n"] = n, Rcpp::_["p"] = p, Rcpp::_["K"] = K,
      Rcpp::_["deltas"] = Rcpp::NumericVector(deltas.begin(), deltas.end()),
      Rcpp::_["group_pairs"] = pairs, Rcpp::_["shrinkage"] = shrink,
      Rcpp::_["group"] = grp, Rcpp::_["pair"] = pr);
  }

private:
  std::vector<double> xs[2];
  std::vector<unsigned char> st, ct;
  std::vector<double> dg;

  static Rcpp::NumericVector mat(const std::vector<double>& v, int r, int c) {
    Rcpp::NumericVector out(v.begin(), v.end());
    out.attr("dim") = Rcpp::IntegerVector::create(r, c);
    return out;
  }
  static Rcpp::NumericVector arr3(const std::vector<double>& v, int a, int b, int c) {
    Rcpp::NumericVector out(v.begin(), v.end());
    out.attr("dim") = Rcpp::IntegerVector::create(a, b, c);
    return out;
  }
};

#endif
