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
//   sum, sum of squares, counts of D > 0, D < 0 and D > c, of |D| > c, and
//   of (structural difference OR |D| > c), the "union" event; counts of
//   structural differences (and, for the horseshoe, cut-model differences,
//   which are also the structural part of the union).
// Per pair and node i, counted per iteration: the number of differential
//   edges at i (structural, |D| > c, union; summed over iterations), the
//   number of iterations with at least one, and the sum of |D_ij| over j.
//
// Diagnostics:
//   batch counts: every count used for edge selection is also kept per batch
//     of consecutive iterations, as the sum over batches of count^2 / length,
//     for batch-means Monte Carlo standard errors;
//   traces: per iteration, edges per group, differential edges per pair, and
//     the log-likelihood per group;
//   WAIC sums (when the data are given): per observation, log sum exp of the
//     log-likelihood and its sum and sum of squares.

#ifndef MULTIGGM_EDGE_PROB_ACC_H
#define MULTIGGM_EDGE_PROB_ACC_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

// A count, optionally also kept per batch of iterations.
struct BatchCount {
  std::vector<double> tot, cur, sq;
  void init(int n, bool batched) {
    tot.assign(n, 0.0);
    if (batched) { cur.assign(n, 0.0); sq.assign(n, 0.0); }
  }
  inline void inc(int i) {
    tot[i] += 1.0;
    if (!cur.empty()) cur[i] += 1.0;
  }
  void flush(double len) {
    for (size_t i = 0; i < cur.size(); i++) { sq[i] += cur[i] * cur[i] / len; cur[i] = 0.0; }
  }
};

class EdgeProbAcc {
public:
  int p, K, npe, ngp, nd;
  bool shrink;
  double n;
  std::vector<double> deltas;
  std::vector<int> gp1, gp2;

  std::vector<double> g_sum[2], g_sumsq[2], g_wsum, g_diag;
  BatchCount g_exceed[2], g_struct, g_cut;
  std::vector<double> d_sum[2], d_sumsq[2], d_pos[2], d_neg[2], d_up[2];
  BatchCount d_exceed[2], d_union[2], d_struct, d_cut;

  // Node level, per pair: [p, ngp] or [p, nd, ngp]
  std::vector<double> nd_rewire[2], nd_deg_exceed[2], nd_any_exceed[2],
    nd_deg_union[2], nd_any_union[2], nd_deg_struct, nd_any_struct,
    nd_deg_cut, nd_any_cut;

  // n_total: number of add() calls the chain will make (for batches and
  // traces); n_batch: number of batches (0 = no batch counts).
  EdgeProbAcc(int p_, int K_, const std::vector<double>& deltas_, bool shrink_,
              int n_total_ = 0, int n_batch_ = 0)
    : p(p_), K(K_), npe(p_ * (p_ - 1) / 2), ngp(K_ * (K_ - 1) / 2),
      nd((int)deltas_.size()), shrink(shrink_), n(0.0), deltas(deltas_),
      n_total(n_total_), t(0), has_data(false), has_waic(false) {
    for (int h = 0; h < K; h++)
      for (int m = h + 1; m < K; m++) { gp1.push_back(h); gp2.push_back(m); }
    n_batch = std::min(n_batch_, n_total);
    if (n_batch < 2) n_batch = 0;
    batch_len = n_batch > 0 ? n_total / n_batch : 0;
    bool bt = n_batch > 0;
    for (int s = 0; s < 2; s++) {
      g_sum[s].assign(npe * K, 0.0);
      g_sumsq[s].assign(npe * K, 0.0);
      g_exceed[s].init(npe * nd * K, bt);
      d_sum[s].assign(npe * ngp, 0.0);
      d_sumsq[s].assign(npe * ngp, 0.0);
      d_pos[s].assign(npe * ngp, 0.0);
      d_neg[s].assign(npe * ngp, 0.0);
      d_up[s].assign(npe * nd * ngp, 0.0);
      d_exceed[s].init(npe * nd * ngp, bt);
      d_union[s].init(npe * nd * ngp, bt);
      xs[s].assign(npe * K, 0.0);
      nd_rewire[s].assign(p * ngp, 0.0);
      nd_deg_exceed[s].assign(p * nd * ngp, 0.0);
      nd_any_exceed[s].assign(p * nd * ngp, 0.0);
      nd_deg_union[s].assign(p * nd * ngp, 0.0);
      nd_any_union[s].assign(p * nd * ngp, 0.0);
    }
    g_struct.init(npe * K, bt);
    g_diag.assign(p * K, 0.0);
    d_struct.init(npe * ngp, bt);
    nd_deg_struct.assign(p * ngp, 0.0);
    nd_any_struct.assign(p * ngp, 0.0);
    if (shrink) {
      g_wsum.assign(npe * K, 0.0);
      g_cut.init(npe * K, bt);
      d_cut.init(npe * ngp, bt);
      nd_deg_cut.assign(p * ngp, 0.0);
      nd_any_cut.assign(p * ngp, 0.0);
    }
    st.assign(npe * K, 0);
    ct.assign(npe * K, 0);
    dg.assign(p * K, 0.0);
    ei.resize(npe); ej.resize(npe);
    int e = 0;
    for (int j = 1; j < p; j++)
      for (int i = 0; i < j; i++) { ei[e] = i; ej[e] = j; e++; }
    tr_edges.assign((size_t)n_total * K, 0.0);
    tr_loglik.assign((size_t)n_total * K, NA_REAL);
    tr_struct.assign((size_t)n_total * ngp, 0.0);
    if (shrink) tr_cut.assign((size_t)n_total * ngp, 0.0);
    for (int s = 0; s < 2; s++) {
      tr_exceed[s].assign((size_t)n_total * nd * ngp, 0.0);
      tr_union[s].assign((size_t)n_total * nd * ngp, 0.0);
    }
    deg_s.assign(p, 0); deg_c.assign(p, 0);
    deg_e.assign(p * nd, 0); deg_u.assign(p * nd, 0);
  }

  // Data for the log-likelihood trace: S_k = X_k' X_k (centered) and n_k.
  // X (optional, n_k x p per group) adds the WAIC sums.
  void set_data(const Rcpp::List& S, const Rcpp::NumericVector& nv,
                Rcpp::Nullable<Rcpp::List> X = R_NilValue) {
    S_dat.assign((size_t)p * p * K, 0.0);
    n_dat.assign(nv.begin(), nv.end());
    for (int k = 0; k < K; k++) {
      Rcpp::NumericMatrix Sk(S[k]);
      std::copy(Sk.begin(), Sk.end(), S_dat.begin() + (size_t)p * p * k);
    }
    has_data = true;
    if (X.isNotNull()) {
      Rcpp::List Xl(X);
      Xt.resize(K); w_lse.resize(K); w_sum.resize(K); w_sumsq.resize(K); nobs.resize(K);
      for (int k = 0; k < K; k++) {
        Rcpp::NumericMatrix Xk(Xl[k]);
        int nk = Xk.nrow();
        nobs[k] = nk;
        Xt[k].assign((size_t)nk * p, 0.0);  // observation-major
        for (int r = 0; r < nk; r++)
          for (int c = 0; c < p; c++) Xt[k][(size_t)r * p + c] = Xk(r, c);
        w_lse[k].assign(nk, R_NegInf);
        w_sum[k].assign(nk, 0.0);
        w_sumsq[k].assign(nk, 0.0);
      }
      has_waic = true;
      ybuf.assign(p, 0.0);
    }
  }

  // omega(i, j, k): precision entry. structural(i, j, k): 0/1 edge indicator
  // (ignored when shrink). weight(i, j, k): w = lambda^2 / (1 + lambda^2) and
  // cut(i, j, k): 0/1 cut-model indicator (used only when shrink).
  template <class FO, class FS, class FW, class FC>
  void add(FO omega, FS structural, FW weight, FC cut) {
    n += 1.0;
    const bool trace = t < n_total;
    for (int k = 0; k < K; k++)
      for (int j = 0; j < p; j++) {
        double v = omega(j, j, k);
        g_diag[j + p * k] += v;
        dg[j + p * k] = std::sqrt(v);
      }

    for (int k = 0; k < K; k++) {
      int e = 0;
      double n_edges = 0.0;
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
            if (c) g_cut.inc(idx);
          } else {
            s = structural(i, j, k) ? 1 : 0;
          }
          st[idx] = s;
          if (s) { g_struct.inc(idx); n_edges += 1.0; }
          for (int sc = 0; sc < 2; sc++) {
            double x = xs[sc][idx];
            g_sum[sc][idx] += x;
            g_sumsq[sc][idx] += x * x;
            double ax = std::abs(x);
            for (int tt = 0; tt < nd; tt++)
              if (ax > deltas[tt]) g_exceed[sc].inc(e + npe * (tt + nd * k));
          }
          e++;
        }
      }
      if (trace) {
        tr_edges[t + (size_t)n_total * k] = n_edges;
        if (has_data) tr_loglik[t + (size_t)n_total * k] = loglik(omega, k);
      }
    }

    for (int g = 0; g < ngp; g++) {
      int h1 = gp1[g], h2 = gp2[g];
      double n_s = 0.0, n_c = 0.0;
      std::fill(deg_s.begin(), deg_s.end(), 0);
      std::fill(deg_c.begin(), deg_c.end(), 0);
      for (int sc = 0; sc < 2; sc++) {
        std::fill(deg_e.begin(), deg_e.end(), 0);
        std::fill(deg_u.begin(), deg_u.end(), 0);
        std::vector<double> n_e(nd, 0.0), n_u(nd, 0.0);
        for (int e = 0; e < npe; e++) {
          int i1 = e + npe * h1, i2 = e + npe * h2, ig = e + npe * g;
          int a = ei[e], b = ej[e];
          bool sdiff = st[i1] != st[i2];
          bool udiff = sdiff;
          if (shrink) udiff = ct[i1] != ct[i2];
          if (sc == 0) {
            if (sdiff) { d_struct.inc(ig); n_s += 1.0; deg_s[a]++; deg_s[b]++; }
            if (shrink && udiff) { d_cut.inc(ig); n_c += 1.0; deg_c[a]++; deg_c[b]++; }
          }
          double D = xs[sc][i1] - xs[sc][i2];
          d_sum[sc][ig] += D;
          d_sumsq[sc][ig] += D * D;
          if (D > 0.0) d_pos[sc][ig] += 1.0;
          else if (D < 0.0) d_neg[sc][ig] += 1.0;
          double aD = std::abs(D);
          nd_rewire[sc][a + p * g] += aD;
          nd_rewire[sc][b + p * g] += aD;
          for (int tt = 0; tt < nd; tt++) {
            int it = e + npe * (tt + nd * g);
            bool ex = aD > deltas[tt];
            if (ex) {
              d_exceed[sc].inc(it); n_e[tt] += 1.0;
              deg_e[a + p * tt]++; deg_e[b + p * tt]++;
              if (D > 0.0) d_up[sc][it] += 1.0;
            }
            if (ex || udiff) {
              d_union[sc].inc(it); n_u[tt] += 1.0;
              deg_u[a + p * tt]++; deg_u[b + p * tt]++;
            }
          }
        }
        for (int tt = 0; tt < nd; tt++) {
          for (int i = 0; i < p; i++) {
            int in = i + p * (tt + nd * g);
            nd_deg_exceed[sc][in] += deg_e[i + p * tt];
            if (deg_e[i + p * tt] > 0) nd_any_exceed[sc][in] += 1.0;
            nd_deg_union[sc][in] += deg_u[i + p * tt];
            if (deg_u[i + p * tt] > 0) nd_any_union[sc][in] += 1.0;
          }
          if (trace) {
            size_t ti = t + (size_t)n_total * (tt + nd * g);
            tr_exceed[sc][ti] = n_e[tt];
            tr_union[sc][ti] = n_u[tt];
          }
        }
      }
      for (int i = 0; i < p; i++) {
        nd_deg_struct[i + p * g] += deg_s[i];
        if (deg_s[i] > 0) nd_any_struct[i + p * g] += 1.0;
        if (shrink) {
          nd_deg_cut[i + p * g] += deg_c[i];
          if (deg_c[i] > 0) nd_any_cut[i + p * g] += 1.0;
        }
      }
      if (trace) {
        tr_struct[t + (size_t)n_total * g] = n_s;
        if (shrink) tr_cut[t + (size_t)n_total * g] = n_c;
      }
    }

    if (n_batch > 0) {
      int b = std::min((int)(t / batch_len), n_batch - 1);
      int end = (b == n_batch - 1) ? n_total : (b + 1) * batch_len;
      if (t + 1 == end) {
        double len = end - b * batch_len;
        for (int s = 0; s < 2; s++) {
          g_exceed[s].flush(len); d_exceed[s].flush(len); d_union[s].flush(len);
        }
        g_struct.flush(len); d_struct.flush(len);
        if (shrink) { g_cut.flush(len); d_cut.flush(len); }
      }
    }
    t++;
  }

  Rcpp::List to_list() const {
    Rcpp::IntegerMatrix pairs(ngp, 2);
    for (int g = 0; g < ngp; g++) { pairs(g, 0) = gp1[g] + 1; pairs(g, 1) = gp2[g] + 1; }
    const char* scale_names[2] = {"precision", "pcor"};
    Rcpp::List grp, pr, node, bg, bp;
    for (int sc = 0; sc < 2; sc++) {
      grp[scale_names[sc]] = Rcpp::List::create(
        Rcpp::_["sum"] = mat(g_sum[sc], npe, K),
        Rcpp::_["sumsq"] = mat(g_sumsq[sc], npe, K),
        Rcpp::_["n_exceed"] = arr3(g_exceed[sc].tot, npe, nd, K));
      pr[scale_names[sc]] = Rcpp::List::create(
        Rcpp::_["sum"] = mat(d_sum[sc], npe, ngp),
        Rcpp::_["sumsq"] = mat(d_sumsq[sc], npe, ngp),
        Rcpp::_["n_pos"] = mat(d_pos[sc], npe, ngp),
        Rcpp::_["n_neg"] = mat(d_neg[sc], npe, ngp),
        Rcpp::_["n_exceed"] = arr3(d_exceed[sc].tot, npe, nd, ngp),
        Rcpp::_["n_union"] = arr3(d_union[sc].tot, npe, nd, ngp),
        Rcpp::_["n_up"] = arr3(d_up[sc], npe, nd, ngp));
      node[scale_names[sc]] = Rcpp::List::create(
        Rcpp::_["rewire_sum"] = mat(nd_rewire[sc], p, ngp),
        Rcpp::_["deg_exceed"] = arr3(nd_deg_exceed[sc], p, nd, ngp),
        Rcpp::_["n_any_exceed"] = arr3(nd_any_exceed[sc], p, nd, ngp),
        Rcpp::_["deg_union"] = arr3(nd_deg_union[sc], p, nd, ngp),
        Rcpp::_["n_any_union"] = arr3(nd_any_union[sc], p, nd, ngp));
      if (n_batch > 0) {
        bg[scale_names[sc]] = Rcpp::List::create(
          Rcpp::_["n_exceed"] = arr3(g_exceed[sc].sq, npe, nd, K));
        bp[scale_names[sc]] = Rcpp::List::create(
          Rcpp::_["n_exceed"] = arr3(d_exceed[sc].sq, npe, nd, ngp),
          Rcpp::_["n_union"] = arr3(d_union[sc].sq, npe, nd, ngp));
      }
    }
    grp["n_struct"] = mat(g_struct.tot, npe, K);
    grp["diag_sum"] = mat(g_diag, p, K);
    pr["n_struct"] = mat(d_struct.tot, npe, ngp);
    node["deg_struct"] = mat(nd_deg_struct, p, ngp);
    node["n_any_struct"] = mat(nd_any_struct, p, ngp);
    if (n_batch > 0) {
      bg["n_struct"] = mat(g_struct.sq, npe, K);
      bp["n_struct"] = mat(d_struct.sq, npe, ngp);
    }
    if (shrink) {
      grp["w_sum"] = mat(g_wsum, npe, K);
      grp["n_cut"] = mat(g_cut.tot, npe, K);
      pr["n_cut"] = mat(d_cut.tot, npe, ngp);
      node["deg_cut"] = mat(nd_deg_cut, p, ngp);
      node["n_any_cut"] = mat(nd_any_cut, p, ngp);
      if (n_batch > 0) {
        bg["n_cut"] = mat(g_cut.sq, npe, K);
        bp["n_cut"] = mat(d_cut.sq, npe, ngp);
      }
    }
    Rcpp::RObject batch = R_NilValue;
    if (n_batch > 0) {
      Rcpp::NumericVector len(n_batch, (double)batch_len);
      len[n_batch - 1] = n_total - (n_batch - 1) * batch_len;
      batch = Rcpp::List::create(Rcpp::_["length"] = len,
                                 Rcpp::_["group"] = bg, Rcpp::_["pair"] = bp);
    }
    return Rcpp::List::create(
      Rcpp::_["n"] = n, Rcpp::_["p"] = p, Rcpp::_["K"] = K,
      Rcpp::_["deltas"] = Rcpp::NumericVector(deltas.begin(), deltas.end()),
      Rcpp::_["group_pairs"] = pairs, Rcpp::_["shrinkage"] = shrink,
      Rcpp::_["group"] = grp, Rcpp::_["pair"] = pr, Rcpp::_["node"] = node,
      Rcpp::_["batch"] = batch);
  }

  // Per-iteration traces, [n_iter, ...].
  Rcpp::List trace_list() const {
    int nt = std::min(t, n_total);
    const char* scale_names[2] = {"precision", "pcor"};
    Rcpp::List ex, un;
    for (int sc = 0; sc < 2; sc++) {
      ex[scale_names[sc]] = arr3(head(tr_exceed[sc], nt, nd * ngp), nt, nd, ngp);
      un[scale_names[sc]] = arr3(head(tr_union[sc], nt, nd * ngp), nt, nd, ngp);
    }
    Rcpp::List out = Rcpp::List::create(
      Rcpp::_["n_edges"] = mat(head(tr_edges, nt, K), nt, K),
      Rcpp::_["loglik"] = has_data ? Rcpp::RObject(mat(head(tr_loglik, nt, K), nt, K))
                                   : Rcpp::RObject(R_NilValue),
      Rcpp::_["n_diff_struct"] = mat(head(tr_struct, nt, ngp), nt, ngp),
      Rcpp::_["n_diff_exceed"] = ex,
      Rcpp::_["n_diff_union"] = un);
    if (shrink) out["n_diff_cut"] = mat(head(tr_cut, nt, ngp), nt, ngp);
    return out;
  }

  // WAIC sums per group and observation, or NULL.
  Rcpp::RObject waic_list() const {
    if (!has_waic) return R_NilValue;
    Rcpp::List out(K);
    for (int k = 0; k < K; k++)
      out[k] = Rcpp::List::create(
        Rcpp::_["n"] = n,
        Rcpp::_["lse"] = Rcpp::NumericVector(w_lse[k].begin(), w_lse[k].end()),
        Rcpp::_["sum"] = Rcpp::NumericVector(w_sum[k].begin(), w_sum[k].end()),
        Rcpp::_["sumsq"] = Rcpp::NumericVector(w_sumsq[k].begin(), w_sumsq[k].end()));
    return out;
  }

private:
  std::vector<double> xs[2];
  std::vector<unsigned char> st, ct;
  std::vector<double> dg;
  std::vector<int> ei, ej, deg_s, deg_c, deg_e, deg_u;
  int n_total, t, n_batch, batch_len;
  std::vector<double> tr_edges, tr_loglik, tr_struct, tr_cut, tr_exceed[2], tr_union[2];
  bool has_data, has_waic;
  std::vector<double> S_dat, n_dat, ybuf, Lbuf;
  std::vector<int> nobs;
  std::vector<std::vector<double> > Xt, w_lse, w_sum, w_sumsq;

  // Gaussian log-likelihood of group k at precision omega(., ., k), from
  // S_k and n_k; also adds the per-observation terms to the WAIC sums.
  template <class FO>
  double loglik(FO omega, int k) {
    Lbuf.assign((size_t)p * p, 0.0);
    double trSO = 0.0;
    const double* Sk = &S_dat[(size_t)p * p * k];
    for (int j = 0; j < p; j++)
      for (int i = 0; i < p; i++) {
        double o = omega(i, j, k);
        trSO += Sk[i + (size_t)p * j] * o;
        if (i >= j) Lbuf[i + (size_t)p * j] = o;
      }
    // In-place lower Cholesky factor, Omega = L L'
    for (int j = 0; j < p; j++) {
      double* Lj = &Lbuf[(size_t)p * j];
      for (int m = 0; m < j; m++) {
        const double* Lm = &Lbuf[(size_t)p * m];
        double ljm = Lm[j];
        if (ljm != 0.0)
          for (int i = j; i < p; i++) Lj[i] -= Lm[i] * ljm;
      }
      if (!(Lj[j] > 0.0)) return NA_REAL;
      double d = std::sqrt(Lj[j]);
      for (int i = j; i < p; i++) Lj[i] /= d;
    }
    double logdet = 0.0;
    for (int j = 0; j < p; j++) logdet += 2.0 * std::log(Lbuf[j + (size_t)p * j]);
    const double c0 = 0.5 * logdet - 0.5 * p * std::log(2.0 * M_PI);
    if (has_waic) {
      for (int r = 0; r < nobs[k]; r++) {
        const double* y = &Xt[k][(size_t)r * p];
        double q = 0.0;
        for (int j = 0; j < p; j++) {
          const double* Lj = &Lbuf[(size_t)p * j];
          double z = 0.0;
          for (int i = j; i < p; i++) z += Lj[i] * y[i];
          q += z * z;
        }
        double l = c0 - 0.5 * q;
        double& a = w_lse[k][r];
        a = (a == R_NegInf) ? l : std::max(a, l) + std::log1p(std::exp(-std::abs(a - l)));
        w_sum[k][r] += l;
        w_sumsq[k][r] += l * l;
      }
    }
    return n_dat[k] * c0 - 0.5 * trSO;
  }

  static std::vector<double> head(const std::vector<double>& v, int nt, int ncol) {
    // v is [n_total, ncol]; keep the first nt rows
    if (ncol == 0 || v.empty()) return std::vector<double>();
    int n_total = (int)(v.size() / ncol);
    if (nt == n_total) return v;
    std::vector<double> out((size_t)nt * ncol);
    for (int c = 0; c < ncol; c++)
      std::copy(v.begin() + (size_t)n_total * c, v.begin() + (size_t)n_total * c + nt,
                out.begin() + (size_t)nt * c);
    return out;
  }
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
