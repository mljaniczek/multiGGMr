#ifndef WANGLI_INTERNAL_H
#define WANGLI_INTERNAL_H

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <map>

// [[Rcpp::depends(RcppArmadillo)]]

// ============================================================================
// Result struct for internal G-Wishart functions (replaces Rcpp::List)
// ============================================================================
struct GWishartResult {
  arma::mat C;
  arma::mat Sig;
  arma::umat adj;
  bool ok;
};

// ============================================================================
// Numeric utilities
// ============================================================================

inline double logsumexp_internal(const arma::vec& x) {
  double xmax = x.max();
  return xmax + std::log(arma::sum(arma::exp(x - xmax)));
}

inline double log_multi_gamma(int p, double a) {
  double out = (double)p * (p - 1) * 0.25 * std::log(M_PI);
  for (int j = 1; j <= p; j++) {
    out += ::lgamma(a + (1.0 - (double)j) / 2.0);
  }
  return out;
}

inline double log_gwishart_complete_const(double b, const arma::mat& D) {
  int p = (int)D.n_rows;
  double nu = b + p - 1.0;
  double sign = 0.0, logdetD = 0.0;
  arma::log_det(logdetD, sign, D);
  return (nu * p / 2.0) * std::log(2.0) - (nu / 2.0) * logdetD + log_multi_gamma(p, nu / 2.0);
}

inline double log_gwishart_complete_pdf(const arma::mat& K, double b, const arma::mat& D) {
  double sign = 0.0, logdetK = 0.0;
  arma::log_det(logdetK, sign, K);
  double trDK = arma::trace(D * K);
  double logZ = log_gwishart_complete_const(b, D);
  return ((b - 2.0) / 2.0) * logdetK - 0.5 * trDK - logZ;
}

// ============================================================================
// Safe linear algebra helpers
// ============================================================================

inline bool safe_chol_try(arma::mat& L, const arma::mat& Ain,
                          const char* which = "lower",
                          double jitter0 = 1e-10, int max_tries = 10) {
  arma::mat A = 0.5 * (Ain + Ain.t());
  bool ok = arma::chol(L, A, which);
  if (ok) return true;
  double eps = jitter0;
  for (int t = 0; t < max_tries; ++t) {
    arma::mat Aj = A + eps * arma::eye(A.n_rows, A.n_cols);
    ok = arma::chol(L, Aj, which);
    if (ok) return true;
    eps *= 10.0;
  }
  return false;
}

inline bool safe_inv_sympd_try(arma::mat& Ainv, const arma::mat& A_in,
                               double jitter0 = 1e-10, int max_tries = 10) {
  arma::mat A = 0.5 * (A_in + A_in.t());
  if (A.n_rows == 1u && A.n_cols == 1u) {
    double a = A(0, 0);
    if (!(a > 0.0)) a += jitter0;
    if (!(a > 0.0)) return false;
    Ainv.set_size(1u, 1u);
    Ainv(0, 0) = 1.0 / a;
    return true;
  }
  bool ok = arma::inv_sympd(Ainv, A);
  if (ok) return true;
  double jitter = jitter0;
  for (int t = 0; t < max_tries; ++t) {
    arma::mat Aj = A + jitter * arma::eye(A.n_rows, A.n_cols);
    ok = arma::inv_sympd(Ainv, Aj);
    if (ok) return true;
    jitter *= 10.0;
  }
  return false;
}

// General symmetric matrix inverse (does NOT require positive-definiteness).
// Used for Woodbury-update intermediates (C_old - K_new) and (Delta - Sig_bb)
// which are differences of PD matrices and can be indefinite.
inline bool safe_inv_general_try(arma::mat& Ainv, const arma::mat& A_in,
                                  double jitter0 = 1e-10, int max_tries = 10) {
  arma::mat A = 0.5 * (A_in + A_in.t());
  if (A.n_rows == 1u && A.n_cols == 1u) {
    double a = A(0, 0);
    if (std::abs(a) < 1e-15) return false;
    Ainv.set_size(1u, 1u);
    Ainv(0, 0) = 1.0 / a;
    return true;
  }
  bool ok = arma::inv(Ainv, A);
  if (ok) return true;
  // Try with small jitter if singular
  double jitter = jitter0;
  for (int t = 0; t < max_tries; ++t) {
    arma::mat Aj = A + jitter * arma::eye(A.n_rows, A.n_cols);
    ok = arma::inv(Ainv, Aj);
    if (ok) return true;
    jitter *= 10.0;
  }
  return false;
}

inline arma::mat safe_inv_sympd(const arma::mat& A_in,
                                double jitter0 = 1e-10, int max_tries = 10) {
  arma::mat Ainv;
  bool ok = safe_inv_sympd_try(Ainv, A_in, jitter0, max_tries);
  if (!ok) Rcpp::stop("safe_inv_sympd: matrix not SPD even after diagonal jitter");
  return Ainv;
}

inline bool wishart_rnd_try(arma::mat& W_out, const arma::mat& Sigma, double df,
                            double jitter0 = 1e-10, int max_tries = 10) {
  int p = (int)Sigma.n_rows;
  arma::mat L;
  if (!safe_chol_try(L, Sigma, "lower", jitter0, max_tries)) return false;
  arma::mat A(p, p, arma::fill::zeros);
  for (int i = 0; i < p; i++) {
    A(i, i) = std::sqrt(R::rchisq(df - i));
    for (int j = 0; j < i; j++) {
      A(i, j) = R::rnorm(0.0, 1.0);
    }
  }
  arma::mat LA = L * A;
  W_out = LA * LA.t();
  W_out = 0.5 * (W_out + W_out.t());
  return true;
}

// ============================================================================
// Internal log_iwishart_invA_const (0-based, no Rcpp export)
// ============================================================================

inline double log_iwishart_invA_const_internal(double df, const arma::mat& S) {
  int p = (int)S.n_rows;
  double sign = 0.0, logdetS = 0.0;
  arma::log_det(logdetS, sign, S);
  return (df + p - 1.0) / 2.0 * (logdetS - p * std::log(2.0))
         - log_multi_gamma(p, (df + p - 1.0) / 2.0);
}

// ============================================================================
// Internal log_J (0-based)
// ============================================================================

inline double log_J_internal(double h, const arma::mat& B, double a11) {
  double B22 = B(1, 1);
  double B11 = B(0, 0);
  double B12 = B(0, 1);
  arma::mat S(1, 1);
  S(0, 0) = B22;
  return 0.5 * std::log(2.0 * M_PI / B22)
         - log_iwishart_invA_const_internal(h, S)
         + (h - 1.0) / 2.0 * std::log(a11)
         - 0.5 * (B11 - (B12 * B12) / B22) * a11;
}

// ============================================================================
// Internal log_H (0-based indices ii, jj)
// ============================================================================

inline double compute_log_H(double b_prior, const arma::mat& D_prior,
                            double n, const arma::mat& S,
                            const arma::mat& C, int ii, int jj) {
  // (i,j)=0 case
  arma::mat C0 = C;
  C0(ii, jj) = 0;
  C0(jj, ii) = 0;
  int e = jj;
  arma::rowvec C_12 = C0.row(e);
  C_12.shed_col(e);
  arma::mat C_22 = C0;
  C_22.shed_row(e);
  C_22.shed_col(e);
  arma::mat invC_22 = safe_inv_sympd(C_22);
  double c = arma::as_scalar(C_12 * invC_22 * C_12.t());
  arma::mat C0_ij(2, 2, arma::fill::zeros);
  C0_ij(0, 0) = C(ii, ii);
  C0_ij(1, 1) = c;

  // (i,j)=1 case
  arma::uvec evec(2);
  evec(0) = ii;
  evec(1) = jj;
  arma::mat C_12b = C.rows(evec);
  C_12b.shed_cols(evec);
  arma::mat C_22b = C;
  arma::uvec rem = arma::sort(evec, "descend");
  C_22b.shed_row(rem(0));
  C_22b.shed_col(rem(0));
  C_22b.shed_row(rem(1));
  C_22b.shed_col(rem(1));
  arma::mat invC_22b = safe_inv_sympd(C_22b);
  arma::mat Ce = C_12b * invC_22b * C_12b.t();
  arma::mat A = C.submat(evec, evec) - Ce;
  double a11 = A(0, 0);
  arma::mat C1_ij = Ce;

  double b_post = b_prior + n;
  arma::mat D_post = D_prior + S;

  arma::mat Ssc(1, 1);
  Ssc(0, 0) = D_post(jj, jj);
  double term1 = -log_iwishart_invA_const_internal(b_post, Ssc);
  double term2 = -log_J_internal(b_post, D_post.submat(evec, evec), a11);
  double term3 = (n + b_prior - 2.0) / 2.0 * std::log(a11);
  double trterm = arma::trace(D_post.submat(evec, evec) * (C0_ij - C1_ij));
  double term4 = -0.5 * trterm;
  return term1 + term2 + term3 + term4;
}

// ============================================================================
// Internal log_GWishart_NOij_pdf (0-based indices ii, jj)
// ============================================================================

inline double compute_log_GWishart_NOij_pdf(double b_prior, const arma::mat& D_prior,
                                            arma::mat C, int ii, int jj, int edgeij) {
  if (edgeij == 0) {
    C(ii, jj) = 0;
    C(jj, ii) = 0;
    arma::rowvec C_12 = C.row(jj);
    C_12.shed_col(jj);
    arma::mat C_22 = C;
    C_22.shed_row(jj);
    C_22.shed_col(jj);
    arma::mat invC_22 = safe_inv_sympd(C_22);
    double c = arma::as_scalar(C_12 * invC_22 * C_12.t());
    arma::mat C0_ij(2, 2, arma::fill::zeros);
    C0_ij(0, 0) = C(ii, ii);
    C0_ij(1, 1) = c;

    double sign = 0.0, logdetC = 0.0;
    arma::log_det(logdetC, sign, C);
    double logJoint = (b_prior - 2.0) / 2.0 * logdetC - 0.5 * arma::trace(D_prior * C);

    arma::uvec evec(2);
    evec(0) = ii;
    evec(1) = jj;
    arma::mat D2 = D_prior.submat(evec, evec);
    double logK2 = log_gwishart_complete_pdf(C0_ij, b_prior, D2);

    arma::mat V = safe_inv_sympd(D2);
    double D_priorii = 1.0 / V(0, 0);
    arma::mat K1(1, 1);
    K1(0, 0) = C0_ij(0, 0);
    arma::mat D1(1, 1);
    D1(0, 0) = D_priorii;
    double logKii = log_gwishart_complete_pdf(K1, b_prior + 1.0, D1);

    return logJoint + logKii - logK2;
  } else {
    arma::uvec evec(2);
    evec(0) = ii;
    evec(1) = jj;
    arma::mat C_12b = C.rows(evec);
    C_12b.shed_cols(evec);
    arma::mat C_22b = C;
    arma::uvec rem = arma::sort(evec, "descend");
    C_22b.shed_row(rem(0));
    C_22b.shed_col(rem(0));
    C_22b.shed_row(rem(1));
    C_22b.shed_col(rem(1));
    arma::mat invC_22 = safe_inv_sympd(C_22b);
    arma::mat Ce = C_12b * invC_22 * C_12b.t();
    arma::mat A = C.submat(evec, evec) - Ce;

    double sign = 0.0, logdetC = 0.0;
    arma::log_det(logdetC, sign, C);
    double logJoint = (b_prior - 2.0) / 2.0 * logdetC - 0.5 * arma::trace(D_prior * C);

    arma::mat D2 = D_prior.submat(evec, evec);
    double logK2 = log_gwishart_complete_pdf(A, b_prior, D2);

    arma::mat V = safe_inv_sympd(D2);
    double D_priorii = 1.0 / V(0, 0);
    arma::mat K1(1, 1);
    K1(0, 0) = A(0, 0);
    arma::mat D1(1, 1);
    D1(0, 0) = D_priorii;
    double logKii = log_gwishart_complete_pdf(K1, b_prior + 1.0, D1);

    return logJoint + logKii - logK2;
  }
}

// ============================================================================
// Internal maximal cliques (Bron-Kerbosch)
// ============================================================================

inline void bronk_internal(const arma::umat& adj, std::vector<int>& Rset,
                           std::vector<int> P, std::vector<int> X,
                           std::vector<std::vector<int>>& cliques) {
  if (P.empty() && X.empty()) {
    cliques.push_back(Rset);
    return;
  }
  int u = -1;
  size_t best = 0;
  std::vector<int> PX = P;
  PX.insert(PX.end(), X.begin(), X.end());
  for (int cand : PX) {
    size_t cnt = 0;
    for (int v : P)
      if (adj(cand, v)) cnt++;
    if (cnt > best) { best = cnt; u = cand; }
  }
  std::vector<int> P_without_Nu;
  for (int v : P) {
    if (u == -1 || !adj(u, v)) P_without_Nu.push_back(v);
  }
  for (int v : P_without_Nu) {
    Rset.push_back(v);
    std::vector<int> P2, X2;
    for (int w : P) if (adj(v, w)) P2.push_back(w);
    for (int w : X) if (adj(v, w)) X2.push_back(w);
    bronk_internal(adj, Rset, P2, X2, cliques);
    Rset.pop_back();
    P.erase(std::remove(P.begin(), P.end(), v), P.end());
    X.push_back(v);
  }
}

inline arma::umat maximal_cliques_internal(const arma::umat& adj0) {
  int p = (int)adj0.n_rows;
  std::vector<int> Rset;
  std::vector<int> P(p), X;
  for (int i = 0; i < p; i++) P[i] = i;
  std::vector<std::vector<int>> cliques;
  bronk_internal(adj0, Rset, P, X, cliques);
  int nc = (int)cliques.size();
  arma::umat M(p, nc, arma::fill::zeros);
  for (int c = 0; c < nc; c++) {
    for (int v : cliques[c]) M(v, c) = 1;
  }
  return M;
}

// ============================================================================
// Internal GWishart_BIPS (0-based, returns struct)
// ============================================================================

inline GWishartResult gwishart_BIPS_internal(double bG, const arma::mat& DG,
                                             const arma::umat& adj, arma::mat C,
                                             int burnin, int nmc) {
  int p = (int)DG.n_rows;
  GWishartResult fail_result;
  fail_result.C = C;
  fail_result.adj = adj;
  fail_result.ok = false;

  arma::umat adj0 = adj;
  adj0.diag().zeros();
  arma::umat cliqueMatrix = maximal_cliques_internal(adj0);
  int nc = (int)cliqueMatrix.n_cols;

  // Enforce graph sparsity: zero off-diagonal entries where no edge exists,
  // but PRESERVE diagonal (diagonal entries are always present in precision
  // matrices; zeroing them makes C singular and causes BIPS to fail).
  arma::umat adj_with_diag = adj;
  adj_with_diag.diag().ones();
  C = C % arma::conv_to<arma::mat>::from(adj_with_diag);
  C = 0.5 * (C + C.t());

  arma::mat Sig(p, p, arma::fill::zeros);

  for (int iter = 1; iter <= burnin + nmc; iter++) {
    for (int ci = 0; ci < nc; ci++) {
      arma::uvec cliqueid = arma::find(cliqueMatrix.col(ci) == 1);
      int cs = (int)cliqueid.n_elem;

      arma::mat Sigma;
      if (!safe_inv_sympd_try(Sigma, DG.submat(cliqueid, cliqueid))) return fail_result;

      arma::mat A;
      if (!wishart_rnd_try(A, Sigma, bG + cs - 1.0)) return fail_result;

      arma::mat C_12 = C.rows(cliqueid);
      C_12.shed_cols(cliqueid);

      arma::mat C_22 = C;
      for (int k = (int)cliqueid.n_elem - 1; k >= 0; k--) {
        C_22.shed_row(cliqueid(k));
        C_22.shed_col(cliqueid(k));
      }

      arma::mat invC_22;
      if (!safe_inv_sympd_try(invC_22, C_22)) return fail_result;

      arma::mat C121 = C_12 * invC_22 * C_12.t();
      C121 = 0.5 * (C121 + C121.t());

      arma::mat K_c = A + C121;
      K_c = 0.5 * (K_c + K_c.t());

      C.submat(cliqueid, cliqueid) = K_c;
    }

    if (iter > burnin) {
      if (!safe_inv_sympd_try(Sig, C)) return fail_result;
    }
  }

  GWishartResult result;
  result.C = C;
  result.Sig = Sig;
  result.adj = adj;
  result.ok = true;
  return result;
}

// ============================================================================
// Internal GWishart_NOij_Gibbs (0-based indices ii, jj)
// ============================================================================

inline GWishartResult gwishart_NOij_gibbs_internal(double bG, const arma::mat& DG,
                                                    arma::umat adj, arma::mat C,
                                                    int ii, int jj, int edgeij,
                                                    int burnin, int nmc) {
  int p = (int)DG.n_rows;
  GWishartResult fail_result;
  fail_result.C = C;
  fail_result.adj = adj;
  fail_result.ok = false;

  C = 0.5 * (C + C.t());
  arma::umat adjt = adj.t();
  arma::umat sum_adj = adj + adjt;
  sum_adj.transform([](arma::uword x) { return (x > 0u) ? 1u : 0u; });
  adj = sum_adj;
  adj.diag().zeros();

  if (edgeij == 0) {
    C(ii, jj) = 0;
    C(jj, ii) = 0;
    double A = R::rgamma(bG / 2.0, 2.0 / DG(jj, jj));
    arma::rowvec C_12 = C.row(jj);
    C_12.shed_col(jj);
    arma::mat C_22 = C;
    C_22.shed_row(jj);
    C_22.shed_col(jj);
    arma::mat invC_22;
    if (!safe_inv_sympd_try(invC_22, C_22)) return fail_result;
    double c = arma::as_scalar(C_12 * invC_22 * C_12.t());
    C(jj, jj) = A + c;
  } else {
    arma::uvec reorder(p);
    int idx = 0;
    for (int k = 0; k < p; k++)
      if (k != ii && k != jj) reorder(idx++) = k;
    reorder(idx++) = ii;
    reorder(idx++) = jj;
    arma::mat C_re = C.submat(reorder, reorder);
    arma::mat R;
    if (!safe_chol_try(R, C_re, "upper", 1e-10, 10)) return fail_result;

    double b = R(p - 2, p - 2);
    double m_post = -b * DG(ii, jj) / DG(jj, jj);
    double sig_post = 1.0 / std::sqrt(DG(jj, jj));
    adj(ii, jj) = 1;
    adj(jj, ii) = 1;
    R(p - 2, p - 1) = R::rnorm(m_post, sig_post);
    R(p - 1, p - 1) = std::sqrt(R::rgamma(bG / 2.0, 2.0 / DG(jj, jj)));
    arma::vec col_end = R.col(p - 1);
    arma::mat cols = R.cols(p - 2, p - 1);
    arma::vec C_updated = cols.t() * col_end;
    C(ii, jj) = C_updated(0);
    C(jj, ii) = C_updated(0);
    C(jj, jj) = C_updated(1);
  }

  arma::mat Sig;
  if (!safe_inv_sympd_try(Sig, C)) return fail_result;

  arma::uvec deg = arma::conv_to<arma::uvec>::from(arma::sum(adj, 1));
  arma::uvec isolated = arma::find(deg == 0);
  for (unsigned int t = 0; t < isolated.n_elem; t++) {
    int node = isolated(t);
    double Kc = R::rgamma(bG / 2.0, 2.0 / DG(node, node));
    C(node, node) = Kc;
    Sig(node, node) = 1.0 / Kc;
  }

  for (int iter = 1; iter <= burnin + nmc; iter++) {
    for (int a = 0; a < p - 1; a++) {
      for (int bidx = a + 1; bidx < p; bidx++) {
        if (adj(a, bidx) == 1 && !(a == ii && bidx == jj)) {
          arma::uvec cliqueid(2);
          cliqueid(0) = a;
          cliqueid(1) = bidx;
          arma::mat Sigma;
          if (!safe_inv_sympd_try(Sigma, DG.submat(cliqueid, cliqueid))) return fail_result;
          arma::mat Aw;
          if (!wishart_rnd_try(Aw, Sigma, bG + 1.0)) return fail_result;

          arma::mat C_12 = C.rows(cliqueid);
          C_12.shed_cols(cliqueid);

          arma::mat Sig_12 = Sig.rows(cliqueid);
          Sig_12.shed_cols(cliqueid);
          arma::mat Sig_22 = Sig;
          Sig_22.shed_rows(cliqueid);
          Sig_22.shed_cols(cliqueid);
          arma::mat invSig_11;
          if (!safe_inv_sympd_try(invSig_11, Sig.submat(cliqueid, cliqueid))) return fail_result;
          invSig_11 = 0.5 * (invSig_11 + invSig_11.t());
          arma::mat invC_22 = Sig_22 - Sig_12.t() * invSig_11 * Sig_12;

          arma::mat K_c = Aw + C_12 * invC_22 * C_12.t();
          K_c = 0.5 * (K_c + K_c.t());

          arma::mat Delta;
          // NB: C_old - K_c is a difference of PD matrices, NOT necessarily PD.
          // Must use general inverse, not inv_sympd.
          if (!safe_inv_general_try(Delta, C.submat(cliqueid, cliqueid) - K_c)) return fail_result;
          C.submat(cliqueid, cliqueid) = K_c;

          arma::mat Sig_bb = Sig.submat(cliqueid, cliqueid);
          arma::mat aa;
          // NB: Delta - Sig_bb is also not necessarily PD.
          if (!safe_inv_general_try(aa, Delta - Sig_bb)) return fail_result;
          aa = 0.5 * (aa + aa.t());
          Sig = Sig + Sig.cols(cliqueid) * aa * Sig.rows(cliqueid);
        }
      }
    }
  }

  GWishartResult result;
  result.C = C;
  result.Sig = Sig;
  result.adj = adj;
  result.ok = true;
  return result;
}

// ============================================================================
// Internal calc_mrf_logC (arma types, no Rcpp wrapper overhead)
// ============================================================================

struct MRFBinaryCacheEntry {
  arma::mat B;   // (2^K) x K
  arma::vec s;   // row sums
};

// Use a function-local static to avoid ODR issues with multiple TUs
inline const MRFBinaryCacheEntry& get_mrf_binary_cache(int K) {
  static std::map<int, MRFBinaryCacheEntry> cache;
  auto it = cache.find(K);
  if (it != cache.end()) return it->second;

  std::size_t n = static_cast<std::size_t>(1) << K;
  MRFBinaryCacheEntry entry;
  entry.B.set_size(n, K);
  entry.B.zeros();
  for (std::size_t r = 0; r < n; ++r) {
    std::size_t v = r;
    for (int c = 0; c < K; ++c) {
      entry.B(r, c) = static_cast<double>(v & 1U);
      v >>= 1U;
    }
  }
  entry.s = arma::sum(entry.B, 1);
  auto ins = cache.emplace(K, std::move(entry));
  return ins.first->second;
}

// Scalar version: compute logC for a single nu value
inline double calc_mrf_logC_scalar(const arma::mat& Theta, double nu_val) {
  int K = (int)Theta.n_rows;
  const MRFBinaryCacheEntry& cache = get_mrf_binary_cache(K);
  arma::vec q = arma::sum((cache.B * Theta) % cache.B, 1);
  arma::vec tmp = nu_val * cache.s + q;
  return logsumexp_internal(tmp);
}

// Vector version: compute logC for a vector of nu values
inline arma::vec calc_mrf_logC_vec(const arma::mat& Theta, const arma::vec& nu_vec) {
  int K = (int)Theta.n_rows;
  const MRFBinaryCacheEntry& cache = get_mrf_binary_cache(K);
  arma::vec q = arma::sum((cache.B * Theta) % cache.B, 1);
  arma::vec out(nu_vec.n_elem);
  for (arma::uword i = 0; i < nu_vec.n_elem; ++i) {
    arma::vec tmp = nu_vec(i) * cache.s + q;
    out(i) = logsumexp_internal(tmp);
  }
  return out;
}

#endif // WANGLI_INTERNAL_H
