// Rcpp-exported wrappers around internal 0-based functions.
// All implementation is in wangli_internal.h.

#include "wangli_internal.h"

// [[Rcpp::export]]
arma::mat wishrnd_cpp(const arma::mat& Sigma, double df) {
  arma::mat W;
  bool ok = wishart_rnd_try(W, Sigma, df);
  if (!ok) Rcpp::stop("wishrnd_cpp: chol failed (Sigma not SPD even after jitter)");
  return W;
}

// [[Rcpp::export]]
double log_iwishart_invA_const_cpp(double df, const arma::mat& S) {
  return log_iwishart_invA_const_internal(df, S);
}

// [[Rcpp::export]]
double log_J_cpp(double h, const arma::mat& B, double a11) {
  return log_J_internal(h, B, a11);
}

// 1-based wrapper: converts i,j from R's 1-based to 0-based
// [[Rcpp::export]]
double log_H_cpp(double b_prior, const arma::mat& D_prior, double n,
                 const arma::mat& S, const arma::mat& C, int i, int j) {
  return compute_log_H(b_prior, D_prior, n, S, C, i - 1, j - 1);
}

// 1-based wrapper
// [[Rcpp::export]]
double log_GWishart_NOij_pdf_cpp(double b_prior, const arma::mat& D_prior,
                                  arma::mat C, int i, int j, int edgeij) {
  return compute_log_GWishart_NOij_pdf(b_prior, D_prior, C, i - 1, j - 1, edgeij);
}

// [[Rcpp::export]]
arma::umat maximal_cliques_cpp(const arma::umat& adj0) {
  return maximal_cliques_internal(adj0);
}

// [[Rcpp::export]]
Rcpp::List GWishart_BIPS_maximumClique_cpp(double bG, const arma::mat& DG,
                                            const arma::umat& adj, arma::mat C,
                                            int burnin, int nmc) {
  GWishartResult res = gwishart_BIPS_internal(bG, DG, adj, C, burnin, nmc);
  return Rcpp::List::create(
    Rcpp::_["C"] = res.C,
    Rcpp::_["Sig"] = res.ok ? res.Sig : arma::mat(),
    Rcpp::_["adj"] = res.adj,
    Rcpp::_["ok"] = res.ok
  );
}

// 1-based wrapper for i,j
// [[Rcpp::export]]
Rcpp::List GWishart_NOij_Gibbs_cpp(double bG, const arma::mat& DG,
                                    arma::umat adj, arma::mat C,
                                    int i, int j, int edgeij,
                                    int burnin, int nmc) {
  GWishartResult res = gwishart_NOij_gibbs_internal(bG, DG, adj, C,
                                                     i - 1, j - 1, edgeij,
                                                     burnin, nmc);
  return Rcpp::List::create(
    Rcpp::_["C"] = res.C,
    Rcpp::_["Sig"] = res.ok ? res.Sig : arma::mat(),
    Rcpp::_["adj"] = res.adj,
    Rcpp::_["ok"] = res.ok
  );
}
