
#include <RcppArmadillo.h>
#include <map>
// [[Rcpp::depends(RcppArmadillo)]]

inline double logsumexp(const arma::vec& x) {
  double xmax = x.max();
  return xmax + std::log(arma::sum(arma::exp(x - xmax)));
}

// Cache of all binary matrices and their row sums by K.
// This avoids regenerating the 2^K x K matrix on every call.
struct BinaryCacheEntry {
  arma::mat B;   // (2^K) x K
  arma::vec s;   // row sums
};

static std::map<int, BinaryCacheEntry> BIN_CACHE;

static const BinaryCacheEntry& get_binary_cache(const int K) {
  auto it = BIN_CACHE.find(K);
  if (it != BIN_CACHE.end()) return it->second;

  std::size_t n = static_cast<std::size_t>(1) << K;
  BinaryCacheEntry entry;
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

  auto ins = BIN_CACHE.emplace(K, std::move(entry));
  return ins.first->second;
}

//' Compute log normalizing constant for the edge-sharing MRF prior
//'
//' Computes \eqn{\log C(\Theta, \nu)} where
//' \deqn{C(\Theta,\nu)=\sum_{g\in\{0,1\}^K} \exp\{\nu \mathbf{1}^\top g + g^\top \Theta g\}.}
//'
//' @param Theta A K x K symmetric matrix.
//' @param nu Numeric vector of length 1 or more.
//' @return Numeric vector of \eqn{\log C(\Theta,\nu)} for each entry of \code{nu}.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calc_mrf_logC_cpp(const arma::mat& Theta,
                                      const Rcpp::NumericVector& nu) {
  int K = Theta.n_rows;
  if (Theta.n_cols != K) Rcpp::stop("Theta must be square.");
  if (K <= 0) Rcpp::stop("Theta must have positive dimension.");
  if (K > 30) Rcpp::stop("K too large for exact enumeration; 2^K would overflow.");

  const BinaryCacheEntry& cache = get_binary_cache(K);
  const arma::mat& B = cache.B;
  const arma::vec& s = cache.s;

  // q_r = g_r' Theta g_r for each row r
  arma::vec q = arma::sum((B * Theta) % B, 1);

  Rcpp::NumericVector out(nu.size());
  arma::vec tmp(B.n_rows);

  for (R_xlen_t i = 0; i < nu.size(); ++i) {
    tmp = static_cast<double>(nu[i]) * s + q;
    out[i] = logsumexp(tmp);
  }
  return out;
}
