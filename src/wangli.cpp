
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// -------- utilities --------

static inline double log_multi_gamma(int p, double a){
  // log Γ_p(a) = (p(p-1)/4) log(pi) + sum_{j=1}^p lgamma(a + (1-j)/2)
  double out = (double)p*(p-1)*0.25*std::log(M_PI);
  for(int j=1;j<=p;j++){
    out += ::lgamma(a + (1.0 - (double)j)/2.0);
  }
  return out;
}

static inline double log_gwishart_complete_const(double b, const arma::mat& D){
  // For complete graph with p = dim(D).
  // Density: p(K) ∝ |K|^{(b-2)/2} exp(-0.5 tr(D K))
  // Equivalent to Wishart(ν=b+p-1, Σ=inv(D))
  int p = (int)D.n_rows;
  double nu = b + p - 1.0;
  double sign=0.0, logdetD=0.0;
  arma::log_det(logdetD, sign, D);
  // sign should be +1 if PD
  double out = (nu * p / 2.0) * std::log(2.0) - (nu/2.0)*logdetD + log_multi_gamma(p, nu/2.0);
  return out;
}

static inline double log_gwishart_complete_pdf(const arma::mat& K, double b, const arma::mat& D){
  double sign=0.0, logdetK=0.0;
  arma::log_det(logdetK, sign, K);
  double trDK = arma::trace(D*K);
  double logZ = log_gwishart_complete_const(b, D);
  return ((b-2.0)/2.0)*logdetK - 0.5*trDK - logZ;
}

//Another helper
static inline bool safe_chol_try(arma::mat& L,
                                 const arma::mat& Ain,
                                 const char* which = "lower",
                                 double jitter0 = 1e-10,
                                 int max_tries = 10) {
  arma::mat A = 0.5 * (Ain + Ain.t()); // enforce symmetry

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


// Armadillo's inv_sympd() throws if the matrix is not numerically SPD.
// In MCMC updates it is common to encounter matrices that are theoretically SPD
// but fail due to floating point tolerances (near-zero eigenvalues).
// This helper attempts a symmetric-PD inverse, adding a small diagonal jitter
// if needed. It fails fast with an informative error if inversion is impossible.

// Attempt a symmetric positive definite inverse without throwing.
// Returns true on success, false if not SPD even after adding diagonal jitter.
static inline bool safe_inv_sympd_try(arma::mat& Ainv,
                                      const arma::mat& A_in,
                                      double jitter0 = 1e-10,
                                      int max_tries = 10) {
  arma::mat A = 0.5 * (A_in + A_in.t());

  // Fast path for 1x1
  if (A.n_rows == 1u && A.n_cols == 1u) {
    double a = A(0,0);
    if (!(a > 0.0)) a += jitter0;
    if (!(a > 0.0)) return false;
    Ainv.set_size(1u,1u);
    Ainv(0,0) = 1.0 / a;
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

// Backwards-compatible wrapper: stop() if inversion fails.
static inline arma::mat safe_inv_sympd(const arma::mat& A_in,
                                      double jitter0 = 1e-10,
                                      int max_tries = 10){
  arma::mat Ainv;
  bool ok = safe_inv_sympd_try(Ainv, A_in, jitter0, max_tries);
  if(!ok){
    Rcpp::stop("safe_inv_sympd: matrix not SPD even after diagonal jitter");
  }
  return Ainv;
}


static inline bool wishart_rnd_try(arma::mat& W_out,
                                   const arma::mat& Sigma,
                                   double df,
                                   double jitter0 = 1e-10,
                                   int max_tries = 10) {

  int p = (int)Sigma.n_rows;

  arma::mat L;
  if (!safe_chol_try(L, Sigma, "lower", jitter0, max_tries)) {
    return false;
  }

  arma::mat A(p, p, arma::fill::zeros);
  for (int i = 0; i < p; i++) {
    A(i,i) = std::sqrt(R::rchisq(df - i));
    for (int j = 0; j < i; j++) {
      A(i,j) = R::rnorm(0.0, 1.0);
    }
  }

  arma::mat LA = L * A;
  arma::mat W  = LA * LA.t();
  W_out = 0.5 * (W + W.t());
  return true;
}


// [[Rcpp::export]]
arma::mat wishrnd_cpp(const arma::mat& Sigma, double df){
  arma::mat W;
  bool ok = wishart_rnd_try(W, Sigma, df);
  if(!ok){
    Rcpp::stop("wishrnd_cpp: chol failed (Sigma not SPD even after jitter)");
  }
  return W;
}


// [[Rcpp::export]]
double log_iwishart_invA_const_cpp(double df, const arma::mat& S){
  int p = (int)S.n_rows;
  double sign=0.0, logdetS=0.0;
  arma::log_det(logdetS, sign, S);
  double out = (df + p - 1.0)/2.0 * (logdetS - p*std::log(2.0)) - log_multi_gamma(p, (df+p-1.0)/2.0);
  return out;
}

// [[Rcpp::export]]
double log_J_cpp(double h, const arma::mat& B, double a11){
  // Wang/Li log_J (uses scalar InvWishart constant)
  // B expected 2x2 (but used only elements)
  double B22 = B(1,1); // careful: in MATLAB B(2,2). Here 0-based => (1,1)
  double B11 = B(0,0);
  double B12 = B(0,1);
  // scalar inverse-wishart constant with p=1: S is 1x1
  arma::mat S(1,1); S(0,0)=B22;
  double term = 0.5*std::log(2.0*M_PI / B22)
    - log_iwishart_invA_const_cpp(h, S)
    + (h-1.0)/2.0 * std::log(a11)
    - 0.5*(B11 - (B12*B12)/B22)*a11;
  return term;
}

// [[Rcpp::export]]
double log_H_cpp(double b_prior, const arma::mat& D_prior, double n, const arma::mat& S, const arma::mat& C, int i, int j){
  // i,j are 1-based indices from R
  int ii=i-1, jj=j-1;
  // int p = (int)C.n_rows; // unused

  // (i,j)=0 case: set edge to 0, compute c
  arma::mat C0 = C;
  C0(ii,jj)=0; C0(jj,ii)=0;
  int e = jj; // 0-based
  arma::rowvec C_12 = C0.row(e);
  C_12.shed_col(e);
  arma::mat C_22 = C0;
  C_22.shed_row(e); C_22.shed_col(e);
  arma::mat invC_22 = safe_inv_sympd(C_22);
  double c = as_scalar(C_12 * invC_22 * C_12.t());
  arma::mat C0_ij(2,2,fill::zeros);
  C0_ij(0,0) = C(ii,ii);
  C0_ij(1,1) = c;

  // (i,j)=1 case
  arma::uvec evec(2); evec(0)=ii; evec(1)=jj;
  arma::mat C_12b = C.rows(evec);
  C_12b.shed_cols(evec);
  arma::mat C_22b = C;
  // remove rows/cols ii and jj (note order: remove larger index first)
  arma::uvec rem = arma::sort(evec, "descend");
  C_22b.shed_row(rem(0)); C_22b.shed_col(rem(0));
  C_22b.shed_row(rem(1)); C_22b.shed_col(rem(1));
  arma::mat invC_22b = safe_inv_sympd(C_22b);
  arma::mat Ce = C_12b * invC_22b * C_12b.t();
  arma::mat A = C.submat(evec,evec) - Ce;
  double a11 = A(0,0);
  arma::mat C1_ij = Ce;

  double b_post = b_prior + n;
  arma::mat D_post = D_prior + S;

  // scalar iwishart const on D_post(j,j)
  arma::mat Ssc(1,1); Ssc(0,0)=D_post(jj,jj);
  double term1 = -log_iwishart_invA_const_cpp(b_post, Ssc);
  double term2 = -log_J_cpp(b_post, D_post.submat(evec,evec), a11);
  double term3 = (n + b_prior - 2.0)/2.0 * std::log(a11);
  double trterm = arma::trace(D_post.submat(evec,evec) * (C0_ij - C1_ij));
  double term4 = -0.5 * trterm;
  return term1 + term2 + term3 + term4;
}

// [[Rcpp::export]]
double log_GWishart_NOij_pdf_cpp(double b_prior, const arma::mat& D_prior, arma::mat C, int i, int j, int edgeij){
  int ii=i-1, jj=j-1;
  // int p = (int)C.n_rows; // unused

  if(edgeij==0){
    C(ii,jj)=0; C(jj,ii)=0;
    // remove j
    arma::rowvec C_12 = C.row(jj);
    C_12.shed_col(jj);
    arma::mat C_22 = C;
    C_22.shed_row(jj); C_22.shed_col(jj);
    arma::mat invC_22 = safe_inv_sympd(C_22);
    double c = as_scalar(C_12 * invC_22 * C_12.t());
    arma::mat C0_ij(2,2,fill::zeros);
    C0_ij(0,0)=C(ii,ii);
    C0_ij(1,1)=c;

    double sign=0.0, logdetC=0.0;
    arma::log_det(logdetC, sign, C);
    double logJoint = (b_prior-2.0)/2.0*logdetC - 0.5*arma::trace(D_prior*C);

    // cliqueid = [i j]
    arma::uvec evec(2); evec(0)=ii; evec(1)=jj;
    arma::mat A = C0_ij;
    arma::mat D2 = D_prior.submat(evec,evec);
    // adjacency complete for clique => complete wishart with p=2
    double logK2 = log_gwishart_complete_pdf(A, b_prior, D2);

    // Kii: p=1 with b_prior+1 and D_priorii = 1/V(1,1) where V=inv(D2)
    arma::mat V = safe_inv_sympd(D2);
    double D_priorii = 1.0 / V(0,0);
    arma::mat K1(1,1); K1(0,0)=A(0,0);
    arma::mat D1(1,1); D1(0,0)=D_priorii;
    double logKii = log_gwishart_complete_pdf(K1, b_prior+1.0, D1);

    return logJoint + logKii - logK2;
  } else {
    // edgeij==1 case: similar but uses A=C(e,e)-Ce
    arma::uvec evec(2); evec(0)=ii; evec(1)=jj;
    arma::mat C_12b = C.rows(evec);
    C_12b.shed_cols(evec);
    arma::mat C_22b = C;
    arma::uvec rem = arma::sort(evec, "descend");
    C_22b.shed_row(rem(0)); C_22b.shed_col(rem(0));
    C_22b.shed_row(rem(1)); C_22b.shed_col(rem(1));
    arma::mat invC_22 = safe_inv_sympd(C_22b);
    arma::mat Ce = C_12b * invC_22 * C_12b.t();
    arma::mat A = C.submat(evec,evec) - Ce;

    double sign=0.0, logdetC=0.0;
    arma::log_det(logdetC, sign, C);
    double logJoint = (b_prior-2.0)/2.0*logdetC - 0.5*arma::trace(D_prior*C);

    arma::mat D2 = D_prior.submat(evec,evec);
    double logK2 = log_gwishart_complete_pdf(A, b_prior, D2);

    arma::mat V = safe_inv_sympd(D2);
    double D_priorii = 1.0 / V(0,0);
    arma::mat K1(1,1); K1(0,0)=A(0,0);
    arma::mat D1(1,1); D1(0,0)=D_priorii;
    double logKii = log_gwishart_complete_pdf(K1, b_prior+1.0, D1);

    return logJoint + logKii - logK2;
  }
}

// ---- maximal cliques (Bron–Kerbosch with pivot) ----
static void bronk(const arma::umat& adj, std::vector<int>& Rset,
                  std::vector<int> P, std::vector<int> X,
                  std::vector< std::vector<int> >& cliques){
  if(P.empty() && X.empty()){
    cliques.push_back(Rset);
    return;
  }
  // choose pivot u from P ∪ X with max neighbors in P
  int u=-1; size_t best=0;
  std::vector<int> PX = P; PX.insert(PX.end(), X.begin(), X.end());
  for(int cand: PX){
    size_t cnt=0;
    for(int v: P) if(adj(cand,v)) cnt++;
    if(cnt>best){best=cnt; u=cand;}
  }
  // iterate v in P \ N(u)
  std::vector<int> P_without_Nu;
  for(int v: P){
    if(u==-1 || !adj(u,v)) P_without_Nu.push_back(v);
  }
  for(int v: P_without_Nu){
    Rset.push_back(v);
    std::vector<int> P2, X2;
    for(int w: P) if(adj(v,w)) P2.push_back(w);
    for(int w: X) if(adj(v,w)) X2.push_back(w);
    bronk(adj, Rset, P2, X2, cliques);
    Rset.pop_back();
    // move v from P to X
    P.erase(std::remove(P.begin(), P.end(), v), P.end());
    X.push_back(v);
  }
}

// [[Rcpp::export]]
arma::umat maximal_cliques_cpp(const arma::umat& adj0){
  // adj0: p x p 0/1 with diagonal 0
  int p = (int)adj0.n_rows;
  std::vector<int> Rset;
  std::vector<int> P(p), X;
  for(int i=0;i<p;i++) P[i]=i;
  std::vector< std::vector<int> > cliques;
  bronk(adj0, Rset, P, X, cliques);
  int nc = (int)cliques.size();
  arma::umat M(p, nc, fill::zeros);
  for(int c=0;c<nc;c++){
    for(int v: cliques[c]) M(v,c)=1;
  }
  return M;
}

// [[Rcpp::export]]

Rcpp::List GWishart_BIPS_maximumClique_cpp(double bG, const arma::mat& DG, const arma::umat& adj, arma::mat C, int burnin, int nmc){
  // Port of Wang/Li GWishart_BIPS_maximumClique with robust SPD handling.
  int p = (int)DG.n_rows;

  arma::umat adj_orig = adj;
  arma::mat  C_orig  = C;

  arma::umat adj0 = adj;
  adj0.diag().zeros();
  arma::umat cliqueMatrix = maximal_cliques_cpp(adj0);
  int nc = (int)cliqueMatrix.n_cols;

  // Respect sparsity pattern implied by adj (off-diagonal). Keep diagonal as-is.
  C = C % arma::conv_to<arma::mat>::from(adj);
  C = 0.5*(C + C.t());

  arma::mat Sig(p,p,fill::zeros);

  for(int iter=1; iter<=burnin+nmc; iter++){
    for(int ci=0; ci<nc; ci++){
      arma::uvec cliqueid = arma::find(cliqueMatrix.col(ci)==1);
      int cs = (int)cliqueid.n_elem;

      arma::mat Sigma;
      if(!safe_inv_sympd_try(Sigma, DG.submat(cliqueid, cliqueid))) {
        return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false);
      }

      arma::mat A;
      if (!wishart_rnd_try(A, Sigma, bG + cs - 1.0)) {
        return Rcpp::List::create(
          Rcpp::_["C"]   = C_orig,
          Rcpp::_["Sig"] = arma::mat(),
          Rcpp::_["adj"] = adj_orig,
          Rcpp::_["ok"]  = false
        );
      }

      arma::mat C_12 = C.rows(cliqueid);
      C_12.shed_cols(cliqueid);

      arma::mat C_22 = C;
      for(int k=(int)cliqueid.n_elem-1; k>=0; k--){
        C_22.shed_row(cliqueid(k));
        C_22.shed_col(cliqueid(k));
      }

      arma::mat invC_22;
      if(!safe_inv_sympd_try(invC_22, C_22)) {
        return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false);
      }

      arma::mat C121 = C_12 * invC_22 * C_12.t();
      C121 = 0.5*(C121 + C121.t());

      arma::mat K_c = A + C121;
      K_c = 0.5*(K_c + K_c.t());

      C.submat(cliqueid, cliqueid) = K_c;
    }

    if(iter > burnin){
      if(!safe_inv_sympd_try(Sig, C)) {
        return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false);
      }
    }
  }

  return Rcpp::List::create(_["C"]=C, _["Sig"]=Sig, _["adj"]=adj_orig, _["ok"]=true);
}


// [[Rcpp::export]]
Rcpp::List GWishart_NOij_Gibbs_cpp(double bG, const arma::mat& DG, arma::umat adj, arma::mat C, int i, int j, int edgeij, int burnin, int nmc){
  // Port of Wang/Li GWishart_NOij_Gibbs. This function returns last C and Sig after iterations.
  int p = (int)DG.n_rows;
  int ii=i-1, jj=j-1;

  arma::umat adj_orig = adj;
  arma::mat  C_orig  = C;
  // Ensure symmetry in C and adj
  C = 0.5*(C + C.t());
  arma::umat adjt = adj.t();
  arma::umat sum = adj + adjt;
  sum.transform( [](arma::uword x) { return (x > 0u) ? 1u : 0u; } );
  adj = sum;
  adj.diag().zeros();

  // Apply constraint edgeij on (i,j)
  if(edgeij==0){
    C(ii,jj)=0; C(jj,ii)=0;
    arma::mat Sigma;
    if(!safe_inv_sympd_try(Sigma,  arma::mat(1,1,fill::value( DG(jj,jj) )) )) { return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false); }
    // sample diagonal C(j,j)
    double A = R::rgamma(bG/2.0, 2.0/DG(jj,jj)); // gamma(shape, scale)
    arma::rowvec C_12 = C.row(jj); C_12.shed_col(jj);
    arma::mat C_22 = C; C_22.shed_row(jj); C_22.shed_col(jj);
    arma::mat invC_22;
    if(!safe_inv_sympd_try(invC_22, C_22)) { return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false); }
    double c = as_scalar(C_12 * invC_22 * C_12.t());
    C(jj,jj) = A + c;
  } else {
    // edge present: reorder columns to have i,j last
    arma::uvec reorder(p);
    int idx=0;
    for(int k=0;k<p;k++) if(k!=ii && k!=jj) reorder(idx++)=k;
    reorder(idx++)=ii; reorder(idx++)=jj;
    arma::mat C_re = C.submat(reorder,reorder);
    arma::mat R;
    if (!safe_chol_try(R, C_re, "upper", 1e-10, 10)) {
      return Rcpp::List::create(
        Rcpp::_["C"]   = C_orig,
        Rcpp::_["Sig"] = arma::mat(),
        Rcpp::_["adj"] = adj_orig,
        Rcpp::_["ok"]  = false
      );
    }

    double b = R(p-2,p-2); // p-1 in matlab => p-2
    double m_post = -b*DG(ii,jj)/DG(jj,jj);
    double sig_post = 1.0/std::sqrt(DG(jj,jj));
    adj(ii,jj)=1; adj(jj,ii)=1;
    R(p-2,p-1) = R::rnorm(m_post, sig_post);
    R(p-1,p-1) = std::sqrt(R::rgamma(bG/2.0, 2.0/DG(jj,jj)));
    // C_updated = R(:,end-1:end)'*R(:,end);
    arma::vec col_end = R.col(p-1);
    arma::mat cols = R.cols(p-2,p-1);
    arma::vec C_updated = cols.t()*col_end;
    C(ii,jj)=C_updated(0); C(jj,ii)=C_updated(0);
    C(jj,jj)=C_updated(1);
  }

  arma::mat Sig;
  if(!safe_inv_sympd_try(Sig, C)) { return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false); }
  // isolated nodes
  // NOTE: qualify with arma::sum to avoid Rcpp sugar `sum()` overloads.
  // For `arma::umat`, `arma::sum(adj, 1)` returns an Armadillo column of uword;
  // convert explicitly to `arma::uvec`.
  arma::uvec deg = arma::conv_to<arma::uvec>::from(arma::sum(adj, 1)); // includes diagonal ones
  arma::uvec isolated = find(deg==0);
  for(unsigned int t=0;t<isolated.n_elem;t++){
    int node = isolated(t);
    double Kc = R::rgamma(bG/2.0, 2.0/DG(node,node));
    C(node,node)=Kc;
    Sig(node,node)=1.0/Kc;
  }

  for(int iter=1; iter<=burnin+nmc; iter++){
    for(int a=0;a<p-1;a++){
      for(int bidx=a+1;bidx<p;bidx++){
        if(adj(a,bidx)==1 && !(a==ii && bidx==jj)){
          arma::uvec cliqueid(2); cliqueid(0)=a; cliqueid(1)=bidx;
          arma::mat Sigma;
          if(!safe_inv_sympd_try(Sigma, DG.submat(cliqueid,cliqueid))) { return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false); }
          arma::mat A;
          if (!wishart_rnd_try(A, Sigma, bG + 1.0)) {
            return Rcpp::List::create(
              Rcpp::_["C"]   = C_orig,
              Rcpp::_["Sig"] = arma::mat(),
              Rcpp::_["adj"] = adj_orig,
              Rcpp::_["ok"]  = false
            );
          }
          arma::mat C_12 = C.rows(cliqueid);
          C_12.shed_cols(cliqueid);

          arma::mat Sig_12 = Sig.rows(cliqueid);
          Sig_12.shed_cols(cliqueid);
          arma::mat Sig_22 = Sig;
          Sig_22.shed_rows(cliqueid);
          Sig_22.shed_cols(cliqueid);
          arma::mat invSig_11;
          if(!safe_inv_sympd_try(invSig_11, Sig.submat(cliqueid,cliqueid))) { return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false); }
          invSig_11 = 0.5*(invSig_11 + invSig_11.t());
          arma::mat invC_22 = Sig_22 - Sig_12.t()*invSig_11*Sig_12;

          arma::mat K_c = A + C_12*invC_22*C_12.t();
          K_c = 0.5*(K_c + K_c.t());

          arma::mat Delta;
          if(!safe_inv_sympd_try(Delta,  C.submat(cliqueid,cliqueid) - K_c )) { return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false); }
          C.submat(cliqueid,cliqueid) = K_c;

          arma::mat Sig_bb = Sig.submat(cliqueid,cliqueid);
          arma::mat aa;
          if(!safe_inv_sympd_try(aa,  Delta - Sig_bb )) { return Rcpp::List::create(_["C"]=C_orig, _["Sig"]=arma::mat(), _["adj"]=adj_orig, _["ok"]=false); }
          aa = 0.5*(aa + aa.t());
          Sig = Sig + Sig.cols(cliqueid)*aa*Sig.rows(cliqueid);
        }
      }
    }
  }

  return Rcpp::List::create(_["C"]=C, _["Sig"]=Sig, _["adj"]=adj, _["ok"]=true);
}
