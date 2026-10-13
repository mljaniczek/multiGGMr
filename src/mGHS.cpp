// This file is part of the mGHS sampler by Claudio Busatto
// (https://github.com/cbusatto/mGHS, commit 1c9fd70), licensed GPL (>= 2).
// Modified for multiGGM; every change is listed in inst/mghs/UPSTREAM.md.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>     
#include <vector>

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]] 

#include "boost/multi_array.hpp"

#include "utils_rgamma3p.h"
#include "utils_mGHS.h"
#include "utils_rIW.h"

typedef boost::multi_array<double, 1> vector_type;
typedef boost::multi_array<double, 2> matrix_type;
typedef boost::multi_array<double, 3> array_type;


// multiGGM: running sums of a group difference D for every edge and group
// pair, kept over all post-burn-in iterations.
struct DiffAcc
{
  int npe, ngp, nd;
  std::vector<double> deltas, sum, sumsq, n_pos, n_neg, n_exceed;
  
  DiffAcc (int npe_, int ngp_, const std::vector<double>& deltas_)
    : npe(npe_), ngp(ngp_), nd(deltas_.size()), deltas(deltas_),
      sum(npe_ * ngp_, 0.0), sumsq(npe_ * ngp_, 0.0), n_pos(npe_ * ngp_, 0.0),
      n_neg(npe_ * ngp_, 0.0), n_exceed(npe_ * deltas_.size() * ngp_, 0.0) {}
  
  void add (int g, int e, double x)
  {
    const int k = e + npe * g;
    sum[k] += x;
    sumsq[k] += x * x;
    if (x > 0.0) n_pos[k] += 1.0;
    else if (x < 0.0) n_neg[k] += 1.0;
    const double ax = std::abs(x);
    for (int t = 0; t < nd; t++)
      if (ax > deltas[t]) n_exceed[e + npe * (t + nd * g)] += 1.0;
  }
  
  Rcpp::List to_list () const
  {
    Rcpp::NumericVector s(sum.begin(), sum.end()), s2(sumsq.begin(), sumsq.end());
    Rcpp::NumericVector np(n_pos.begin(), n_pos.end()), nn(n_neg.begin(), n_neg.end());
    Rcpp::NumericVector ne(n_exceed.begin(), n_exceed.end());
    Rcpp::IntegerVector d2 = Rcpp::IntegerVector::create(npe, ngp);
    s.attr("dim") = d2; s2.attr("dim") = d2; np.attr("dim") = d2; nn.attr("dim") = d2;
    ne.attr("dim") = Rcpp::IntegerVector::create(npe, nd, ngp);
    return Rcpp::List::create(Rcpp::Named("sum") = s, Rcpp::Named("sumsq") = s2,
                              Rcpp::Named("n_pos") = np, Rcpp::Named("n_neg") = nn,
                              Rcpp::Named("n_exceed") = ne);
  }
};


// reason: 1 local-scale guard (alpha), 5 local-scale guard (beta),
// 2 local-scale draw, 3 global-scale guard, 4 global-scale draw
static Rcpp::List stuck_result (const int& b, const int& reason)
{
  return Rcpp::List::create(Rcpp::Named("stuck") = 1, Rcpp::Named("stuck_iter") = b + 1,
                            Rcpp::Named("stuck_reason") = reason,
                            Rcpp::Named("reject_tries") = mghs_reject_tries);
}


// [[Rcpp::export]]
Rcpp::List mghs_engine_cpp (const int& burnin, const int& nsave, const int& thin, const Eigen::VectorXd& n, const Rcpp::List& S, const int& p, const Eigen::VectorXd hyp_ta, const bool& save_draws, const Rcpp::NumericVector& diff_deltas, const bool& verbose, const int& print_every) 
{
  // multiGGM: derived from upstream mGHS_chain(). The sampling steps are
  // unchanged, so for thin = 1 the draws match upstream for the same seed.
  // Changes: thinning; draws stored as an R array [p, p, K, nsave]; running
  // sums over every post-burn-in iteration (posterior means and group
  // differences); stuck detection after each three-parameter gamma draw;
  // user interrupts; R-style progress output.
  //
  // INPUT
  // burnin, nsave, thin: iterations are burnin + nsave * thin
  // n: K-dim vector containing sample size of each group
  // S: K-dim list with matrix S[[k]] = t(Yk) %*% Yk
  // p: dimensionality of the problem
  // hyp_ta: hyperparameters a and b of cut model
  // save_draws: store precision draws (C_save); running sums are always kept
  // diff_deltas: thresholds c for the running counts of |D| > c
  
  typedef array_type::index_range range;
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // initialization
  
  // const int p = S.size(), pm1 = p-1;
  const int pm1 = p-1;
  const int K = n.size(), Km1 = K-1;
  const int B = burnin + nsave * thin;
  const int bn = burnin;
  mghs_reject_tries = 0.0;
  
  int h = 0, hs = 0, b = 0, j = 0, js = 0, c = 0, cs = 0;
  
  // working objects
  const double IW_df = 0.5 * p * pm1;
  double gamma = 0.0;
  double alpha_l = 0.0, beta_l = 0.0;
  double alpha_t = 0.0, beta_t = 0.0;
  double s_22 = 0.0;
  double mth = 0.0, mth2 = 0.0;
  double tmp1 = 0.0, tmp2 = 0.0, sv = 0.0;
  double a_MH = 0.0;
  
  Eigen::VectorXd delta_12_inv(Km1); delta_12_inv.setZero();
  Eigen::VectorXd diag_V(K); diag_V.setZero();
  Eigen::VectorXd mu(K); mu.setOnes();
  Eigen::VectorXd e_ij(K); e_ij.setZero();
  Eigen::VectorXd v(pm1); v.setZero();
  Eigen::VectorXd temp_1(pm1); temp_1.setZero();
  Eigen::VectorXd temp_2(pm1); temp_2.setZero();
  Eigen::VectorXd temp_3(pm1); temp_3.setZero();
  Eigen::VectorXd s_12(pm1); s_12.setZero();
  Eigen::VectorXd sigma_12(pm1); sigma_12.setZero();
  Eigen::VectorXd delta_12(pm1); delta_12.setZero();
  
  Eigen::MatrixXd Rinv(Km1, Km1); Rinv.setZero();
  Eigen::MatrixXd rtR(Km1, K); rtR.setZero();
  Eigen::MatrixXd IW_scale(K, K); IW_scale.setZero();
  Eigen::MatrixXd Psi(K, K); Psi.setZero();
  Eigen::MatrixXd R_s(K, K); R_s.setZero();
  Eigen::MatrixXd inv_Omega_11(pm1, pm1); inv_Omega_11.setZero();
  
  // save values
  Eigen::VectorXd tau(K); tau.setOnes();
  Eigen::VectorXd zeta(K); zeta.setOnes();
  
  array_type Omega(boost::extents[p][p][K]);
  array_type Sigma(boost::extents[p][p][K]);
  for (h = 0; h < K; h++)
  {
    for (j = 0; j < p; j++)
    {
      Omega[j][j][h] = 1.0;
      Sigma[j][j][h] = 1.0;
    }
  }
  
  array_type Lambda(boost::extents[p][p][K]);
  std::fill(Lambda.origin(), Lambda.origin() + K*p*p, 1.0);
  
  array_type Eta(boost::extents[p][p][K]);
  std::fill(Eta.origin(), Eta.origin() + K*p*p, 1.0);
  
  Eigen::MatrixXd R(K, K); R.setIdentity();
  array_type S1(boost::extents[p][p][K]);
  for (h = 0; h < K; h++)
    S1[boost::indices[range()][range()][h]] = fill_S(S[h], p);
  
  // output
  

  double ta_s = 0.0, l2 = 0.0;
  
  Eigen::VectorXd l1(K); l1.setZero();
  Eigen::VectorXd qta2((int)IW_df); qta2.setZero();
  Eigen::VectorXd ta(K);
  
  for (h = 0; h < K; h++) 
    ta(h) = 0.5;
  
  Eigen::MatrixXd Z((int)IW_df, K); Z.setZero();
  Eigen::MatrixXd qta1((int)IW_df, K); qta1.setZero();
  
  for (h = 0; h < K; h++)
  {
    l1(h) = IW_df * std::log(0.5);
    for (j = 0; j < IW_df; j++)
    {
      qta1(j, h) = 0.5;
      Z(j, h) = R::rbinom(1, 0.5);
    }
  } 
  
  // multiGGM storage
  Rcpp::NumericVector C_save(save_draws ? p * p * K * nsave : 0);
  Rcpp::NumericVector R_save(K * K * nsave), tau_save(K * nsave), ta_save(K * nsave);
  Rcpp::NumericVector C_sum(p * p * K), pcor_sum(p * p * K), z_sum(p * p * K);
  Eigen::MatrixXd dsq(p, K);
  const int npe = p * pm1 / 2;
  std::vector<int> gp1, gp2;
  for (h = 0; h < K; h++)
    for (hs = h + 1; hs < K; hs++) { gp1.push_back(h); gp2.push_back(hs); }
  const int ngp = gp1.size();
  std::vector<double> deltas(diff_deltas.begin(), diff_deltas.end());
  DiffAcc acc_prec(npe, ngp, deltas), acc_pcor(npe, ngp, deltas);
  int n_post = 0, e = 0, g = 0, h1 = 0, h2 = 0, s_idx = 0;
  double d1 = 0.0, d2 = 0.0;
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // start algorithm
  // ::::::::::::::::::::::::::::::::::::::::::::
  
  for (b = 0; b < B; b++)
  {
    for (h = 0; h < K; h++)
    {
      for (j = 0; j < p; j++)
      {
        // ::::::::::::::::::::::::::::::::::::::::::::
        // sample gamma and v
        
        s_22 = S1[j][j][h];
        gamma = R::rgamma(0.5 * n(h) + 1.0, 2.0 / s_22);
        
        cs = 0;
        for (js = 0; js < p; js++)
        {
          if (js != j) 
          {
            c = 0;
            for (hs = 0; hs < K; hs++)
            {
              if (hs != h)
              {
                delta_12_inv(c) = Omega[js][j][hs] / std::sqrt(tau(hs) * Lambda[js][j][hs]);
                c += 1;
              }
            }
            
            tmp1 = tau(h) * Lambda[js][j][h];
            delta_12(cs) = tmp1;
            
            v(cs) = R::rnorm(0, 1);
            temp_1(cs) = (rtR.col(h).transpose() * delta_12_inv).value();
            temp_2(cs) = std::sqrt(tmp1) * temp_1(cs);
            
            s_12(cs) = S1[js][j][h];
            sigma_12(cs) = Sigma[js][j][h];
            
            cs += 1;
          }
        }
        
        delta_12 = inverse(delta_12, pm1) / mu(h);
        inv_Omega_11 = delete_rc_from_boost1(Sigma[boost::indices[range(0, p)][range(0, p)][h]], j, pm1);
        inv_Omega_11 -= sigma_12 * (sigma_12.transpose() / Sigma[j][j][h]);
        
        v = rmvnorm(v, j, s_22, s_12, inv_Omega_11, delta_12, temp_2, pm1);
        
        // ::::::::::::::::::::::::::::::::::::::::::::
        // sample lambda and eta
        
        mth = mu(h) * tau(h);
        mth2 = std::sqrt(tau(h)) * mu(h);
        
        temp_2 = inv_Omega_11 * v;
        temp_3 = temp_2 / gamma;
        
        c = 0;
        for (js = 0; js < p; js++)
        {
          if (js != j) 
          {
            tmp1 = v(c);
            Omega[js][j][h] = tmp1;
            Omega[j][js][h] = tmp1;
            
            alpha_l = std::sqrt(0.5 * tmp1 * tmp1 / mth + 1.0 / Eta[js][j][h]);
            beta_l = temp_1(c) * tmp1 / mth2;
            
            if (!std::isfinite(alpha_l) || !std::isfinite(beta_l) || alpha_l > 1e+06 || beta_l > 1e+06) {
              return stuck_result(b, (!std::isfinite(alpha_l) || alpha_l > 1e+06) ? 1 : 5);
            }
            
            tmp1 = rg3p_c1(alpha_l, beta_l);
            if (!(tmp1 > 0.0) || !std::isfinite(1.0 / (tmp1 * tmp1))) return stuck_result(b, 2);
            tmp1 = 1.0 / (tmp1 * tmp1);
            Lambda[js][j][h] = tmp1;  
            Lambda[j][js][h] = tmp1;  
            
            tmp1 = 1.0 / R::rgamma(1, 1.0 / (1.0 + 1.0 / tmp1));
            Eta[js][j][h] = tmp1;
            Eta[j][js][h] = tmp1;
            
            tmp1 = -temp_3(c);
            Sigma[js][j][h] = tmp1;
            Sigma[j][js][h] = tmp1;
            
            c += 1;
          }
        }
        
        inv_Omega_11 += temp_2 * temp_3.transpose();
        
        Omega[j][j][h] = gamma + v.transpose() * temp_2;
        Sigma[j][j][h] = 1.0 / gamma;
        
        Sigma[boost::indices[range(0, p)][range(0, p)][h]] = get_Sigma(Sigma[boost::indices[range(0, p)][range(0, p)][h]], inv_Omega_11, j, pm1, p);
      }
    }
    
    // ::::::::::::::::::::::::::::::::::::::::::::
    // sample tau and zeta
    
    diag_V.setZero();
    
    for (h = 0; h < K; h++)
    {
      alpha_t = 0.0;
      beta_t = 0.0;
      
      sv = std::pow(Omega[0][0][h], 2);
      // sv = 0.0;
      
      for (j = 1; j < p; j++)
      {
        sv += std::pow(Omega[j][j][h], 2);
        
        for (js = 0; js < j; js++)
        {
          tmp2 = Omega[js][j][h];
          sv += tmp2 * tmp2;   
          
          tmp1 = tmp2 / std::sqrt(Lambda[js][j][h]);
          alpha_t += tmp1 * tmp1;
          
          c = 0;
          for (hs = 0; hs < K; hs++)
          {
            if (hs != h)
            {
              delta_12_inv(c) = Omega[js][j][hs] / std::sqrt(tau(hs) * Lambda[js][j][hs]);
              c += 1;
            }
          }
          
          beta_t += tmp1 * (rtR.col(h).transpose() * delta_12_inv).value();
        }
      }
      
      diag_V(h) = 1.0 / std::sqrt(sv);
      alpha_t = std::sqrt(0.5 / mu(h) * alpha_t + 1.0 / zeta(h));
      
      if (!std::isfinite(alpha_t) || !std::isfinite(beta_t / mu(h)) || alpha_t > 1e+06 || beta_t / mu(h) > 1e+06) {
        return stuck_result(b, 3);
      }
      
      tau(h) = 1.0 / std::pow(rg3p_approx(alpha_t, beta_t / mu(h), (int)IW_df), 2);
      if (!(tau(h) > 0.0) || !std::isfinite(tau(h))) return stuck_result(b, 4);
      // tau(h) = 1.0 / std::pow(rg3p(alpha_t, beta_t / mu(h), (int)IW_df), 2);
      zeta(h) = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / tau(h)));
    }
    
    // ::::::::::::::::::::::::::::::::::::::::::::
    // sample R
    
    IW_scale.setZero();
    for (j = 1; j < p; j++)
    {
      for (js = 0; js < j; js++)
      {
        for (h = 0; h < K; h++)
          e_ij(h) = diag_V(h) * Omega[js][j][h] / std::sqrt(tau(h) * Lambda[js][j][h]);
        IW_scale += e_ij * e_ij.transpose();
      }
    }
    
    Psi = rIW(IW_df, IW_scale, K);
    e_ij = Psi.diagonal();
    e_ij = sinverse(e_ij, K);
    R_s = e_ij.asDiagonal() * Psi * e_ij.asDiagonal();
    
    a_MH = exp((K+1) * 0.5 * (std::log(R_s.determinant()) - std::log(R.determinant())));
    
    if (R::runif(0, 1) < a_MH)
    {
      R = R_s;
      
      for (h = 0; h < K; h++)
      {
        Rinv = delete_rc_from_eigen(R, h, Km1);
        Rinv = Rinv.inverse();
        rtR.col(h) = delete_elem_from_eigen_vec(R.col(h), h, Km1);
        rtR.col(h) = rtR.col(h).transpose() * Rinv;
        
        mu(h) = 1.0 - rtR.col(h).transpose() * delete_elem_from_eigen_vec(R.col(h), h, Km1);
      }
    }
    
    // ::::::::::::::::::::::::::::::::::::::::::::
    // sample z and t with cut-model
    
    for (h = 0; h < K; h++)
    {
      ta_s = R::rbeta(hyp_ta(0), hyp_ta(1));
      l1(h) = 0.0; l2 = 0.0;
      
      cs = 0;
      for (j = 0; j < p-1; j++)
      {
        for (js = j+1; js < p; js++)
        {
          tmp1 = Lambda[js][j][h];
          qta1(cs, h) = R::pnorm(tmp1 / (1 + tmp1) - ta(h), 0.0, 1.0, true, false);
          qta2(cs) = R::pnorm(tmp1 / (1 + tmp1) - ta_s, 0.0, 1.0, true, false);
          
          tmp1 = qta1(cs, h);
          Z(cs, h) = R::rbinom(1, tmp1);
          
          if (Z(cs, h) == 1)
          {
            l1(h) += std::log(tmp1);
            l2 += std::log(qta2(cs));
          } else 
          {
            l1(h) += std::log(1.0 - tmp1);
            l2 += std::log(1.0 - qta2(cs));
          }
          cs += 1;
        }
      }
      
      a_MH = l2 + R::dbeta(ta_s, hyp_ta(0), hyp_ta(1), true) - l1(h) - R::dbeta(ta(h), hyp_ta(0), hyp_ta(1), true);
      if (std::log(R::runif(0, 1)) < a_MH)
      {
        ta(h) = ta_s;
        l1(h) = l2;
        qta1.col(h) = qta2;
      }
    }
    
    // ::::::::::::::::::::::::::::::::::::::::::::
    // save values
    
    
    // l1_out.col(b) = l1;
    
    if (b >= bn)
    {
      n_post += 1;
      for (h = 0; h < K; h++)
        for (j = 0; j < p; j++)
          dsq(j, h) = std::sqrt(Omega[j][j][h]);
      
      for (h = 0; h < K; h++)
      {
        for (j = 0; j < p; j++)
        {
          for (js = 0; js < p; js++)
          {
            c = js + p * j + p * p * h;
            C_sum[c] += Omega[js][j][h];
            pcor_sum[c] += (js == j) ? 1.0 : -Omega[js][j][h] / (dsq(js, h) * dsq(j, h));
          }
        }
        // Z rows are ordered row-wise: (0,1), (0,2), ..., (0,p-1), (1,2), ...
        c = 0;
        for (j = 0; j < pm1; j++)
        {
          for (js = j + 1; js < p; js++)
          {
            z_sum[j + p * js + p * p * h] += Z(c, h);
            z_sum[js + p * j + p * p * h] += Z(c, h);
            c += 1;
          }
        }
      }
      
      // group differences, edges in R's upper.tri() order
      for (g = 0; g < ngp; g++)
      {
        h1 = gp1[g]; h2 = gp2[g];
        e = 0;
        for (j = 1; j < p; j++)
        {
          for (js = 0; js < j; js++)
          {
            d1 = Omega[js][j][h1];
            d2 = Omega[js][j][h2];
            acc_prec.add(g, e, d1 - d2);
            acc_pcor.add(g, e, -d1 / (dsq(js, h1) * dsq(j, h1)) + d2 / (dsq(js, h2) * dsq(j, h2)));
            e += 1;
          }
        }
      }
      
      if ((b - bn) % thin == 0)
      {
        s_idx = (b - bn) / thin;
        for (h = 0; h < K; h++)
        {
          tau_save[h + K * s_idx] = tau(h);
          ta_save[h + K * s_idx] = ta(h);
          for (hs = 0; hs < K; hs++)
            R_save[h + K * hs + K * K * s_idx] = R(h, hs);
          if (save_draws)
            for (j = 0; j < p; j++)
              for (js = 0; js < p; js++)
                C_save[js + p * j + p * p * h + p * p * K * s_idx] = Omega[js][j][h];
        }
      }
    }
    
    Rcpp::checkUserInterrupt();
    if (verbose && print_every > 0 && (b + 1) % print_every == 0)
      Rcpp::Rcout << "mghs iteration " << (b + 1) << " / " << B << "\n";
  }
  
  // ::::::::::::::::::::::::::::::::::::::::::::
  // output
  // ::::::::::::::::::::::::::::::::::::::::::::
  
  for (c = 0; c < p * p * K; c++)
  {
    C_sum[c] /= n_post;
    pcor_sum[c] /= n_post;
    z_sum[c] /= n_post;
  }
  Rcpp::IntegerVector dim3 = Rcpp::IntegerVector::create(p, p, K);
  C_sum.attr("dim") = dim3;
  pcor_sum.attr("dim") = dim3;
  z_sum.attr("dim") = dim3;
  R_save.attr("dim") = Rcpp::IntegerVector::create(K, K, nsave);
  tau_save.attr("dim") = Rcpp::IntegerVector::create(K, nsave);
  ta_save.attr("dim") = Rcpp::IntegerVector::create(K, nsave);
  Rcpp::RObject C_out = R_NilValue;
  if (save_draws)
  {
    C_save.attr("dim") = Rcpp::IntegerVector::create(p, p, K, nsave);
    C_out = C_save;
  }
  Rcpp::IntegerMatrix pairs(ngp, 2);
  for (g = 0; g < ngp; g++) { pairs(g, 0) = gp1[g] + 1; pairs(g, 1) = gp2[g] + 1; }
  
  return Rcpp::List::create(
    Rcpp::Named("stuck") = 0,
    Rcpp::Named("n_post") = n_post,
    Rcpp::Named("reject_tries") = mghs_reject_tries,
    Rcpp::Named("C_save") = C_out,
    Rcpp::Named("R_save") = R_save,
    Rcpp::Named("tau_save") = tau_save,
    Rcpp::Named("ta_save") = ta_save,
    Rcpp::Named("C_mean") = C_sum,
    Rcpp::Named("pcor_mean") = pcor_sum,
    Rcpp::Named("cut_z_mean") = z_sum,
    Rcpp::Named("diff_acc") = Rcpp::List::create(
      Rcpp::Named("precision") = acc_prec.to_list(),
      Rcpp::Named("pcor") = acc_pcor.to_list(),
      Rcpp::Named("deltas") = diff_deltas,
      Rcpp::Named("group_pairs") = pairs));
}

// end file
