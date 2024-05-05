#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;


// [[Rcpp:interface(cpp,r)]]
// [[Rcpp::export]]
arma::field<arma::mat> rmniw1(
    const arma::mat& A,     // KxN
    const arma::mat& V,     // KxK
    const arma::mat& S,     // NxN
    const double&    nu     // scalar
) {
  mat Sigma     = iwishrnd(S, nu);
  mat X_tmp     = mat(size(A), fill::randn);
  mat X         = A + chol(V).t() * X_tmp * chol(Sigma);
  
  field<mat> out(2);
  out(0)        = X;
  out(1)        = Sigma;
  return out;
} // END rmniw1


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double sample_m (
    const arma::mat&    aux_A,    // KxN
    const arma::mat&    aux_V,    // KxK
    const double&       aux_s,   // scalar
    const double&       aux_w,   // scalar
    const Rcpp::List&   prior
) {
  
  const int N       = aux_A.n_cols;
  mat prior_S_inv   = as<mat>(prior["S_inv"]);
  mat prior_S       = inv_sympd(prior_S_inv);
  double prior_mu_m = as<double>(prior["mu_m"]);
  double prior_sigma2_m = as<double>(prior["sigma2_m"]);
  
  double precision_tmp  = (1 / prior_sigma2_m);
  double mean_tmp       = (prior_mu_m / prior_sigma2_m);
  for (int n = 0; n < N; n++) {
    precision_tmp  += 1 / (aux_s * aux_w * prior_S(n, n) * aux_V(n, n));
    mean_tmp       += aux_A(n, n) / (aux_s * aux_w * prior_S(n, n) * aux_V(n, n));
  }
  double sigma2_m_bar   = 1 / precision_tmp;
  double mu_m_bar       = sigma2_m_bar * mean_tmp;
  
  double out        = randn( distr_param(mu_m_bar, pow(sigma2_m_bar, 0.5)) );
  return out;
} // END sample_m


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double sample_w (
    const arma::mat&    aux_V,    // KxK
    const Rcpp::List&   prior
) {
  
  int K             = aux_V.n_cols;
  mat prior_W       = as<mat>(prior["W"]);
  double prior_s_w  = as<double>(prior["s_w"]);
  double prior_a_w  = as<double>(prior["a_w"]);
  double prior_eta  = as<double>(prior["eta"]);
  
  mat aux_V_inv     = inv_sympd(aux_V);
  double s_w_bar    = prior_s_w + 0.5 * trace(aux_V_inv * prior_W);
  double a_w_bar    = prior_a_w + 0.5 * K * prior_eta;
  
  double out        = randg( distr_param(a_w_bar, s_w_bar) );
  return out;
} // END sample_w


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
double sample_s (
    const arma::mat&    aux_A,      // KxN
    const arma::mat&    aux_V,      // KxK
    const arma::mat&    aux_Sigma,  // NxN
    const double&       aux_m,      // scalar
    const Rcpp::List&   prior
) {
  
  int N             = aux_A.n_cols;
  int K             = aux_A.n_rows;
  mat prior_M       = as<mat>(prior["M"]);
  mat prior_S_inv   = as<mat>(prior["S_inv"]);
  mat prior_S_Sigma_inv = as<mat>(prior["S_Sigma_inv"]);
  double prior_s_s  = as<double>(prior["s_s"]);
  double prior_nu_s = as<double>(prior["nu_s"]);
  double prior_mu_Sigma = as<double>(prior["mu_Sigma"]);
  
  mat quad_tmp1     = (aux_A - aux_m * prior_M) * prior_S_inv * trans(aux_A - aux_m * prior_M);
  
  double s_s_bar    = prior_s_s + trace(solve(aux_V, quad_tmp1)) + trace(prior_S_Sigma_inv * aux_Sigma);
  double nu_s_bar   = prior_nu_s + K * N + N * prior_mu_Sigma;
  
  double out        = s_s_bar / chi2rnd( nu_s_bar );
  return out;
} // END sample_s


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat sample_Sigma (
    const arma::cube&   aux_Sigma_c_inv,  // NxNxC
    const double&       aux_s,            // scalar
    const double&       aux_nu,           // scalar
    const Rcpp::List&   prior
) {
  
  int C             = aux_Sigma_c_inv.n_slices;
  int N             = aux_Sigma_c_inv.n_rows;
  mat prior_S_Sigma_inv = as<mat>(prior["S_Sigma_inv"]);
  double prior_mu_Sigma = as<double>(prior["mu_Sigma"]);
  
  mat sum_aux_Sigma_c_inv(N, N);
  for (int c = 0; c < C; c++) {
    sum_aux_Sigma_c_inv += aux_Sigma_c_inv.slice(c);
  }
  
  mat S_Sigma_bar   = (prior_S_Sigma_inv / aux_s) + sum_aux_Sigma_c_inv;
  double mu_bar     = C * aux_nu + prior_mu_Sigma;
  
  mat out           = wishrnd( S_Sigma_bar, mu_bar );
  return out;
} // END sample_Sigma


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> sample_AV (
    const arma::cube&   aux_A_c,          // KxNxC
    const arma::cube&   aux_Sigma_c_inv,  // NxNxC
    const double&       aux_s,            // scalar
    const double&       aux_m,            // scalar
    const double&       aux_w,            // scalar
    const Rcpp::List&   prior
) {
  
  int C             = aux_A_c.n_slices;
  int N             = aux_A_c.n_cols;
  int K             = aux_A_c.n_rows;
  
  mat prior_S_inv   = as<mat>(prior["S_inv"]);
  mat prior_M       = as<mat>(prior["M"]);
  mat prior_W       = as<mat>(prior["W"]);
  double prior_eta  = as<double>(prior["eta"]);
  
  mat sum_Sc_inv(N, N);
  mat sum_Sc_invAt(N, K);
  mat sum_ASc_invAt(K, K);
  for (int c = 0; c < C; c++) {
    sum_Sc_inv     += aux_Sigma_c_inv.slice(c);
    mat Sc_invAt    = aux_Sigma_c_inv.slice(c) * aux_A_c.slice(c).t();
    sum_Sc_invAt   += Sc_invAt;
    sum_ASc_invAt  += aux_A_c.slice(c) * Sc_invAt;
  } // END c loop
  
  mat S_bar_inv     = (prior_S_inv / aux_s) + sum_Sc_inv;
  mat S_bar         = inv_sympd(S_bar_inv);
  mat M_bar_trans   = S_bar * ( (aux_m / aux_s) * (prior_S_inv * prior_M.t()) + sum_Sc_invAt);
  mat W_bar         = (aux_w * prior_W) + (pow(aux_m, 2) / aux_s) * (prior_M * prior_S_inv * prior_M.t())
    + sum_ASc_invAt - M_bar_trans.t() * S_bar_inv * M_bar_trans;
  double eta_bar    = C * N + prior_eta;  
  
  arma::field<arma::mat> aux_AV = rmniw1( M_bar_trans, S_bar, W_bar, eta_bar );
  aux_AV(0)         = trans(aux_AV(0));
  return aux_AV;
} // END sample_AV



// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::field<arma::mat> sample_A_c_Sigma_c (
    const arma::mat&    Y_c,              // T_cxN
    const arma::mat&    X_c,              // T_cxK
    const arma::mat&    aux_A,            // KxN
    const arma::mat&    aux_V,            // KxK
    const arma::mat&    aux_Sigma,        // NxN
    const double&       aux_nu            // scalar
) {
  int T_c           = Y_c.n_rows;
  
  mat aux_V_inv     = inv_sympd( aux_V );
  mat V_bar_inv     = X_c.t() * X_c + aux_V_inv;
  mat V_bar         = inv_sympd( V_bar_inv );
  mat A_bar         = V_bar * ( X_c.t() * Y_c + aux_V_inv * aux_A );
  mat Sigma_bar     = aux_Sigma + Y_c.t() * Y_c + aux_A.t() * aux_V_inv * aux_A - A_bar.t() * V_bar_inv * A_bar;
  double nu_bar     = T_c + aux_nu;
  
  arma::field<arma::mat> aux_A_c_Sigma_c = rmniw1( A_bar, V_bar, Sigma_bar, nu_bar );
  return aux_A_c_Sigma_c;
} // END sample_A_c_Sigma_c
