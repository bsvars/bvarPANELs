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
