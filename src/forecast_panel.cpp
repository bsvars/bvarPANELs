
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List forecast_bvarPANEL (
    arma::field<arma::cube>&  posterior_A_c_cpp,      // (S)(K, N, C)
    arma::field<arma::cube>&  posterior_Sigma_c_cpp,  // (S)(N, N, C)
    Rcpp::List&               X_c,                    // (C)(T_c, K)
    const int                 horizon
) {
  
  const int       N = posterior_A_c_cpp(0).n_cols;
  const int       S = posterior_A_c_cpp.n_elem;
  const int       K = posterior_A_c_cpp(0).n_rows;
  const int       C = posterior_A_c_cpp(0).n_slices;
  const int       p = (K - 1) / N;
  
  field<cube>     forecasts(C);                       // of (horizon, N, S) cubes
  rowvec          one(1, fill::value(1));
  
  for (int c=0; c<C; c++) {
    
    mat     XXcc    = as<mat>(X_c[c]);
    rowvec  x_t     = XXcc.tail_rows(1).cols(0, K - 2);
    cube    forecasts_c(horizon, N, S);
    
    for (int s=0; s<S; s++) {
      
      vec Xt        = trans(join_rows(x_t, one));
      mat Sigma_cs  = posterior_Sigma_c_cpp(s).slice(c);
      mat A_cs      = trans(posterior_A_c_cpp(s).slice(c));
    
      for (int h=0; h<horizon; h++) {
        
        forecasts_c.slice(s).row(h) = trans(mvnrnd( A_cs * Xt, Sigma_cs ));
        
        if ( p == 1 ) {
          Xt          = trans(join_rows(forecasts_c.slice(s).row(h), one));
        } else {
          Xt          = trans(join_rows(forecasts_c.slice(s).row(h), Xt.subvec(N, K - 2).t(), one));
        }
        
      } // END h loop
    } // END s loop
    
    forecasts(c) = forecasts_c;
    
  } // END c loop
  
  return List::create(
    _["forecasts_cpp"]  = forecasts
  );
} // END forecast_bsvarPANEL



// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec mvnrnd_cond (
    arma::vec x,        // Nx1 with NAs or without
    arma::vec mu,       // Nx1 mean vector
    arma::mat Sigma     // NxN covariance matrix
) {
  int   N         = x.n_elem;
  uvec  ind       = find_finite(x);
  uvec  ind_nan   = find_nan(x);
  mat   aj        = eye(N, N);
  
  vec   x2        = x(ind); 
  
  vec   mu1       = mu(ind_nan);
  vec   mu2       = mu(ind);
  mat   Sigma11   = Sigma(ind_nan, ind_nan);
  mat   Sigma12   = Sigma(ind_nan, ind);
  mat   Sigma22   = Sigma(ind, ind);
  mat   Sigma22_inv = inv_sympd(Sigma22);
  
  vec   mu_cond     = mu1 + Sigma12 * Sigma22_inv * (x2 - mu2);
  mat   Sigma_cond  = Sigma11 - Sigma12 * Sigma22_inv * Sigma12.t();
  
  vec   draw = mvnrnd( mu_cond, Sigma_cond);
  
  vec   out = aj.cols(ind_nan) * draw + aj.cols(ind) * x2;
  return out;
} // END mvnrnd_cond




// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
Rcpp::List forecast_conditional_bvarPANEL (
    arma::field<arma::cube>&  posterior_A_c_cpp,      // (S)(K, N, C)
    arma::field<arma::cube>&  posterior_Sigma_c_cpp,  // (S)(N, N, C)
    Rcpp::List&               X_c,                    // (C)(T_c, K)
    Rcpp::List&               cond_forecasts,         // (C)(horizon, N)
    const int                 horizon
) {
  
  const int       N = posterior_A_c_cpp(0).n_cols;
  const int       S = posterior_A_c_cpp.n_elem;
  const int       K = posterior_A_c_cpp(0).n_rows;
  const int       C = posterior_A_c_cpp(0).n_slices;
  const int       p = (K - 1) / N;
  
  field<cube>     forecasts(C);                       // of (horizon, N, S) cubes
  rowvec          one(1, fill::value(1));
  
  for (int c=0; c<C; c++) {
    
    mat     XXcc    = as<mat>(X_c[c]);
    mat     cond_fc = as<mat>(cond_forecasts[c]);
    rowvec  x_t     = XXcc.tail_rows(1).cols(0, K - 2);
    cube    forecasts_c(horizon, N, S);
    
    for (int s=0; s<S; s++) {
      
      vec Xt        = trans(join_rows(x_t, one));
      mat Sigma_cs  = posterior_Sigma_c_cpp(s).slice(c);
      mat A_cs      = trans(posterior_A_c_cpp(s).slice(c));
      
      for (int h=0; h<horizon; h++) {
        
        forecasts_c.slice(s).row(h) = trans(
          mvnrnd_cond( cond_fc.row(h).t(), A_cs * Xt, Sigma_cs )
        );
        
        if ( p == 1 ) {
          Xt          = trans(join_rows(forecasts_c.slice(s).row(h), one));
        } else {
          Xt          = trans(join_rows(forecasts_c.slice(s).row(h), Xt.subvec(N, K - 2).t(), one));
        }
        
      } // END h loop
    } // END s loop
    
    forecasts(c) = forecasts_c;
    
  } // END c loop
  
  return List::create(
    _["forecasts_cpp"]  = forecasts
  );
} // END forecast_conditional_bvarPANEL