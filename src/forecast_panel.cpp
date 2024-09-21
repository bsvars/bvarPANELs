
#include <RcppArmadillo.h>
#include <bsvars.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
arma::vec mvnrnd_cond_rest (
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
Rcpp::List forecast_bvarPANEL (
    arma::field<arma::cube>&  posterior_A_c_cpp,      // (S)(K, N, C)
    arma::field<arma::cube>&  posterior_Sigma_c_cpp,  // (S)(N, N, C)
    Rcpp::List&               X_c,                    // (C)(T_c, K)
    Rcpp::List&               cond_forecasts,         // (C)(horizon, N)
    Rcpp::List&               exog_forecasts,         // (C)(horizon, d)
    const int                 horizon
) {
  
  const int       S = posterior_A_c_cpp.n_elem;
  const int       N = posterior_A_c_cpp(0).n_cols;
  const int       K = posterior_A_c_cpp(0).n_rows;
  const int       C = posterior_A_c_cpp(0).n_slices;
  
  mat     EXcc    = as<mat>(exog_forecasts[0]);
  const int       d = EXcc.n_cols;
  
  field<cube>     forecasts(C);                       // of (horizon, N, S) cubes
  
  for (int c=0; c<C; c++) {
    
    mat     XXcc    = as<mat>(X_c[c]);
    mat     EXcc    = as<mat>(exog_forecasts[c]);
    bool    do_exog = EXcc.is_finite();
    mat     cond_fc = as<mat>(cond_forecasts[c]);
    
    rowvec  x_t;
    if ( do_exog ) {
      x_t = XXcc.tail_rows(1).cols(0, K - 1 - d);
    } else {
      x_t = XXcc.tail_rows(1).cols(0, K - 1);
    }
    
    vec     Xt(K);
    cube    forecasts_c(horizon, N, S);
    
    for (int s=0; s<S; s++) {
      
      if ( do_exog ) {
        Xt          = trans(join_rows(x_t, EXcc.row(0)));
      } else {
        Xt          = trans(x_t);
      }
      
      mat Sigma_cs  = posterior_Sigma_c_cpp(s).slice(c);
      mat A_cs      = trans(posterior_A_c_cpp(s).slice(c));
      
      for (int h=0; h<horizon; h++) {
        
        vec   cond_fc_h   = trans(cond_fc.row(h));
        uvec  nonf_el     = find_nonfinite(cond_fc_h);
        int   nonf_no     = nonf_el.n_elem;
        
        if ( nonf_no == N ) {
          forecasts_c.slice(s).row(h) = trans(
            mvnrnd( A_cs * Xt, Sigma_cs )
          );
        } else {
          forecasts_c.slice(s).row(h) = trans(
            bsvars::mvnrnd_cond( cond_fc_h, A_cs * Xt, Sigma_cs )   // does not work if cond_fc_h is all nan
          );
        }
        
        if ( h != horizon - 1 ) {
          if ( do_exog ) {
            Xt          = trans(join_rows(forecasts_c.slice(s).row(h), Xt.subvec(N, K - 1 - d).t(), EXcc.row(h + 1)));
          } else {
            Xt          = trans(join_rows(forecasts_c.slice(s).row(h), Xt.subvec(N, K - 1).t()));
          }
        }
        
      } // END h loop
    } // END s loop
    
    forecasts(c) = forecasts_c;
    
  } // END c loop
  
  return List::create(
    _["forecasts_cpp"]  = forecasts
  );
} // END forecast_conditional_bvarPANEL