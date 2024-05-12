
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
