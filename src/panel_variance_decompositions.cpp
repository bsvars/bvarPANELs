
#include <RcppArmadillo.h>
#include <bsvars.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp:interface(cpp,r)]]
// [[Rcpp::export]]
arma::cube Sigma2B_c (
    arma::cube&   posterior_Sigma_c,    // (N, N, S)
    const bool    lower = true
) {
  
  const int   S = posterior_Sigma_c.n_slices;
  const int   N = posterior_Sigma_c.n_rows;
  
  cube        posterior_B_c(N, N, S);
  
  for (int s = 0; s < S; s++) {
    if ( lower ) {
      posterior_B_c.slice(s) = chol(posterior_Sigma_c.slice(s), "lower");
    } else {
      posterior_B_c.slice(s) = chol(posterior_Sigma_c.slice(s), "upper");
    }
  } // END s loop
    
  return posterior_B_c;
} // END Sigma2B_c



arma::cube flip_cube_rows_cols (
  arma::cube&   x   // (N, K, S)
) {
  
  const int   N = x.n_rows;
  const int   K = x.n_cols;
  const int   S = x.n_slices;
  
  cube        y(K, N, S);
  
  for (int s = 0; s < S; s++) {
    y.slice(s)  = trans(x.slice(s));
  }
  return y;
} // END flip_cube_rows_cols



// [[Rcpp:interface(cpp,r)]]
// [[Rcpp::export]]
arma::field<arma::cube> panel_variance_decompositions (
    arma::field<arma::cube>&  posterior_Sigma,    // (S)(N, N, C)
    arma::field<arma::cube>&  posterior_A,        // (S)(K, N, C)
    arma::cube&               global_Sigma,       // (N, N, S)
    arma::cube&               global_A,           // (K, N, S)
    const int     horizon,
    const int     p,
    const bool    lower = true
) {
  
  const int   S = posterior_Sigma.n_rows;
  const int   C = posterior_Sigma(0).n_slices;
  
  field<cube> fevds(C + 1, S);
  
  for (int s=0; s<S; s++) {
    cube        posterior_B_s     = Sigma2B_c ( posterior_Sigma(s), lower );
    
    cube A_s_aperm = flip_cube_rows_cols ( posterior_A(s) );
    field<cube> posterior_irf_s   = bsvars::bsvars_ir (
                                      posterior_B_s,
                                      A_s_aperm,
                                      horizon,
                                      p,
                                      false
                                    );
    
    field<cube> posterior_fevd_s  = bsvars::bsvars_fevd_homosk ( posterior_irf_s );
    for (int c = 0; c < C; c++) {
      fevds(c,s) = posterior_fevd_s(c);
    } // END s loop
  } // END c loop
  
  // Global fevds
  cube        global_B    = Sigma2B_c ( global_Sigma, lower );
  cube        global_A_aperm = flip_cube_rows_cols ( global_A );
  field<cube> global_irf  = bsvars::bsvars_ir (
                              global_B,
                              global_A_aperm,
                              horizon,
                              p,
                              true
                            );
  field<cube> global_fevd = bsvars::bsvars_fevd_homosk ( global_irf );
  
  for (int s=0; s<S; s++) {
    fevds(C,s) = global_fevd(s);
  } // END s loop
  
  return fevds;
} // END panel_variance_decompositions
