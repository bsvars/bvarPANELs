
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


// [[Rcpp:interface(cpp,r)]]
// [[Rcpp::export]]
arma::field<arma::cube> panel_variance_decompositions (
    arma::field<arma::cube>&  posterior_Sigma,    // (C)(N, N, S)
    arma::field<arma::cube>&  posterior_A,        // (C)(K, N, S)
    arma::cube&               global_Sigma,       // (N, N, S)
    arma::cube&               global_A,           // (K, N, S)
    const int     horizon,
    const int     p,
    const bool    lower = true
) {
  
  const int   C = posterior_Sigma.n_elem;
  const int   S = posterior_Sigma(0).n_slices;
  
  field<cube> fevds(C + 1,S);
  
  for (int c = 0; c < C; c++) {
    
    cube        posterior_B_c     = Sigma2B_c ( posterior_Sigma(c), lower );
    field<cube> posterior_irf_c   = bsvars::bsvars_ir (
                                      posterior_B_c,
                                      posterior_A(c),
                                      horizon,
                                      p,
                                      true
                                    );
    field<cube> posterior_fevd_c  = bsvars::bsvars_fevd_homosk ( posterior_irf_c );
    
    for (int s=0; s<S; s++) {
      fevds(c,s) = posterior_fevd_c(s);
    } // END s loop
  } // END c loop
  
  // Global fevds
  cube        global_B    = Sigma2B_c ( global_Sigma, lower );
  field<cube> global_irf  = bsvars::bsvars_ir (
                              global_B,
                              global_A,
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
