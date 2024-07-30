#include <RcppArmadillo.h>
#include "progress.hpp"
#include "bsvars.h"

#include "sample_mniw.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp:interface(cpp,r)]]
// [[Rcpp::export]]
Rcpp::List bvarPANEL(
    const int&                    S,          // No. of posterior draws
    const Rcpp::List&             Y,          // a C-list of T_cxN elements
    const Rcpp::List&             X,          // a C-list of T_cxK elements
    const Rcpp::List&             prior,      // a list of priors
    const Rcpp::List&             starting_values, 
    const int                     thin, // introduce thinning
    const bool                    show_progress,
    const arma::vec&              adptive_alpha_gamma, // 2x1 vector with target acceptance rate and step size
    bool                          simple = false
) {
  
  // Progress bar setup
  vec prog_rep_points = arma::round(arma::linspace(0, S, 50));
  
  std::string oo      = "";
  if ( thin != 1 ) {
    oo                = bsvars::ordinal(thin) + " ";
  }
  
  if (show_progress) {
    Rcout << "**************************************************|" << endl;
    Rcout << "bvarPANELs: Forecasting with Bayesian Hierarchical|" << endl;
    Rcout << "            Panel Vector Autoregressions          |" << endl;
    Rcout << "**************************************************|" << endl;
    // Rcout << " Gibbs sampler for the SVAR-SV model              |" << endl;
    // Rcout << "**************************************************|" << endl;
    Rcout << " Progress of the MCMC simulation for " << S << " draws" << endl;
    Rcout << "    Every " << oo << "draw is saved via MCMC thinning" << endl;
    Rcout << " Press Esc to interrupt the computations" << endl;
    Rcout << "**************************************************|" << endl;
  }
  Progress p(50, show_progress);
  
  cube    aux_A_c     = as<cube>(starting_values["A_c"]);
  cube    aux_Sigma_c = as<cube>(starting_values["Sigma_c"]);
  mat     aux_A       = as<mat>(starting_values["A"]);
  mat     aux_V       = as<mat>(starting_values["V"]);
  mat     aux_Sigma   = as<mat>(starting_values["Sigma"]);
  double  aux_nu      = as<double>(starting_values["nu"]);
  double  aux_m       = as<double>(starting_values["m"]);
  double  aux_w       = as<double>(starting_values["w"]);
  double  aux_s       = as<double>(starting_values["s"]);
  
  const int C         = aux_A_c.n_slices;
  const int N         = aux_A.n_cols;
  const int K         = aux_A.n_rows;
  
  const int   SS    = floor(S / thin);
  
  field<cube> posterior_A_c_cpp(SS);
  field<cube> posterior_Sigma_c_cpp(SS);
  cube        posterior_A(K, N, SS);
  cube        posterior_V(K, K, SS);
  cube        posterior_Sigma(N, N, SS);
  vec         posterior_nu(SS);
  vec         posterior_m(SS);
  vec         posterior_w(SS);
  vec         posterior_s(SS);
  
  field<mat> y(C);
  field<mat> x(C);
  cube aux_Sigma_c_inv(N, N, C);
  for (int c=0; c<C; c++) {
    y(c)                  = as<mat>(Y[c]);
    x(c)                  = as<mat>(X[c]);
    aux_Sigma_c_inv.slice(c) = inv_sympd( aux_Sigma_c.slice(c) );
  } // END c loop
  
  mat   aux_V_which = aux_V;
  int   ss = 0;
  
  // the initial value for the adaptive_scale is set to the negative inverse of 
  // Hessian for the posterior log_kenel for nu
  double adaptive_scale = cov_nu(aux_nu, C, N);
  vec aux_nu_tmp(2);

  for (int s=0; s<S; s++) {
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    // sample aux_m, aux_w, aux_s - this is OK
    aux_m       = sample_m( aux_A, aux_V, aux_s, aux_w, prior );
    
    if ( !simple ) {
      aux_w       = sample_w( aux_V, prior );
      aux_V_which = aux_V;
    } else {
      aux_w       = sample_w_simple( aux_A_c, aux_Sigma_c_inv, aux_A, aux_V, aux_m, aux_s, prior );
      aux_V_which = aux_w * aux_V;
    }
    
    aux_s       = sample_s( aux_A, aux_V_which, aux_Sigma, aux_m, prior );
    
    // sample aux_nu
    if ( !simple ) {
      aux_nu_tmp      = sample_nu ( aux_nu, adaptive_scale, aux_Sigma_c, aux_Sigma_c_inv, aux_Sigma, prior, s, adptive_alpha_gamma );
      aux_nu          = aux_nu_tmp(0);
      adaptive_scale  = aux_nu_tmp(1);
    }
    
    // sample aux_Sigma - this is OK!
    aux_Sigma   = sample_Sigma( aux_Sigma_c_inv, aux_s, aux_nu, prior );
    
    // sample aux_A, aux_V
    field<mat> tmp_AV     = sample_AV( aux_A_c, aux_Sigma_c_inv, aux_s, aux_m, aux_w, prior, simple );
    aux_A       = tmp_AV(0);  
    if ( !simple ) {
      aux_V       = tmp_AV(1);
      aux_V_which = aux_V;
    }
    
    // sample aux_A_c, aux_Sigma_c
    for (int c=0; c<C; c++) {
      field<mat> tmp_A_c_Sigma_c  = sample_A_c_Sigma_c( y(c), x(c), aux_A, aux_V_which, aux_Sigma, aux_nu );
      aux_A_c.slice(c)            = tmp_A_c_Sigma_c(0);
      aux_Sigma_c.slice(c)        = tmp_A_c_Sigma_c(1);
      aux_Sigma_c_inv.slice(c)    = inv_sympd( aux_Sigma_c.slice(c) );
    } // END c loop
    
    if (s % thin == 0) {
      posterior_A_c_cpp(ss)     = aux_A_c;
      posterior_Sigma_c_cpp(ss) = aux_Sigma_c;
      posterior_A.slice(ss)     = aux_A;
      posterior_V.slice(ss)     = aux_V;
      posterior_Sigma.slice(ss) = aux_Sigma;
      posterior_nu(ss)          = aux_nu;
      posterior_m(ss)           = aux_m;
      posterior_w(ss)           = aux_w;
      posterior_s(ss)           = aux_s;
      
      ss++;
    }
  } // END s loop
  
  return List::create(
    _["last_draw"]  = List::create(
      _["A_c"]      = aux_A_c,
      _["Sigma_c"]  = aux_Sigma_c,
      _["A"]        = aux_A,
      _["V"]        = aux_V,
      _["Sigma"]    = aux_Sigma,
      _["nu"]       = aux_nu,
      _["m"]        = aux_m,
      _["w"]        = aux_w,
      _["s"]        = aux_s
    ),
    _["posterior"]  = List::create(
      _["A_c_cpp"]  = posterior_A_c_cpp,
      _["Sigma_c_cpp"]  = posterior_Sigma_c_cpp,
      _["A"]        = posterior_A,
      _["V"]        = posterior_V,
      _["Sigma"]    = posterior_Sigma,
      _["nu"]       = posterior_nu,
      _["m"]        = posterior_m,
      _["w"]        = posterior_w,
      _["s"]        = posterior_s
    )
  );
} // END bvarPANEL
