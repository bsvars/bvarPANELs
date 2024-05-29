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
    const int                     thin = 100, // introduce thinning
    const bool                    show_progress = true
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
  
  int   ss = 0;
  
  for (int s=0; s<S; s++) {
    // Rcout << "Iteration: " << s << endl;
    
    // Increment progress bar
    if (any(prog_rep_points == s)) p.increment();
    // Check for user interrupts
    if (s % 200 == 0) checkUserInterrupt();
    
    // sample aux_m, aux_w, aux_s
    // Rcout << "  sample m" << endl;
    aux_m       = sample_m( aux_A, aux_V, aux_s, aux_w, prior );
    
    // Rcout << "  sample w" << endl;
    aux_w       = sample_w( aux_V, prior );
    
    // Rcout << "  sample s" << endl;
    aux_s       = sample_s( aux_A, aux_V, aux_Sigma, aux_m, prior );
    
    // sample aux_nu
    // Rcout << "  sample nu" << endl;
    aux_nu      = sample_nu( aux_nu, aux_Sigma_c, aux_Sigma_c_inv, aux_Sigma, prior );
    
    // sample aux_Sigma
    // Rcout << "  sample Sigma" << endl;
    aux_Sigma   = sample_Sigma( aux_Sigma_c_inv, aux_s, aux_nu, prior );
    
    // sample aux_A, aux_V
    // Rcout << "  sample AV" << endl;
    field<mat> tmp_AV     = sample_AV( aux_A_c, aux_Sigma_c_inv, aux_s, aux_m, aux_w, prior );
    aux_A       = tmp_AV(0);  
    aux_V       = tmp_AV(1);
    
    // sample aux_A_c, aux_Sigma_c
    // Rcout << "  sample A_c Sigma_c" << endl;
    // Rcout << "  aux_nu: " << aux_nu << endl;
    for (int c=0; c<C; c++) {
      field<mat> tmp_A_c_Sigma_c  = sample_A_c_Sigma_c( y(c), x(c), aux_A, aux_V, aux_Sigma, aux_nu );
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
