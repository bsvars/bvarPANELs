#ifndef _FORECAST_PANEL_H_
#define _FORECAST_PANEL_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


Rcpp::List forecast_bvarPANEL (
    arma::field<arma::cube>&  posterior_A_c_cpp,      // (S)(K, N, C)
    arma::field<arma::cube>&  posterior_Sigma_c_cpp,  // (S)(N, N, C)
    Rcpp::List&               X_c,                    // (C)(T_c, K)
    const int                 horizon
);


arma::vec mvnrnd_cond (
    arma::vec x,        // Nx1 with NAs or without
    arma::vec mu,       // Nx1 mean vector
    arma::mat Sigma     // NxN covariance matrix
);


Rcpp::List forecast_conditional_bvarPANEL (
    arma::field<arma::cube>&  posterior_A_c_cpp,      // (S)(K, N, C)
    arma::field<arma::cube>&  posterior_Sigma_c_cpp,  // (S)(N, N, C)
    Rcpp::List&               X_c,                    // (C)(T_c, K)
    Rcpp::List&               cond_forecasts,         // (C)(horizon, N)
    const int                 horizon
);


#endif  // _FORECAST_PANEL_H_
