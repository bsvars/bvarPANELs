#ifndef _FORECAST_PANEL_H_
#define _FORECAST_PANEL_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


arma::vec mvnrnd_truncated (
    arma::vec     mu,       // Nx1 mean vector
    arma::mat     Sigma,    // NxN covariance matrix
    arma::vec     LB,       // Nx1 lower bounds for truncation
    arma::vec     UB        // Nx1 upper bounds for truncation
);


arma::vec mvnrnd_cond_truncated (
    arma::vec x,        // Nx1 with NAs or without
    arma::vec mu,       // Nx1 mean vector
    arma::mat Sigma,    // NxN covariance matrix
    arma::vec LB,       // Nx1 lower bounds for truncation
    arma::vec UB        // Nx1 upper bounds for truncation
);


Rcpp::List forecast_bvarPANEL (
    arma::field<arma::cube>&  posterior_A_c_cpp,      // (S)(K, N, C)
    arma::field<arma::cube>&  posterior_Sigma_c_cpp,  // (S)(N, N, C)
    Rcpp::List&               X_c,                    // (C)(T_c, K)
    Rcpp::List&               cond_forecasts,         // (C)(horizon, N)
    Rcpp::List&               exog_forecasts,         // (C)(horizon, d)
    const int                 horizon,
    arma::vec                 LB,                     // Nx1 lower bounds for truncation
    arma::vec                 UB,                     // Nx1 upper bounds for truncation
    const bool                show_progress
);


#endif  // _FORECAST_PANEL_H_
