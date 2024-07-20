#ifndef _PANEL_VARIANCE_DECOMPOSITIONS_H_
#define _PANEL_VARIANCE_DECOMPOSITIONS_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


arma::cube Sigma2B_c (
    arma::cube&   posterior_Sigma_c,    // (N, N, S)
    const bool    lower = true
);


arma::cube flip_cube_rows_cols (
    arma::cube&   x   // (N, K, S)
);


arma::field<arma::cube> panel_variance_decompositions (
    arma::field<arma::cube>&  posterior_Sigma,    // (C)(N, N, S)
    arma::field<arma::cube>&  posterior_A,        // (C)(K, N, S)
    arma::cube&               global_Sigma,       // (N, N, S)
    arma::cube&               global_A,           // (K, N, S)
    const int     horizon,
    const int     p,
    const bool    lower = true
);


#endif  // _PANEL_VARIANCE_DECOMPOSITIONS_H_
