#ifndef _SAMPLE_MNIW_H_
#define _SAMPLE_MNIW_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


arma::field<arma::mat> rmniw1(
    const arma::mat& A,     // KxN
    const arma::mat& V,     // KxK
    const arma::mat& S,     // NxN
    const double&    nu     // scalar
);


double sample_m (
    const arma::mat&    aux_A,    // KxN
    const arma::mat&    aux_V,    // KxK
    const double&       aux_s,   // scalar
    const double&       aux_w,   // scalar
    const Rcpp::List&   prior
);


double sample_w (
    const arma::mat&    aux_V,    // KxK
    const Rcpp::List&   prior
);

#endif  // _SAMPLE_MNIW_H_
