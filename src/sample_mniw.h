#ifndef _SAMPLE_MNIW_H_
#define _SAMPLE_MNIW_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

arma::mat rmn_cpp(
    const arma::mat& A,     // KxN
    const arma::mat& S,     // NxN
    const arma::mat& V      // KxK
);


arma::mat riw1_cpp (
    const arma::mat& S,     // NxN
    const double&    nu     // scalar
);


arma::mat riw2_cpp (
    const arma::mat& S,     // NxN
    const double&    nu     // scalar
);

#endif  // _SAMPLE_MNIW_H_
