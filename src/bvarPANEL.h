#ifndef _BVARPANEL_H_
#define _BVARPANEL_H_

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


Rcpp::List bvarPANEL(
    const int&                    S,          // No. of posterior draws
    const Rcpp::List&             prior,      // a list of priors
    const Rcpp::List&             starting_values, 
    const int                     thin = 100, // introduce thinning
    const bool                    show_progress = true
);


#endif  // _BVARPANEL_H_
