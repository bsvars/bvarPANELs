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


double sample_s (
    const arma::mat&    aux_A,      // KxN
    const arma::mat&    aux_V,      // KxK
    const arma::mat&    aux_Sigma,  // NxN
    const double&       aux_m,      // scalar
    const Rcpp::List&   prior
);


double log_kernel_nu (
    const double&       aux_nu,           // scalar
    const arma::cube&   aux_Sigma_c_cpp,  // NxNxC
    const arma::cube&   aux_Sigma_c_inv,  // NxNxC
    const arma::mat&    aux_Sigma,        // NxN
    const double&       prior_lambda,     // scalar
    const int&          C,                // scalar
    const int&          N,                // scalar
    const int&          K                 // scalar
);


double mcmc_accpetance_rate1 (
    arma::vec& mcmc
);


double sample_nu (
    const double&       aux_nu,           // scalar
    const arma::vec&    posterior_nu,     // sx1
    const arma::cube&   aux_Sigma_c_cpp,  // NxNxC
    const arma::cube&   aux_Sigma_c_inv,  // NxNxC
    const arma::mat&    aux_Sigma,        // NxN
    const Rcpp::List&   prior,
    const int&          iteration,        // MCMC iteration passed
    arma::vec&          scale,            // (Sx1) adaptive scaling
    const arma::vec&    rate_target_start_initial
);


arma::mat sample_Sigma (
    const arma::cube&   aux_Sigma_c_inv,  // NxNxC
    const double&       aux_s,            // scalar
    const double&       aux_nu,           // scalar
    const Rcpp::List&   prior
);


arma::field<arma::mat> sample_AV (
    const arma::cube&   aux_A_c_cpp,      // KxNxC
    const arma::cube&   aux_Sigma_c_inv,  // NxNxC
    const double&       aux_s,            // scalar
    const double&       aux_m,            // scalar
    const double&       aux_w,            // scalar
    const Rcpp::List&   prior
);


arma::field<arma::mat> sample_A_c_Sigma_c (
    const arma::mat&    Y_c,              // T_cxN
    const arma::mat&    X_c,              // T_cxK
    const arma::mat&    aux_A,            // KxN
    const arma::mat&    aux_V,            // KxK
    const arma::mat&    aux_Sigma,        // NxN
    const double&       aux_nu            // scalar
);


#endif  // _SAMPLE_MNIW_H_
