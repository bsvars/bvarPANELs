#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;


// [[Rcpp:interface(cpp,r)]]
// [[Rcpp::export]]
arma::field<arma::mat> rmniw1(
    const arma::mat& A,     // KxN
    const arma::mat& V,     // KxK
    const arma::mat& S,     // NxN
    const double&    nu     // scalar
) {
  mat Sigma     = iwishrnd(S, nu);
  mat X_tmp     = mat(size(A), fill::randn);
  mat X         = A + chol(S).t() * X_tmp * chol(V);
  
  field<mat> out(2);
  out(0)        = X;
  out(1)        = Sigma;
  return out;
} // END rmniw1
