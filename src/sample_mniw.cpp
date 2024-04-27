#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;


// [[Rcpp:interface(cpp,r)]]
// [[Rcpp::export]]
arma::mat rmn_cpp(
    const arma::mat& A,     // KxN
    const arma::mat& S,     // NxN
    const arma::mat& V      // KxK
) {
  
  mat X = mat(size(A), fill::randn);
  return A + chol(S).t() * X * chol(V);
} // END rmn_cpp


// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat riw1_cpp (
    const arma::mat& S,     // NxN
    const double&    nu     // scalar
) {
  // Based on Macroeconometrics notes Lecture 10
  
  int N           = S.n_cols;
  
  mat s_chol      = chol(S, "lower");
  
  mat Q(N, N, fill::zeros);
  Q.diag()        = sqrt(pow(chi2rnd(nu - N + 1, N), -1));
  
  for (int i = 0; i < (N - 1); i++) {
    Q.submat(i + 1, i, N - 1, i) = randn(N - i - 1);
  }
  mat Q_inv       = inv(trimatu(Q));
  
  return s_chol * Q_inv.t() * Q_inv * s_chol.t();
} // END riw1_cpp




// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat riw2_cpp (
    const arma::mat& S,     // NxN
    const double&    nu     // scalar
) {
  // Based on algorithm B.4.4. from Appendinx B by Bauwens, Lubrano, Richard (1999) Bayesian Inference in Dynamic Econometric Models, Oxford Uni Press
  
  int N           = S.n_cols;
  
  mat s_chol      = chol(S, "lower");
  
  mat Q(N, N, fill::zeros);
  Q.diag()        = sqrt(pow(chi2rnd(nu - N + 1, N), -1));
  
  for (int i = 0; i < (N - 1); i++) {
    Q.submat(i + 1, i, N - 1, i) = randn(N - i - 1);
  }
  mat Q_inv       = inv(trimatu(Q));
  
  return s_chol * Q_inv.t() * Q_inv * s_chol.t();
} // END riw2_cpp
