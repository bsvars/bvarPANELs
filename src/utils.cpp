
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::interfaces(cpp)]]
// [[Rcpp::export]]
std::string ordinal(
    int n
) {
  std::string suffix;
  if (n % 10 == 1 && n % 100 != 11) {
    suffix = "st";
  } else if (n % 10 == 2 && n % 100 != 12) {
    suffix = "nd";
  } else if (n % 10 == 3 && n % 100 != 13) {
    suffix = "rd";
  } else {
    suffix = "th";
  }
  return std::to_string(n) + suffix;
} // END ordinal
