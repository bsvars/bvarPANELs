// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rmn_cpp
arma::mat rmn_cpp(const arma::mat& A, const arma::mat& S, const arma::mat& V);
RcppExport SEXP _bvarPANELs_rmn_cpp(SEXP ASEXP, SEXP SSEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(rmn_cpp(A, S, V));
    return rcpp_result_gen;
END_RCPP
}
// riw1_cpp
arma::mat riw1_cpp(const arma::mat& S, const double& nu);
RcppExport SEXP _bvarPANELs_riw1_cpp(SEXP SSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(riw1_cpp(S, nu));
    return rcpp_result_gen;
END_RCPP
}
// riw2_cpp
arma::mat riw2_cpp(const arma::mat& S, const double& nu);
RcppExport SEXP _bvarPANELs_riw2_cpp(SEXP SSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(riw2_cpp(S, nu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bvarPANELs_rmn_cpp", (DL_FUNC) &_bvarPANELs_rmn_cpp, 3},
    {"_bvarPANELs_riw1_cpp", (DL_FUNC) &_bvarPANELs_riw1_cpp, 2},
    {"_bvarPANELs_riw2_cpp", (DL_FUNC) &_bvarPANELs_riw2_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_bvarPANELs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
