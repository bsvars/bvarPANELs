// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rmniw1
arma::field<arma::mat> rmniw1(const arma::mat& A, const arma::mat& V, const arma::mat& S, const double& nu);
RcppExport SEXP _bvarPANELs_rmniw1(SEXP ASEXP, SEXP VSEXP, SEXP SSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rmniw1(A, V, S, nu));
    return rcpp_result_gen;
END_RCPP
}
// sample_m
double sample_m(const arma::mat& aux_A, const arma::mat& aux_V, const double& aux_s, const double& aux_w, const Rcpp::List& prior);
RcppExport SEXP _bvarPANELs_sample_m(SEXP aux_ASEXP, SEXP aux_VSEXP, SEXP aux_sSEXP, SEXP aux_wSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_A(aux_ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_V(aux_VSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_s(aux_sSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_w(aux_wSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_m(aux_A, aux_V, aux_s, aux_w, prior));
    return rcpp_result_gen;
END_RCPP
}
// sample_w
double sample_w(const arma::mat& aux_V, const Rcpp::List& prior);
RcppExport SEXP _bvarPANELs_sample_w(SEXP aux_VSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_V(aux_VSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_w(aux_V, prior));
    return rcpp_result_gen;
END_RCPP
}
// sample_s
double sample_s(const arma::mat& aux_A, const arma::mat& aux_V, const arma::mat& aux_Sigma, const double& aux_m, const Rcpp::List& prior);
RcppExport SEXP _bvarPANELs_sample_s(SEXP aux_ASEXP, SEXP aux_VSEXP, SEXP aux_SigmaSEXP, SEXP aux_mSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_A(aux_ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_V(aux_VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_Sigma(aux_SigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_m(aux_mSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_s(aux_A, aux_V, aux_Sigma, aux_m, prior));
    return rcpp_result_gen;
END_RCPP
}
// sample_Sigma
arma::mat sample_Sigma(const arma::cube& aux_Sigma_c_inv, const double& aux_s, const double& aux_nu, const Rcpp::List& prior);
RcppExport SEXP _bvarPANELs_sample_Sigma(SEXP aux_Sigma_c_invSEXP, SEXP aux_sSEXP, SEXP aux_nuSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_Sigma_c_inv(aux_Sigma_c_invSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_s(aux_sSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_nu(aux_nuSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_Sigma(aux_Sigma_c_inv, aux_s, aux_nu, prior));
    return rcpp_result_gen;
END_RCPP
}
// sample_AV
arma::field<arma::mat> sample_AV(const arma::cube& aux_A_c, const arma::cube& aux_Sigma_c_inv, const double& aux_s, const double& aux_m, const double& aux_w, const Rcpp::List& prior);
RcppExport SEXP _bvarPANELs_sample_AV(SEXP aux_A_cSEXP, SEXP aux_Sigma_c_invSEXP, SEXP aux_sSEXP, SEXP aux_mSEXP, SEXP aux_wSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_A_c(aux_A_cSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_Sigma_c_inv(aux_Sigma_c_invSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_s(aux_sSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_m(aux_mSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_w(aux_wSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_AV(aux_A_c, aux_Sigma_c_inv, aux_s, aux_m, aux_w, prior));
    return rcpp_result_gen;
END_RCPP
}
// sample_A_c_Sigma_c
arma::field<arma::mat> sample_A_c_Sigma_c(const arma::mat& Y_c, const arma::mat& X_c, const arma::mat& aux_A, const arma::mat& aux_V, const arma::mat& aux_Sigma, const double& aux_nu);
RcppExport SEXP _bvarPANELs_sample_A_c_Sigma_c(SEXP Y_cSEXP, SEXP X_cSEXP, SEXP aux_ASEXP, SEXP aux_VSEXP, SEXP aux_SigmaSEXP, SEXP aux_nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y_c(Y_cSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X_c(X_cSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_A(aux_ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_V(aux_VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_Sigma(aux_SigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_nu(aux_nuSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_A_c_Sigma_c(Y_c, X_c, aux_A, aux_V, aux_Sigma, aux_nu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bvarPANELs_rmniw1", (DL_FUNC) &_bvarPANELs_rmniw1, 4},
    {"_bvarPANELs_sample_m", (DL_FUNC) &_bvarPANELs_sample_m, 5},
    {"_bvarPANELs_sample_w", (DL_FUNC) &_bvarPANELs_sample_w, 2},
    {"_bvarPANELs_sample_s", (DL_FUNC) &_bvarPANELs_sample_s, 5},
    {"_bvarPANELs_sample_Sigma", (DL_FUNC) &_bvarPANELs_sample_Sigma, 4},
    {"_bvarPANELs_sample_AV", (DL_FUNC) &_bvarPANELs_sample_AV, 6},
    {"_bvarPANELs_sample_A_c_Sigma_c", (DL_FUNC) &_bvarPANELs_sample_A_c_Sigma_c, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_bvarPANELs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
