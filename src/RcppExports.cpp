// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/bvarPANELs.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bvarPANEL
Rcpp::List bvarPANEL(const int& S, const Rcpp::List& Y, const Rcpp::List& X, const Rcpp::List& prior, const Rcpp::List& starting_values, const int thin, const bool show_progress, const arma::vec& rate_target_start_initial);
RcppExport SEXP _bvarPANELs_bvarPANEL(SEXP SSEXP, SEXP YSEXP, SEXP XSEXP, SEXP priorSEXP, SEXP starting_valuesSEXP, SEXP thinSEXP, SEXP show_progressSEXP, SEXP rate_target_start_initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type starting_values(starting_valuesSEXP);
    Rcpp::traits::input_parameter< const int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const bool >::type show_progress(show_progressSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rate_target_start_initial(rate_target_start_initialSEXP);
    rcpp_result_gen = Rcpp::wrap(bvarPANEL(S, Y, X, prior, starting_values, thin, show_progress, rate_target_start_initial));
    return rcpp_result_gen;
END_RCPP
}
// forecast_bvarPANEL
Rcpp::List forecast_bvarPANEL(arma::field<arma::cube>& posterior_A_c_cpp, arma::field<arma::cube>& posterior_Sigma_c_cpp, Rcpp::List& X_c, const int horizon);
static SEXP _bvarPANELs_forecast_bvarPANEL_try(SEXP posterior_A_c_cppSEXP, SEXP posterior_Sigma_c_cppSEXP, SEXP X_cSEXP, SEXP horizonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::field<arma::cube>& >::type posterior_A_c_cpp(posterior_A_c_cppSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::cube>& >::type posterior_Sigma_c_cpp(posterior_Sigma_c_cppSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type X_c(X_cSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    rcpp_result_gen = Rcpp::wrap(forecast_bvarPANEL(posterior_A_c_cpp, posterior_Sigma_c_cpp, X_c, horizon));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bvarPANELs_forecast_bvarPANEL(SEXP posterior_A_c_cppSEXP, SEXP posterior_Sigma_c_cppSEXP, SEXP X_cSEXP, SEXP horizonSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bvarPANELs_forecast_bvarPANEL_try(posterior_A_c_cppSEXP, posterior_Sigma_c_cppSEXP, X_cSEXP, horizonSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// mvnrnd_cond
arma::vec mvnrnd_cond(arma::vec x, arma::vec mu, arma::mat Sigma);
static SEXP _bvarPANELs_mvnrnd_cond_try(SEXP xSEXP, SEXP muSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvnrnd_cond(x, mu, Sigma));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bvarPANELs_mvnrnd_cond(SEXP xSEXP, SEXP muSEXP, SEXP SigmaSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bvarPANELs_mvnrnd_cond_try(xSEXP, muSEXP, SigmaSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// forecast_conditional_bvarPANEL
Rcpp::List forecast_conditional_bvarPANEL(arma::field<arma::cube>& posterior_A_c_cpp, arma::field<arma::cube>& posterior_Sigma_c_cpp, Rcpp::List& X_c, Rcpp::List& cond_forecasts, const int horizon);
static SEXP _bvarPANELs_forecast_conditional_bvarPANEL_try(SEXP posterior_A_c_cppSEXP, SEXP posterior_Sigma_c_cppSEXP, SEXP X_cSEXP, SEXP cond_forecastsSEXP, SEXP horizonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::field<arma::cube>& >::type posterior_A_c_cpp(posterior_A_c_cppSEXP);
    Rcpp::traits::input_parameter< arma::field<arma::cube>& >::type posterior_Sigma_c_cpp(posterior_Sigma_c_cppSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type X_c(X_cSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type cond_forecasts(cond_forecastsSEXP);
    Rcpp::traits::input_parameter< const int >::type horizon(horizonSEXP);
    rcpp_result_gen = Rcpp::wrap(forecast_conditional_bvarPANEL(posterior_A_c_cpp, posterior_Sigma_c_cpp, X_c, cond_forecasts, horizon));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _bvarPANELs_forecast_conditional_bvarPANEL(SEXP posterior_A_c_cppSEXP, SEXP posterior_Sigma_c_cppSEXP, SEXP X_cSEXP, SEXP cond_forecastsSEXP, SEXP horizonSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_bvarPANELs_forecast_conditional_bvarPANEL_try(posterior_A_c_cppSEXP, posterior_Sigma_c_cppSEXP, X_cSEXP, cond_forecastsSEXP, horizonSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error("%s", CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
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
// log_kernel_nu
double log_kernel_nu(const double& aux_nu, const arma::cube& aux_Sigma_c_cpp, const arma::cube& aux_Sigma_c_inv, const arma::mat& aux_Sigma, const double& prior_lambda, const int& C, const int& N, const int& K);
RcppExport SEXP _bvarPANELs_log_kernel_nu(SEXP aux_nuSEXP, SEXP aux_Sigma_c_cppSEXP, SEXP aux_Sigma_c_invSEXP, SEXP aux_SigmaSEXP, SEXP prior_lambdaSEXP, SEXP CSEXP, SEXP NSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type aux_nu(aux_nuSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_Sigma_c_cpp(aux_Sigma_c_cppSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_Sigma_c_inv(aux_Sigma_c_invSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_Sigma(aux_SigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type prior_lambda(prior_lambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(log_kernel_nu(aux_nu, aux_Sigma_c_cpp, aux_Sigma_c_inv, aux_Sigma, prior_lambda, C, N, K));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_accpetance_rate1
double mcmc_accpetance_rate1(arma::vec& mcmc);
RcppExport SEXP _bvarPANELs_mcmc_accpetance_rate1(SEXP mcmcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type mcmc(mcmcSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_accpetance_rate1(mcmc));
    return rcpp_result_gen;
END_RCPP
}
// sample_nu
double sample_nu(const double& aux_nu, const arma::vec& posterior_nu, const arma::cube& aux_Sigma_c_cpp, const arma::cube& aux_Sigma_c_inv, const arma::mat& aux_Sigma, const Rcpp::List& prior, const int& iteration, arma::vec& scale, const arma::vec& rate_target_start_initial);
RcppExport SEXP _bvarPANELs_sample_nu(SEXP aux_nuSEXP, SEXP posterior_nuSEXP, SEXP aux_Sigma_c_cppSEXP, SEXP aux_Sigma_c_invSEXP, SEXP aux_SigmaSEXP, SEXP priorSEXP, SEXP iterationSEXP, SEXP scaleSEXP, SEXP rate_target_start_initialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type aux_nu(aux_nuSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type posterior_nu(posterior_nuSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_Sigma_c_cpp(aux_Sigma_c_cppSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_Sigma_c_inv(aux_Sigma_c_invSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type aux_Sigma(aux_SigmaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const int& >::type iteration(iterationSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rate_target_start_initial(rate_target_start_initialSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_nu(aux_nu, posterior_nu, aux_Sigma_c_cpp, aux_Sigma_c_inv, aux_Sigma, prior, iteration, scale, rate_target_start_initial));
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
arma::field<arma::mat> sample_AV(const arma::cube& aux_A_c_cpp, const arma::cube& aux_Sigma_c_inv, const double& aux_s, const double& aux_m, const double& aux_w, const Rcpp::List& prior);
RcppExport SEXP _bvarPANELs_sample_AV(SEXP aux_A_c_cppSEXP, SEXP aux_Sigma_c_invSEXP, SEXP aux_sSEXP, SEXP aux_mSEXP, SEXP aux_wSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_A_c_cpp(aux_A_c_cppSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type aux_Sigma_c_inv(aux_Sigma_c_invSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_s(aux_sSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_m(aux_mSEXP);
    Rcpp::traits::input_parameter< const double& >::type aux_w(aux_wSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_AV(aux_A_c_cpp, aux_Sigma_c_inv, aux_s, aux_m, aux_w, prior));
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

// validate (ensure exported C++ functions exist before calling them)
static int _bvarPANELs_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("Rcpp::List(*forecast_bvarPANEL)(arma::field<arma::cube>&,arma::field<arma::cube>&,Rcpp::List&,const int)");
        signatures.insert("arma::vec(*mvnrnd_cond)(arma::vec,arma::vec,arma::mat)");
        signatures.insert("Rcpp::List(*forecast_conditional_bvarPANEL)(arma::field<arma::cube>&,arma::field<arma::cube>&,Rcpp::List&,Rcpp::List&,const int)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _bvarPANELs_RcppExport_registerCCallable() { 
    R_RegisterCCallable("bvarPANELs", "_bvarPANELs_forecast_bvarPANEL", (DL_FUNC)_bvarPANELs_forecast_bvarPANEL_try);
    R_RegisterCCallable("bvarPANELs", "_bvarPANELs_mvnrnd_cond", (DL_FUNC)_bvarPANELs_mvnrnd_cond_try);
    R_RegisterCCallable("bvarPANELs", "_bvarPANELs_forecast_conditional_bvarPANEL", (DL_FUNC)_bvarPANELs_forecast_conditional_bvarPANEL_try);
    R_RegisterCCallable("bvarPANELs", "_bvarPANELs_RcppExport_validate", (DL_FUNC)_bvarPANELs_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_bvarPANELs_bvarPANEL", (DL_FUNC) &_bvarPANELs_bvarPANEL, 8},
    {"_bvarPANELs_forecast_bvarPANEL", (DL_FUNC) &_bvarPANELs_forecast_bvarPANEL, 4},
    {"_bvarPANELs_mvnrnd_cond", (DL_FUNC) &_bvarPANELs_mvnrnd_cond, 3},
    {"_bvarPANELs_forecast_conditional_bvarPANEL", (DL_FUNC) &_bvarPANELs_forecast_conditional_bvarPANEL, 5},
    {"_bvarPANELs_rmniw1", (DL_FUNC) &_bvarPANELs_rmniw1, 4},
    {"_bvarPANELs_sample_m", (DL_FUNC) &_bvarPANELs_sample_m, 5},
    {"_bvarPANELs_sample_w", (DL_FUNC) &_bvarPANELs_sample_w, 2},
    {"_bvarPANELs_sample_s", (DL_FUNC) &_bvarPANELs_sample_s, 5},
    {"_bvarPANELs_log_kernel_nu", (DL_FUNC) &_bvarPANELs_log_kernel_nu, 8},
    {"_bvarPANELs_mcmc_accpetance_rate1", (DL_FUNC) &_bvarPANELs_mcmc_accpetance_rate1, 1},
    {"_bvarPANELs_sample_nu", (DL_FUNC) &_bvarPANELs_sample_nu, 9},
    {"_bvarPANELs_sample_Sigma", (DL_FUNC) &_bvarPANELs_sample_Sigma, 4},
    {"_bvarPANELs_sample_AV", (DL_FUNC) &_bvarPANELs_sample_AV, 6},
    {"_bvarPANELs_sample_A_c_Sigma_c", (DL_FUNC) &_bvarPANELs_sample_A_c_Sigma_c, 6},
    {"_bvarPANELs_RcppExport_registerCCallable", (DL_FUNC) &_bvarPANELs_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_bvarPANELs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
