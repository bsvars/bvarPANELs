// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_bvarPANELs_RCPPEXPORTS_H_GEN_
#define RCPP_bvarPANELs_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace bvarPANELs {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("bvarPANELs", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("bvarPANELs", "_bvarPANELs_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in bvarPANELs");
            }
        }
    }

    inline Rcpp::List forecast_bvarPANEL(arma::field<arma::cube>& posterior_A_c_cpp, arma::field<arma::cube>& posterior_Sigma_c_cpp, Rcpp::List& X_c, const int horizon) {
        typedef SEXP(*Ptr_forecast_bvarPANEL)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_forecast_bvarPANEL p_forecast_bvarPANEL = NULL;
        if (p_forecast_bvarPANEL == NULL) {
            validateSignature("Rcpp::List(*forecast_bvarPANEL)(arma::field<arma::cube>&,arma::field<arma::cube>&,Rcpp::List&,const int)");
            p_forecast_bvarPANEL = (Ptr_forecast_bvarPANEL)R_GetCCallable("bvarPANELs", "_bvarPANELs_forecast_bvarPANEL");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_forecast_bvarPANEL(Shield<SEXP>(Rcpp::wrap(posterior_A_c_cpp)), Shield<SEXP>(Rcpp::wrap(posterior_Sigma_c_cpp)), Shield<SEXP>(Rcpp::wrap(X_c)), Shield<SEXP>(Rcpp::wrap(horizon)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline arma::vec mvnrnd_cond(arma::vec x, arma::vec mu, arma::mat Sigma) {
        typedef SEXP(*Ptr_mvnrnd_cond)(SEXP,SEXP,SEXP);
        static Ptr_mvnrnd_cond p_mvnrnd_cond = NULL;
        if (p_mvnrnd_cond == NULL) {
            validateSignature("arma::vec(*mvnrnd_cond)(arma::vec,arma::vec,arma::mat)");
            p_mvnrnd_cond = (Ptr_mvnrnd_cond)R_GetCCallable("bvarPANELs", "_bvarPANELs_mvnrnd_cond");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mvnrnd_cond(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(mu)), Shield<SEXP>(Rcpp::wrap(Sigma)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline Rcpp::List forecast_conditional_bvarPANEL(arma::field<arma::cube>& posterior_A_c_cpp, arma::field<arma::cube>& posterior_Sigma_c_cpp, Rcpp::List& X_c, Rcpp::List& cond_forecasts, const int horizon) {
        typedef SEXP(*Ptr_forecast_conditional_bvarPANEL)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_forecast_conditional_bvarPANEL p_forecast_conditional_bvarPANEL = NULL;
        if (p_forecast_conditional_bvarPANEL == NULL) {
            validateSignature("Rcpp::List(*forecast_conditional_bvarPANEL)(arma::field<arma::cube>&,arma::field<arma::cube>&,Rcpp::List&,Rcpp::List&,const int)");
            p_forecast_conditional_bvarPANEL = (Ptr_forecast_conditional_bvarPANEL)R_GetCCallable("bvarPANELs", "_bvarPANELs_forecast_conditional_bvarPANEL");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_forecast_conditional_bvarPANEL(Shield<SEXP>(Rcpp::wrap(posterior_A_c_cpp)), Shield<SEXP>(Rcpp::wrap(posterior_Sigma_c_cpp)), Shield<SEXP>(Rcpp::wrap(X_c)), Shield<SEXP>(Rcpp::wrap(cond_forecasts)), Shield<SEXP>(Rcpp::wrap(horizon)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

}

#endif // RCPP_bvarPANELs_RCPPEXPORTS_H_GEN_
