#  #####################################################################################
#  R package bvarPANELs by Tomasz Woźniak and Miguel Sanchez-Martinez
#  Copyright (C) 2024 International Labour Organization
#
#  This file is part of the R package bvarPANELs: Forecasting with Bayesian 
#  Hierarchical Panel Vector Autoregressions
#
#  The R package bvarPANELs is free software: you can redistribute it
#  and/or modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation, either version 3 or
#  any later version of the License.
#
#  The R package bvarPANELs is distributed in the hope that it will be
#  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with the R package bvarPANELs If that is not the case, please
#  refer to <http://www.gnu.org/licenses/>.
#  #####################################################################################
#
#' @title Forecasting with Bayesian Hierarchical Panel Vector Autoregressions
#'
#' @description Provides Bayesian estimation and forecasting of dynamic panel 
#' data using Bayesian Hierarchical Panel Vector Autoregressions. The model 
#' includes country-specific VARs that share a global prior distribution. Under 
#' this prior expected value, each country's system follows a global VAR with 
#' country-invariant parameters. Further flexibility is provided by the 
#' hierarchical prior structure that retains the Minnesota prior interpretation 
#' for the global VAR and features estimated prior covariance matrices, 
#' shrinkage, and persistence levels. Bayesian forecasting is developed for 
#' models including exogenous variables, allowing conditional forecasts given 
#' the future trajectories of some variables and restricted forecasts assuring 
#' that rates are forecasted to stay positive and less than 100. The package 
#' implements the model specification, estimation, and forecasting routines, 
#' facilitating coherent workflows and reproducibility. Beautiful plots, 
#' informative summary functions, and extensive documentation complement all 
#' this. An extraordinary computational speed is achieved thanks to employing 
#' frontier econometric and numerical techniques and algorithms written in C++. 
#' The 'bvarPANELs' package is aligned regarding objects, workflows, and code 
#' structure with the R packages 'bsvars' by Woźniak (2024) 
#' \doi{10.32614/CRAN.package.bsvars} and 'bsvarSIGNs' by Wang & Woźniak (2024) 
#' \doi{10.32614/CRAN.package.bsvarSIGNs}, and they constitute an integrated 
#' toolset. Copyright: 2024 International Labour Organization.
#' 
#' @details 
#' The package provides a set of functions for predictive analysis with the 
#' Bayesian Hierarchical Panel Vector Autoregression. 
#' 
#' \strong{The Model.} The model specification is initiated using function 
#' \code{\link{specify_bvarPANEL}} that creates an object of class 
#' \code{BVARPANEL} containing the prior specification, starting values for 
#' estimation, data matrices, and the setup of the Monte Carlo Markov 
#' Chain sampling algorithm.
#' 
#' The model features country-specific Vector Autoregressive 
#' (VAR) equation for \code{N} dependent variables with \code{T_c} observations 
#' for each country \code{c}. Its equation is given by
#' \deqn{\mathbf{Y}_c = \mathbf{A}_c\mathbf{X}_c + \mathbf{E}_c}
#' where \eqn{\mathbf{Y}_c} is an \code{T_c x N} matrix of dependent variables for 
#' country \code{c}, \eqn{\mathbf{X}_c} is a \code{T_c x K} matrix of explanatory 
#' variables, \eqn{\mathbf{E}_c} is an \code{T_c x N} matrix of error terms, and 
#' \eqn{\mathbf{A}_c} is an \code{NxK} matrix of country-specific autoregressive slope 
#' coefficients and parameters on deterministic terms in \eqn{\mathbf{X}_c}. The 
#' parameter matrix \eqn{\mathbf{A}_c} includes autoregressive matrices capturing the 
#' effects of the lagged vectors of dependent variables at lags from 1 to \code{p},
#' a constant term and a set of exogenous variables.
#' 
#' The error terms for each of the periods have zero conditional mean and
#' conditional covariance given by the \code{NxN} matrix \eqn{\mathbf{\Sigma}}. The errors 
#' are jointly normally distributed and serially uncorrelated. These 
#' assumptions are summarised using a matrix-variate normal distribution 
#' (see Woźniak, 2016):
#' \deqn{\mathbf{E}_c \sim MN(\mathbf{0}_{T_c x N}, \mathbf{\Sigma}, \mathbf{I}_{T_c})}
#' where the identity matrix \eqn{\mathbf{I}_{T_c}} of order \code{T_c} and joint 
#' normality imply no serial autocorrelation. Matrix \eqn{\mathbf{0}_{T_c x N}}
#' denotes a \code{T_c x N} matrix of zeros.
#' 
#' \strong{Global Prior Distributions.} 
#' The Hierarchical Panel VAR model features a sophisticated hierarchical prior
#' structure that grants the model flexibility, interpretability, and improved 
#' forecasting performance.
#' 
#' The country-specific parameters follow a prior distribution that, at its 
#' mean value, represents a global VAR model with a global autoregressive 
#' parameter matrix \eqn{\mathbf{A}} of dimension \code{KxN} and an \code{NxN} global 
#' covariance matrix \eqn{\mathbf{\Sigma}}:
#' \deqn{\mathbf{Y}_c = \mathbf{AX}_c + \mathbf{E}_c}
#' This global VAR model under the prior mean is represented by the parameters 
#' of the matrix-variate normal inverted Wishart distribution (see Woźniak, 2016)
#' given by:
#' \deqn{\mathbf{A}_c, \mathbf{\Sigma}_c | \mathbf{A}, \mathbf{V}, \mathbf{\Sigma}, \nu \sim MNIW(\mathbf{A}, \mathbf{V}, (N - \nu - 1)\mathbf{\Sigma}, \nu)}
#' where \eqn{V} is a \code{KxK} column-specific covariance matrix, 
#' \eqn{(N - \nu - 1)\mathbf{\Sigma}} is the row-specific matrix, and \eqn{\nu > N+1} is the 
#' degrees-of-freedom parameter. 
#' 
#' All of the parameters of the prior distribution above feature their own prior
#' distributions and are estimated. These prior distributions are given by:
#' \deqn{\mathbf{A} | \mathbf{V}, m, w, s \sim MN (m\underline{\mathbf{M}}, \mathbf{V}, s\underline{\mathbf{S}} ) }
#' with the \code{KxN} mean matrix \eqn{m\underline{\mathbf{M}}}, the \code{KxK}
#' column-specific covariance matrix \eqn{\mathbf{V}}, and the \code{NxN} matrix of
#' row-specific covariance \eqn{s\underline{\mathbf{S}}}.
#' 
#' The global error term covariance matrix, \eqn{\mathbf{\Sigma}}, follows a Wishart 
#' distribution with \code{NxN} scale matrix \eqn{s\underline{\mathbf{S}}_\Sigma} 
#' and shape parameter \eqn{\underline{\mu}_\Sigma}
#' \deqn{\mathbf{\Sigma} | s, \nu \sim W(s\underline{\mathbf{S}}_\Sigma,\underline{\mu}_\Sigma)}
#' 
#' \strong{Other Prior Distributions.} 
#' The column-specific covariance \eqn{\mathbf{V}} follows the inverse-Wishart 
#' distribution with scale \eqn{w\underline{\mathbf{W}}} and shape 
#' \eqn{\underline{\eta}}:
#' \deqn{\mathbf{V} | w \sim IW(w\underline{\mathbf{W}}, \underline{\eta})}
#' The shape parameter \eqn{\nu} follows an exponential distribution with mean 
#' \eqn{\underline\lambda}:
#' \deqn{\nu \sim\exp(\underline\lambda)}
#' Finally, the priors for the remaining scalar hyper-parameters are:
#' \deqn{m \sim N(\underline{\mu}_m, \underline{\sigma}_m^2)}
#' \deqn{w \sim G(\underline{s}_w, \underline{a}_w)}
#' \deqn{s \sim IG2 (\underline{s}_s, \underline{\nu}_s)}
#' 
#' The prior hyper-parameters in this note are grouped into those that are:
#' \describe{
#'  \item{fixed}{and denoted using underscore, such as e.g. 
#'  \eqn{\underline{\mathbf{M}}}, \eqn{\underline{\mu}_\Sigma}, or 
#'  \eqn{\underline{\nu}_s}. These hyper-parameters must be fixed and their 
#'  default values are set by initiating the model specification using function
#'  \code{\link{specify_bvarPANEL}}. These values can be accesses from such 
#'  generated object in its element \code{prior} and can be modified by the user.}
#'  \item{estimated}{not featuring the underscore in the notation, such as e.g.
#'  \eqn{\mathbf{A}}, \eqn{\mathbf{\Sigma}}, or \eqn{m}. These hyper-parameters
#'  are estimated and their posterior draws are available from an object 
#'  generated after the estimation running the function \code{estimate}.}
#' }
#' 
#' \strong{Estimation.}
#' The package implements Bayesian estimation using the Gibbs sampler. This 
#' algorithm provides a sample of random draws from the posterior distribution 
#' of the parameters of the model. The posterior distribution is defined by 
#' Bayes' Rule stating that the posterior
#' distribution of the parameters given data and is proportional to the likelihood 
#' function and the prior distribution of the parameters:
#' \deqn{p(\mathbf{\theta} | \mathbf{Y}) \propto L(\mathbf{\theta}; \mathbf{Y}) p(\mathbf{\theta})}
#' where \eqn{\mathbf{\theta}} collects all the parameters of the model to be 
#' estimated. At each of its iterations a single draw of all of the 
#' parameters of the model, including the estimated hyper-parameters, is obtained. 
#' This Bayesian procedure estimates jointly all the parameters of the model and
#' is implemented in the \code{\link{estimate.BVARPANEL}} and 
#' \code{\link{estimate.PosteriorBVARPANEL}} functions.
#' 
#' \strong{Forecasting.} 
#' The package implements Bayesian forecasting providing the a sample of draws
#' from the joint predictive density defined as the joint density of the future
#' unknown values to be predicted, \eqn{\mathbf{Y}_f}, given data, \eqn{\mathbf{Y}}
#' closely following Karlsson (2013):
#' \deqn{p(\mathbf{Y}_f | \mathbf{Y})}
#' The package offers the possibility of:
#' \describe{
#'   \item{forecasting for models with exogenous variables}{given the provided 
#'   future values of the exogenous variables.}
#'   \item{conditional predictions}{given provided future projections for some 
#'   of the variables.}
#'   \item{trucated forecasts}{for variables that represents rates from the 
#'   interval \eqn{[0,100]}.}
#' }
#' The forecasting is performed using function \code{\link{forecast.PosteriorBVARPANEL}}.
#' 
#' @seealso \code{\link{specify_bvarPANEL}}, \code{\link{estimate.BVARPANEL}}, 
#' \code{\link{forecast.PosteriorBVARPANEL}}
#' 
#' @references
#' Karlsson, S. (2013). Forecasting with Bayesian Vector Autoregression, 
#' in: \emph{Handbook of Economic Forecasting}, Elsevier. volume \bold{2}, 791–897,
#' \doi{10.1016/B978-0-444-62731-5.00015-4}.
#' 
#' Woźniak, T. (2016). Bayesian Vector Autoregressions, 
#' \emph{Australian Economic Review}, \bold{49}, 365-380, 
#' \doi{10.1111/1467-8462.12179}.
#' 
#' @name bvarPANELs-package
#' @aliases bvarPANELs-package bvarPANELs
#' @docType package
#' @useDynLib bvarPANELs, .registration = TRUE
#' @importFrom bsvars estimate forecast compute_variance_decompositions
#' @importFrom Rcpp sourceCpp
#' @importFrom R6 R6Class
#' @importFrom RcppTN rtn dtn
#' @importFrom stats sd quantile
#' @importFrom tmvtnsim rtmvnorm
#' @import RcppProgress
#' @note This package is currently in active development.
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' @keywords package models ts
#' @examples
#' # Basic estimation and forecasting example
#' ############################################################
#' data(ilo_dynamic_panel)                                   # load the data
#' set.seed(123)
#' specification = specify_bvarPANEL$new(ilo_dynamic_panel) # specify the model
#' burn_in       = estimate(specification, S = 10)         # run the burn-in; use say S = 10000
#' posterior     = estimate(burn_in, S = 10)               # estimate the model; use say S = 10000
#' predictive    = forecast(posterior, 2)                  # forecast the future       
#' 
#' # workflow with the pipe |>
#' set.seed(123)
#' ilo_dynamic_panel |>
#'   specify_bvarPANEL$new() |>
#'   estimate(S = 20) |> 
#'   estimate(S = 20) |> 
#'   forecast(horizon = 2) -> predictive
#'   
#' plot(predictive, which_c = "POL")
#' 
#' # Full estimation and forecasting example with 
#' #   exogenous variables, conditional forecasts, and truncation for rates
#' ############################################################
#' data(ilo_dynamic_panel)                                 # load the data
#' data(ilo_exogenous_variables)                           # load the exogenous variables
#' data(ilo_exogenous_forecasts)                           # load the exogenous forecasts
#' data(ilo_conditional_forecasts)                         # load the conditional forecasts
#' set.seed(123)
#' specification = specify_bvarPANEL$new(
#'                   ilo_dynamic_panel,
#'                   exogenous = ilo_exogenous_variables,
#'                   type = c("real", rep("rates", 3))
#'                 )
#' burn_in       = estimate(specification, S = 10)         # run the burn-in; use say S = 10000
#' posterior     = estimate(burn_in, S = 10)               # estimate the model; use say S = 10000
#' predictive    = forecast(
#'                   posterior, 
#'                   horizon = 6,
#'                   exogenous_forecast = ilo_exogenous_forecasts,
#'                   conditional_forecast = ilo_conditional_forecasts
#'                 )
#'                 
#' plot(predictive, which_c = "POL")
#' 
"_PACKAGE"
