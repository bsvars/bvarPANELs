#' @title Forecasting using Hierarchical Panel Vector Autoregressions
#'
#' @description Samples from the joint predictive density of the dependent 
#' variables for all countries at forecast horizons 
#' from 1 to \code{horizon} specified as an argument of the function. 
#' Also implements conditional forecasting based on the provided projections
#' for some of the variables.
#' 
#' @details 
#' The package provides a range of options regarding the forecasting procedure.
#' They are dependent on the model and forecast specifications and include 
#' Bayesian forecasting many periods ahead, conditional forecasting, and 
#' forecasting for models with exogenous variables.
#' 
#' \strong{One-period-ahead predictive density.}
#' The model assumptions provided in the documentation for \code{\link{bvarPANELs}} 
#' determine the country-specific one-period ahead conditional predictive density 
#' for the unknown vector \eqn{\mathbf{y}_{c.t+1}} given the data available at 
#' time \eqn{t} and the parameters of the model. It is multivariate normal with
#' the mean \eqn{\mathbf{A}_c' \mathbf{x}_{c.t+1}} and the covariance matrix 
#' \eqn{\mathbf{\Sigma}_c}
#' \deqn{p(\mathbf{y}_{c.t+1} | \mathbf{x}_{c.t+1}, \mathbf{A}_c, \mathbf{\Sigma}_c) = N_N(\mathbf{A}_c' \mathbf{x}_{c.t+1}, \mathbf{\Sigma}_c)}
#' where \eqn{\mathbf{x}_{c.t+1}} includes the lagged
#' values of \eqn{\mathbf{y}_{c.t+1}}, the constant term, and, potentially,
#' exogenous variables if they were specified by the user. 
#' 
#' \strong{Bayesian predictive density.}
#' The one-period ahead predictive density is used to sample from the joint 
#' predictive density of the unknown future values. This predictive density is
#' defined as a joint density of \eqn{\mathbf{y}_{c.t+h}} at horizons 
#' \eqn{h = 1,\dots,H}, where \eqn{H} corresponds to the value of argument 
#' \code{horizon}, given the data available at time \eqn{t}:
#' \deqn{p( \mathbf{y}_{c.T_c + H}, \dots, \mathbf{y}_{c.T_c + 1} | \mathbf{Y}_c, \mathbf{X}_c) = 
#' \int p(\mathbf{y}_{c.T_c + H}, \dots, \mathbf{y}_{c.T_c + 1} | \mathbf{Y}_c, \mathbf{X}_c, \mathbf{A}_c, \boldsymbol\Sigma_c)
#' p( \mathbf{A}_c, \boldsymbol\Sigma_c | \mathbf{Y}_c, \mathbf{X}_c) d(\mathbf{A}_c, \boldsymbol\Sigma_c)}
#' Therefore, the Bayesian forecast does not depend on the parameter values as
#' the parameters are integrated out with respect to their posterior distribution.
#' Consequently, Bayesian forecasts incorporate the uncertainty with respect to
#' estimation. Sampling from the density is facilitated using the draws from the 
#' posterior density and sequential sampling from the one-period ahead 
#' predictive density.
#' 
#' \strong{Conditional forecasting} of some of the variables given the future 
#' values of the remaining variables is implemented following 
#' Waggoner and Zha (1999) and is based on the conditional normal density given
#' the future projections of some of the variables created basing on the 
#' one-period ahead predictive density.
#' 
#' \strong{Exogenous variables.}
#' Forecasting with models for which specification argument 
#' \code{exogenous_variables} was specified required providing the future values
#' of these exogenous variables in the argument \code{exogenous_forecast} of the
#' \code{\link{forecast.PosteriorBVARPANEL}} function.
#' 
#' @method forecast PosteriorBVARPANEL
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' PosteriorBVARPANEL obtained by running the \code{estimate} function.
#' @param horizon a positive integer, specifying the forecasting horizon.
#' @param exogenous_forecast not used here ATM; included for compatibility with 
#' generic \code{forecast}.
#' @param conditional_forecast a list of length \code{C} containing 
#' \code{horizon x N} matrices with forecasted values for selected variables. 
#' These matrices should only contain \code{numeric} or \code{NA} values. The 
#' entries with \code{NA} values correspond to the values that are forecasted 
#' conditionally on the realisations provided as \code{numeric} values.
#' 
#' @return A list of class \code{ForecastsPANEL} with \code{C} elements containing 
#' the draws from the country-specific predictive density and data in a form of 
#' object class \code{Forecasts} that includes:
#' 
#' \describe{
#'  \item{forecasts}{an \code{horizonxNxS} array with the draws from the 
#'  country-specific predictive density}
#'  \item{Y}{a \code{T_cxN} matrix with the country-specific data}
#' }
#'
#' @references
#' Waggoner, D. F., & Zha, T. (1999) 
#' Conditional forecasts in dynamic multivariate models, 
#' \emph{Review of Economics and Statistics}, \bold{81}(4), 639-651,
#' \doi{10.1162/003465399558508}.
#'
#' @seealso \code{\link{estimate.PosteriorBVARPANEL}}, 
#' \code{\link{summary.ForecastsPANEL}}, \code{\link{plot.ForecastsPANEL}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(ilo_cubic_panel)                                   # load the data
#' data(ilo_exogenous_variables)                           # load the exogenous variables
#' data(ilo_exogenous_forecasts)                           # load the exogenous forecast
#' set.seed(123)
#' 
#' # specify the model
#' specification = specify_bvarPANEL$new(ilo_cubic_panel, exogenous = ilo_exogenous_variables)
#' burn_in       = estimate(specification, 10)             # run the burn-in; use say S = 5000
#' posterior     = estimate(burn_in, 10)                   # estimate the model; use say S = 10000
#' 
#' # forecast 6 years ahead
#' predictive    = forecast(posterior, 6, exogenous_forecast = ilo_exogenous_forecasts)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_cubic_panel |>
#'   specify_bvarPANEL$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20) |> 
#'   forecast(horizon = 2) -> predictive
#' 
#' # conditional forecasting 6 years ahead conditioning on 
#' #  provided future values for the Gross Domestic Product 
#' #  growth rate
#' ############################################################
#' data(ilo_conditional_forecasts)                        # load the conditional forecasts of dgdp
#' specification = specify_bvarPANEL$new(ilo_cubic_panel)    # specify the model
#' burn_in       = estimate(specification, 10)            # run the burn-in; use say S = 5000
#' posterior     = estimate(burn_in, 10)                  # estimate the model; use say S = 10000
#' # forecast 6 years ahead
#' predictive    = forecast(posterior, 6, conditional_forecast = ilo_conditional_forecasts)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_cubic_panel |>
#'   specify_bvarPANEL$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20) |> 
#'   forecast(
#'     horizon = 6, 
#'     conditional_forecast = ilo_conditional_forecasts
#'   ) -> predictive
#' 
#' @export
forecast.PosteriorBVARPANEL = function(
    posterior, 
    horizon = 1, 
    exogenous_forecast = NULL,
    conditional_forecast = NULL
) {
  
  posterior_A_c_cpp     = posterior$posterior$A_c_cpp
  posterior_Sigma_c_cpp = posterior$posterior$Sigma_c_cpp
  X_c             = posterior$last_draw$data_matrices$X
  Y_c             = posterior$last_draw$data_matrices$Y
  N               = dim(Y_c[[1]])[2]
  K               = dim(X_c[[1]])[2]
  C               = length(Y_c)
  c_names         = names(posterior$last_draw$data_matrices$Y)
  
  d               = K - N * posterior$last_draw$p - 1
  if (d == 0 ) {
    # this will not be used for forecasting, but needs to be provided
    exogenous_forecast = list()
    for (c in 1:C) exogenous_forecast[[c]] = matrix(NA, horizon, 1)
  } else {
    stopifnot("Forecasted values of exogenous variables are missing." 
              = (d > 0) & !is.null(exogenous_forecast))
    stopifnot("The matrix of exogenous_forecast does not have a correct number of columns." 
              = unique(simplify2array(lapply(exogenous_forecast, function(x){ncol(x)}))) == d)
    stopifnot("Provide exogenous_forecast for all forecast periods specified by argument horizon." 
              = unique(simplify2array(lapply(exogenous_forecast, function(x){nrow(x)}))) == horizon)
    stopifnot("Argument exogenous has to be a matrix." 
              = all(simplify2array(lapply(exogenous_forecast, function(x){is.matrix(x) & is.numeric(x)}))))
    stopifnot("Argument exogenous cannot include missing values." 
              = unique(simplify2array(lapply(exogenous_forecast, function(x){any(is.na(x))}))) == FALSE)
  }
  
  if ( is.null(conditional_forecast) ) {
    conditional_forecast = list()
    for (c in 1:C) conditional_forecast[[c]] = matrix(NA, horizon, N)
  } else {
    stopifnot("Argument conditional_forecast must be a list with the same countries 
              as in the provided data." 
              = is.list(conditional_forecast) & length(conditional_forecast) == length(Y_c)
    )
    stopifnot("Argument conditional_forecast must be a list with the same countries 
              as in the provided data."
              = all(names(Y_c) == names(conditional_forecast))
    )
    stopifnot("Argument conditional_forecast must be a list with matrices with numeric values."
              = all(sapply(conditional_forecast, function(x) is.matrix(x) & is.numeric(x)))
    )
    stopifnot("All the matrices provided in argument conditional_forecast must have 
              the same number of rows equal to the value of argument horizon."
              = unique(sapply(conditional_forecast, function(x) nrow(x) )) == horizon
    )
    stopifnot("All the matrices provided in argument conditional_forecast must have 
              the same number of columns equal to the number of columns in the used data."
              = unique(sapply(conditional_forecast, function(x) ncol(x) )) == N
    )
  }
  
  # perform forecasting
  fff           = .Call(`_bvarPANELs_forecast_bvarPANEL`, 
                        posterior_A_c_cpp, 
                        posterior_Sigma_c_cpp, 
                        X_c, 
                        conditional_forecast, 
                        exogenous_forecast, 
                        horizon
                       )
                          
  forecasts       = list()
  for (c in 1:C) {
    fore            = list()
    fore_tmp        = aperm(fff$forecasts_cpp[c,1][[1]], c(2,1,3))
    fore$forecasts  = fore_tmp
    fore$Y          = t(Y_c[[c]])
    class(fore)     = "Forecasts"
    forecasts[[c]]  = fore
  }
  names(forecasts)  = c_names
  class(forecasts)  = "ForecastsPANEL"
  
  return(forecasts)
}
