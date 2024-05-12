#' @title Forecasting using Hierarchical Pannel Vector Autoregressions
#'
#' @description Samples from the joint predictive density of the dependent 
#' variables for all countries at forecast horizons 
#' from 1 to \code{horizon} specified as an argument of the function.
#' 
#' @method forecast PosteriorBVARPANEL
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' PosteriorBVARPANEL obtained by running the \code{estimate} function.
#' @param horizon a positive integer, specifying the forecasting horizon.
#' @param exogenous_forecast not used here ATM; included for compatibility with 
#' generic \code{forecast}.
#' 
#' @return A list of class \code{PanelForecasts} containing the
#' draws from the predictive density and data. The output list includes element:
#' 
#' \describe{
#'  \item{forecasts}{an \code{horizonxNxCxS} array with the draws from predictive density}
#'  \item{forecasts_cpp}{an unspecified object passed for computations in **cpp**}
#'  \item{Y}{a \code{C}-element list with \code{T_cxN} matrices with the country-specific data}
#' }
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(ilo_cubic_panel)                                   # load the data
#' set.seed(123)
#' specification = specify_bvarPANEL$new(ilo_cubic_panel)  # specify the model
#' burn_in       = estimate(specification, 10)             # run the burn-in
#' posterior     = estimate(burn_in, 10)                   # estimate the model
#' predictive    = forecast(posterior, 2)                  # forecast 2 years ahead
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
#' @export
forecast.PosteriorBVARPANEL = function(posterior, horizon = 1, exogenous_forecast = NULL) {
  
  
  posterior_A_c_cpp     = posterior$posterior$A_c_cpp
  posterior_Sigma_c_cpp = posterior$posterior$Sigma_c_cpp
  X_c             = posterior$last_draw$data_matrices$X
  Y_c             = posterior$last_draw$data_matrices$Y
  
  fore            = .Call(`_bvarPANELs_forecast_bvarPANEL`, posterior_A_c_cpp, posterior_Sigma_c_cpp, X_c, horizon)
  
  N               = dim(Y_c[[1]])[2]
  S               = dim(posterior_A_c_cpp)[1]
  C               = length(Y_c)
  forecasts      = array(NA, c(horizon, N, S, C))
  for (c in 1:C) {
    forecasts[,,,c] = fore$forecasts_cpp[c,1][[1]]
  }
  forecasts       = aperm(forecasts, c(1,2,4,3))
  fore$forecasts  = forecasts
  fore$Y          = Y_c
  class(fore)     = "PanelForecasts"
  
  return(fore)
}
