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
#' @param conditional_forecast a list of length \code{C} containing 
#' \code{horizon x N} matrices with forecasted values for selected variables. 
#' These matrices should only contain \code{numeric} or \code{NA} values. The 
#' entries with \code{NA} values correspond to the values that are forecasted 
#' conditionally on the realisations provided as \code{numeric} values.
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
#' # conditional forecasting 6 years ahead conditioning on 
#' #  provided future values for the Gross Domestic Product 
#' #  growth rate
#' ############################################################
#' #' data(ilo_conditional_forecast)                        # load the conditional forecasts of dgdp
#' predictive    = forecast(posterior, 6, conditional_forecast = ilo_conditional_forecast)
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
#'     conditional_forecast = ilo_conditional_forecast
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
  
  do_conditional_forecasting = !is.null(conditional_forecast)
  
  if (!do_conditional_forecasting) {
    # perform forecasting
    fore          = .Call(`_bvarPANELs_forecast_bvarPANEL`, 
                          posterior_A_c_cpp, posterior_Sigma_c_cpp, X_c, horizon)
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
    
    # perform conditional forecasting
    fore          = .Call(`_bvarPANELs_forecast_conditional_bvarPANEL`, 
                          posterior_A_c_cpp, posterior_Sigma_c_cpp, X_c, conditional_forecast, horizon)
  }
  
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
