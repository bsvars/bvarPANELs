
#' @title Data containing conditional projections for the logarithm of GDP (gdp)
#' for 189 United Nations countries from 2024 to 2029
#'
#' @description For each of the countries a time series of 6 observations on 
#' GDP growth rates (gdp) formatted so they is provided to generate 
#' conditional forecasts of labour market outcomes given the provided projected 
#' paths of output. 
#' Last data update was implemented on 2024-09-21.
#'
#' @usage data(ilo_conditional_forecasts)
#' 
#' @format A list of 189 \code{ts} objects with time series of 6 observations 
#' on 4 variables:
#' \describe{
#'   \item{gdp}{logarithm of gross domestic product - contains projected values}
#'   \item{UR}{unemployment rate - contains missing values}
#'   \item{EPR}{annual employment rate - contains missing values}
#'   \item{LFPR}{annual labour force participation rate - contains missing values}
#' }
#' 
#' @source
#' International Labour Organization. (2020). ILO modelled estimates database,
#' ILOSTAT [database]. Available from \url{https://ilostat.ilo.org/data/}.
#' 
#' @examples 
#' data(ilo_conditional_forecasts)   # upload the data
#' 
"ilo_conditional_forecasts"