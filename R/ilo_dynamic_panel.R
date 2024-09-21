
#' @title A 4-variable annual system for forecasting labour market outcomes 
#' for 189 United Nations countries from 1991 to 2023
#'
#' @description For each of the countries a time series of 33 observations on 4 
#' variables including the logarithm of Gross Domestic Product (gdp), as well as
#' the labour market outcomes including the unemployment rate (UR), employment 
#' rate (EPR), labour force participation rate (LFPR). The missing observations 
#' are filled using imputation method.
#' Last data update was implemented on 2024-09-21.
#'
#' @usage data(ilo_dynamic_panel)
#' 
#' @format A list of 189 \code{ts} objects with time series of 33 observations 
#' on 4 variables:
#' \describe{
#'   \item{gdp}{logarithm of gross domestic product}
#'   \item{UR}{annual unemployment rate}
#'   \item{EPR}{annual employment rate}
#'   \item{LFPR}{annual labour force participation rate}
#' }
#' 
#' @source 
#' International Labour Organization. (2020). ILO modelled estimates database, 
#' ILOSTAT [database]. Available from \url{https://ilostat.ilo.org/data/}.
#' 
#' @examples 
#' data(ilo_dynamic_panel)   # upload the data
#' 
"ilo_dynamic_panel"