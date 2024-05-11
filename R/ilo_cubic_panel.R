
#' @title A 4-variable annual system for forecasting labour market outcomes 
#' for 189 United Nations countries from 1991 to 2023
#'
#' @description For each of the countries a time series of 33 observations on 4 
#' variables including the unemployment rate (UR), employment rate (EPR), labour
#' force participation rate (LFPR) and the growth rate of GDP (dgdp) is provided.
#' The missing observations are filled using imputation method.
#' Last data update was implemented on 2024-05-11.
#'
#' @usage data(ilo_cubic_panel)
#' 
#' @format A list of 189 \code{ts} objects with time series of 33 observations 
#' on 4 variables:
#' \describe{
#'   \item{UR}{annual unemployment rate}
#'   \item{EPR}{annual employment rate}
#'   \item{LFPR}{annual labour force participation rate}
#'   \item{dgdp}{annual growth rate of gross domestic product}
#' }
#' 
#' @source 
#' International Labour Organization. (2020). ILO modelled estimates database, 
#' ILOSTAT [database]. Available from \url{https://ilostat.ilo.org/data/}.
#' 
#' @examples 
#' data(ilo_cubic_panel)   # upload the data
#' 
"ilo_cubic_panel"