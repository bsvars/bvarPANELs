
#' @title Data containing future observations for 189 United Nations countries 
#' from 2024 to 2029 to be used to forecast with models with 
#' \code{ilo_exogenous_variables}
#'
#' @description For each of the countries a time series of 6 observations on 
#' On the dummies is provided. These future values are all equal to zero. They 
#' provide benchmark for the objects to be used when \code{exogenous_variables} 
#' are used. 
#' Last data update was implemented on 2024-05-11.
#'
#' @usage data(ilo_exogenous_forecasts)
#' 
#' @format A list of 189 \code{ts} objects with time series of 6 observations 
#' on 3 variables:
#' \describe{
#'   \item{2008}{the aftermath of the Global Financial Crisis}
#'   \item{2020}{the COVID pandemic}
#'   \item{2021}{the aftermath of the COVID pandemic}
#' }
#' 
#' @examples 
#' data(ilo_exogenous_forecasts)   # upload the data
#' 
"ilo_exogenous_forecasts"