
#' @title A 3-variable annual system for of dummy observations for 2008, 2020, 
#' and 2021 to be used in the estimation of the Panel VAR model for 
#' 189 United Nations countries from 1991 to 2023
#'
#' @description For each of the countries a time series of 33 observations on 3 
#' dummy variables for the years 2008, 2020, and 2021 is provided. 
#' Last data update was implemented on 2024-06-29.
#'
#' @usage data(ilo_exogenous_variables)
#' 
#' @format A list of 189 \code{ts} objects with time series of 33 observations 
#' on 3 variables:
#' \describe{
#'   \item{2008}{the aftermath of the Global Financial Crisis}
#'   \item{2020}{the COVID pandemic}
#'   \item{2021}{the aftermath of the COVID pandemic}
#' }
#' 
#' @examples 
#' data(ilo_exogenous_variables)   # upload the data
#' 
"ilo_exogenous_variables"
