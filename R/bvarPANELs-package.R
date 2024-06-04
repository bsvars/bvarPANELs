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
#' @description Forecasting a multi-country time series panel data using 
#' Bayesian Vector Autoregressions with a three-level country-global 
#' hierarchical prior structure. Copyright: 2024 International Labour Organization.
#' 
#' @name bvarPANELs-package
#' @aliases bvarPANELs-package bvarPANELs
#' @docType package
#' @useDynLib bvarPANELs, .registration = TRUE
#' @importFrom bsvars estimate forecast
#' @importFrom Rcpp sourceCpp
#' @importFrom R6 R6Class
#' @importFrom RcppTN rtn dtn
#' @import RcppProgress
#' @note This package is currently in active development.
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' @keywords package models ts
#' #' @examples
#' @examples
#' data(ilo_cubic_panel)                                   # load the data
#' set.seed(123)
#' specification = specify_bvarPANEL$new(ilo_cubic_panel)  # specify the model
#' burn_in       = estimate(specification, 20)             # run the burn-in
#' posterior     = estimate(burn_in, 20)                   # estimate the model
#' predictive    = forecast(posterior, 2)                  # forecast 2 years ahead
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_cubic_panel |>
#'   specify_bvarPANEL$new() |>
#'   estimate(S = 20) |> 
#'   estimate(S = 20) |> 
#'   forecast(horizon = 2) -> predictive
"_PACKAGE"
