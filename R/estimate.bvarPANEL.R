
#' @title Bayesian estimation of a Bayesian Hierarchical Panel Vector 
#' Autoregression using Gibbs sampler
#'
#' @description Estimates the Bayesian Hierarchical Panel VAR using the Gibbs 
#' sampler proposed by Sanchez-Martinez & Woźniak (2024).
#' 
#' @details 
#' The homoskedastic SVAR model is given by the reduced form equation:
#' \deqn{Y_c = A_cX_c + E_c}
#' where \eqn{Y_c} is an \code{T_c x N} matrix of dependent variables for 
#' country \code{c}, \eqn{X_c} is a \code{T_c x K} matrix of explanatory 
#' variables, \eqn{E_c} is an \code{T_c x N} matrix of error terms, and 
#' \eqn{A_c} is an \code{NxK} matrix of country-specific autoregressive slope 
#' coefficients and parameters on deterministic terms in \eqn{X_c}.
#' 
#' @param specification an object of class \code{BVARPANEL} generated using the 
#' \code{specify_bvarPANEL$new()} function.
#' @param S a positive integer, the number of posterior draws to be generated
#' @param thin a positive integer, specifying the frequency of MCMC output thinning
#' @param show_progress a logical value, if \code{TRUE} the estimation progress 
#' bar is visible
#' 
#' @return An object of class \code{PosteriorBVARPANEL} containing the Bayesian 
#' estimation output and containing two elements:
#' 
#'  \code{posterior} a list with a collection of \code{S} draws from the 
#'  posterior distribution generated via Gibbs sampler.
#' \code{last_draw} an object of class \code{BVARPANEL} with the last draw of the 
#' current MCMC run as the starting value to be passed to the continuation of 
#' the MCMC estimation using the \code{estimate()} method. 
#'
#' @seealso \code{\link{specify_bvarPANEL}}, \code{\link{specify_posterior_bvarPANEL}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @method estimate BVARPANEL
#' 
#' @examples
#' data(ilo_cubic_panel)                                   # load the data
#' data(ilo_exogenous_variables)                           # load the exogenous variables
#' set.seed(123)
#' # specify the model
#' specification = specify_bvarPANEL$new(ilo_cubic_panel, exogenous = ilo_exogenous_variables)
#' burn_in       = estimate(specification, 10)             # run the burn-in
#' posterior     = estimate(burn_in, 10)                   # estimate the model
#' 
#' @export
estimate.BVARPANEL <- function(
    specification, 
    S, 
    thin = 1L, 
    show_progress = TRUE
) {
  
  # get the inputs to estimation
  prior               = specification$prior$get_prior()
  starting_values     = specification$starting_values$get_starting_values()
  data_matrices       = specification$data_matrices$get_data_matrices()
  adaptiveMH          = specification$adaptiveMH
  
  # estimation
  qqq                 = .Call(`_bvarPANELs_bvarPANEL`, S, data_matrices$Y, data_matrices$X, prior, starting_values, thin, show_progress, adaptiveMH)
  
  specification$starting_values$set_starting_values(qqq$last_draw)
  output              = specify_posterior_bvarPANEL$new(specification, qqq$posterior)
  
  return(output)
} # END estimate.BVARPANEL



#' @inherit estimate.BVARPANEL
#' 
#' @method estimate PosteriorBVARPANEL
#' 
#' @param specification an object of class \code{PosteriorBVARPANEL} generated using 
#' the \code{estimate.BVARPANEL()} function. This setup facilitates the 
#' continuation of the MCMC sampling starting from the last draw of the previous 
#' run.
#' 
#' @export
estimate.PosteriorBVARPANEL <- function(
    specification, 
    S, 
    thin = 1, 
    show_progress = TRUE
) {
  
  # get the inputs to estimation
  prior               = specification$last_draw$prior$get_prior()
  starting_values     = specification$last_draw$starting_values$get_starting_values()
  data_matrices       = specification$last_draw$data_matrices$get_data_matrices()
  adaptiveMH          = specification$last_draw$adaptiveMH
  
  # estimation
  qqq                 = .Call(`_bvarPANELs_bvarPANEL`, S, data_matrices$Y, data_matrices$X, prior, starting_values, thin, show_progress, adaptiveMH)
  
  specification$last_draw$starting_values$set_starting_values(qqq$last_draw)
  output              = specify_posterior_bvarPANEL$new(specification$last_draw, qqq$posterior)
  
  return(output)
} # END estimate.PosteriorBSVAR
