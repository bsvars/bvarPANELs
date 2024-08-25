
#' @title Bayesian estimation of a Bayesian Hierarchical Panel Vector 
#' Autoregression using Gibbs sampler
#'
#' @description Estimates the Bayesian Hierarchical Panel VAR using the Gibbs 
#' sampler proposed by Sanchez-Martinez & Woźniak (2024).
#' 
#' @details 
#' The Bayesian Hierarchical Panel Vector Autoregressive model described in 
#' \code{\link{bvarPANELs}} is estimated using the Gibbs sampler. In this 
#' estimation procedure all the parameters of the model are estimated jointly.
#' The list of parameters of the model includes:
#' \describe{
#'  \item{\eqn{\mathbf{A}_c}}{a \code{KxN} country-specific autoregressive parameters matrix}
#'  \item{\eqn{\mathbf{\Sigma}_c}}{an \code{NxN} country-specific covariance matrix}
#'  \item{\eqn{\mathbf{A}}}{a \code{KxN} global autoregressive parameters matrix}
#'  \item{\eqn{\mathbf{\Sigma}}}{an \code{NxN} global covariance matrix}
#'  \item{\eqn{\mathbf{V}}}{a \code{KxK} covariance matrix of prior for global autoregressive parameters}
#'  \item{\eqn{\nu}}{prior degrees of freedom parameter}
#'  \item{\eqn{m}}{prior average global persistence parameter}
#'  \item{\eqn{w}}{prior scaling parameter}
#'  \item{\eqn{s}}{prior scaling parameter}
#' }
#' 
#' \strong{Gibbs sampler.}
#' Is an algorithm to sample random draws from the posterior distribution of the
#' parameters of the model given the data. The algorithm is briefly explained 
#' on an example of two-parameter model with parameters \eqn{\theta_1} and 
#' \eqn{\theta_2}. In order to sample from the joint posterior distribution 
#' \eqn{p(\theta_1,\theta_2|\mathbf{Y})} the Gibbs sampler proceeds by sampling 
#' from full-conditional posterior distributions of each parameter given all the
#' other parameters, denoted by \eqn{p(\theta_1|\theta_2,\mathbf{Y})} and 
#' \eqn{p(\theta_2|\theta_1,\mathbf{Y})}.
#' 
#' To obtain \code{S} draws from the posterior distribution:
#' \enumerate{
#' \item Set the initial values of the parameters \eqn{\theta_2^{(0)}}
#' \item At each of the \code{s} iteration:
#' \enumerate{
#' \item Sample \eqn{\theta_1^{(s)}} from \eqn{p(\theta_1|\theta_2^{(s-1)},\mathbf{Y})}
#' \item Sample \eqn{\theta_2^{(s)}} from \eqn{p(\theta_2|\theta_1^{(s)},\mathbf{Y})}
#' } 
#' \item Repeat step 2. \code{S} times. Return \eqn{\{\theta_1^{(s)},\theta_2^{(s)}\}_{s=1}^{S}} 
#' as a sample drawn from the posterior distribution \eqn{p(\theta_1,\theta_2|\mathbf{Y})}.
#' }
#' The \code{estimate()} function returns the draws from the posetrior distribution
#' of the parameters of the hierarchical panel VAR model listed above.
#' 
#' \strong{Thinning.} 
#' Thinning is a procedure to reduce the dependence in the returned sample from
#' the posterior distribution. It is obtained by returning every \code{thin} 
#' draw in the final sample. This procedure reduces the number of draws returned 
#' by the \code{estimate()} function.
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
#' @seealso \code{\link{bvarPANELs}}, \code{\link{specify_bvarPANEL}}, \code{\link{specify_posterior_bvarPANEL}}
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
