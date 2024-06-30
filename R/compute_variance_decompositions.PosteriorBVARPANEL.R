




#' @title Computes posterior draws of the forecast error variance decomposition
#' @description For each country, each of the draws from the posterior estimation 
#' of the model is transformed into a draw from the posterior distribution of the forecast 
#' error variance decomposition.
#' 
#' @method compute_variance_decompositions PosteriorBVARPANEL
#' 
#' @param posterior posterior estimation outcome - an object of class 
#' \code{PosteriorBVARPANEL} obtained by running the \code{estimate} function.
#' 
#' @param horizon a positive integer number denoting the forecast horizon for 
#' the forecast error variance decompositions.
#' 
#' @return  An object of class PosteriorFEVDPANEL, that is, an 
#' \code{CxNxNx(horizon+1)xS} array with attribute PosteriorFEVDPANEL 
#' containing \code{S} draws of the forecast error variance decompositions for 
#' each country.
#' 
#' @seealso \code{\link{estimate}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
#' 
#' @examples
#' # upload data
#' data(ilo_cubic_panel)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' specification  = specify_bvarPANEL$new(ilo_cubic_panel, p = 1)
#' 
#' # run the burn-in
#' burn_in        = estimate(specification, 10)
#' 
#' # estimate the model
#' posterior      = estimate(burn_in, 20)
#' 
#' # compute forecast error variance decomposition 4 years ahead
#' fevd           = compute_variance_decompositions(posterior, horizon = 4)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_cubic_panel |>
#'   specify_bvarPANEL$new(p = 1) |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20) |> 
#'   compute_variance_decompositions(horizon = 4) -> fevd
#' 
#' @export
compute_variance_decompositions.PosteriorBVARPANEL <- function(posterior, horizon) {

  posterior_Sigma = posterior$posterior$Sigma_c
  posterior_A     = posterior$posterior$A_c
  posterior_Sg    = posterior$posterior$Sigma
  posterior_Ag    = posterior$posterior$A
  N               = dim(posterior_A)[2]
  C               = dim(posterior_A)[3]
  S               = dim(posterior_A)[4]
  p               = posterior$last_draw$p
  
  fff             = .Call(`_bvarPANELs_panel_variance_decompositions`, posterior_Sigma, posterior_A, posterior_Sg, posterior_Ag, horizon, p, TRUE)
  fevd            = array(NA, c(C + 1, N, N, horizon + 1, S))
  
  for (s in 1:S) {
    for (c in 1:(C + 1)) {
      fevd[c,,,,s] = fff[c,s][[1]]
    }
  }
  
  class(fevd) <- "PosteriorFEVDPANEL"
}