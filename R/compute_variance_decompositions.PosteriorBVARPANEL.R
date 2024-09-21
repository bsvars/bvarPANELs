
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
#' @return  An object of class \code{PosteriorFEVDPANEL}, that is, a list with 
#' \code{C} elements containing \code{NxNx(horizon+1)xS} arrays of class 
#' \code{PosteriorFEVD} with \code{S} draws of country-specific forecast error 
#' variance decompositions.
#' 
#' @seealso \code{\link{estimate.PosteriorBVARPANEL}}, 
#' \code{\link{summary.PosteriorFEVDPANEL}}, 
#' \code{\link{plot.PosteriorFEVDPANEL}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @references 
#' Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
#' 
#' @examples
#' # upload data
#' data(ilo_dynamic_panel)
#' 
#' # specify the model and set seed
#' set.seed(123)
#' specification  = specify_bvarPANEL$new(ilo_dynamic_panel, p = 1)
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
#' ilo_dynamic_panel |>
#'   specify_bvarPANEL$new(p = 1) |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20) |> 
#'   compute_variance_decompositions(horizon = 4) -> fevd
#' 
#' @export
compute_variance_decompositions.PosteriorBVARPANEL <- function(posterior, horizon) {

  posterior_Sigma = posterior$posterior$Sigma_c_cpp
  posterior_A     = posterior$posterior$A_c_cpp
  posterior_Sg    = posterior$posterior$Sigma
  posterior_Ag    = posterior$posterior$A
  N               = dim(posterior_Ag)[2]
  C               = dim(posterior_A[1][[1]])[3]
  S               = dim(posterior_A)[1]
  p               = posterior$last_draw$p
  c_names         = names(posterior$last_draw$data_matrices$Y)
  
  fff             = .Call(`_bvarPANELs_panel_variance_decompositions`, 
                          posterior_Sigma, 
                          posterior_A, 
                          posterior_Sg, 
                          posterior_Ag, 
                          horizon, 
                          p, 
                          TRUE
                    )
  
  fevd            = list()
  for (c in 1:(C + 1)) {
    fevd_c          = array(NA, c(N, N, horizon + 1, S))
    for (s in 1:S) {
      fevd_c[,,,s]  = fff[c, s][[1]]
    }
    na_check        = apply(fevd_c, 4, function(x) any(is.na(x)))
    fevd_c          = fevd_c[,,, !na_check]
    class(fevd_c)   = "PosteriorFEVD"
    fevd[[c]]       = fevd_c
  }
  
  names(fevd) = c(c_names, "global")
  
  class(fevd) <- "PosteriorFEVDPANEL"
  return(fevd)
}