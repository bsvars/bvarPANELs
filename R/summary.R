

#' @title Provides posterior estimation summary for Bayesian Hierarchical Panel 
#' Vector Autoregressions
#'
#' @description Provides posterior mean, standard deviations, as well as 5 and 95 
#' percentiles of the parameters for all \code{C} countries.
#' 
#' @param object an object of class \code{PosteriorBVARPANEL} obtained using the
#' \code{estimate()} function applied to  
#' Vector Autoregressions containing draws from the  posterior distribution of 
#' the parameters. 
#' @param ... additional arguments affecting the summary produced.
#' 
#' @return A list reporting the posterior mean, standard deviations, as well as 5 and 95 
#' percentiles of the country-specific and global parameters.
#' 
#' @method summary PosteriorBVARPANEL
#' 
#' @seealso \code{\link{estimate.BVARPANEL}}, \code{\link{specify_bvarPANEL}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' # upload data
#' data(ilo_cubic_panel)                                   # load the data
#' data(ilo_exogenous_variables)                           # load the exogenous variables
#' 
#' set.seed(123)
#' 
#' # specify the model
#' specification = specify_bvarPANEL$new(ilo_cubic_panel, exogenous = ilo_exogenous_variables)
#' burn_in       = estimate(specification, 10)             # run the burn-in
#' posterior     = estimate(burn_in, 10)                   # estimate the model
#' summary(posterior)
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_cubic_panel |>
#'   specify_bvarPANEL$new(exogenous = ilo_exogenous_variables) |>
#'   estimate(S = 10) |> 
#'   estimate(S = 10) |> 
#'   summary()
#' 
#' @export
summary.PosteriorBVARPANEL = function(
    object,
    ...
) {
  
  S         = dim(object$posterior$A_c)[4]
  C         = dim(object$posterior$A_c)[3]
  N         = dim(object$posterior$A_c)[2]
  K         = dim(object$posterior$A_c)[1]
  p         = object$last_draw$p
  d         = K - N * p
  
  out       = list()
  param     = c("A", "Sigma")
  country_names = names(object$last_draw$data_matrices$Y)
  
  # country-specific parameter summary
  for (c in 1:C) {
    
    out_c         = list()
    out_c$Sigma   = list()
    out_c$A       = list()
    
    for (n in 1:N) {
    
      Sigma_c       = matrix(object$posterior$Sigma_c[n,1:n,c,], nrow = n)
      out_c$Sigma[[n]] = matrix(
        cbind(
          apply(Sigma_c, 1, mean),
          apply(Sigma_c, 1, sd),
          apply(Sigma_c, 1, quantile, probs = 0.05),
          apply(Sigma_c, 1, quantile, probs = 0.95)
        ),
        ncol = 4
      )
      colnames(out_c$Sigma[[n]]) = c("mean", "sd", "5% quantile", "95% quantile")
      rownames(out_c$Sigma[[n]]) = paste0("Sigma[", n, ",", 1:n, "]")  
      
      A_c      = object$posterior$A_c[,n,c,]  
      out_c$A[[n]] = cbind(
        apply(A_c, 1, mean),
        apply(A_c, 1, sd),
        apply(A_c, 1, quantile, probs = 0.05),
        apply(A_c, 1, quantile, probs = 0.95)
      )
      colnames(out_c$A[[n]]) = c("mean", "sd", "5% quantile", "95% quantile")
      
      Anames  = c(
        paste0(
          rep("lag", p * N),
          kronecker((1:p), rep(1, N)),
          rep("_var", p * N),
          kronecker((1:N), rep(1, p))
        ),
        "const"
      )
      if (d > 1) {
        Anames = c(Anames, paste0("exo", 1:(d - 1)))
      }
      rownames(out_c$A[[n]]) = Anames
    } # END n loop
    
    names(out_c$Sigma) = paste0("equation", 1:N)
    names(out_c$A) = paste0("equation", 1:N)
    
    out[[c]]  = out_c
  } # END c loop
  
  names(out) = country_names
  
  
  # global parameter summary
  out_g         = list()
  out_g$A       = list()
  out_g$Sigma   = list()
  out_g$V       = list()
  out_g$hyper   = list()
  
  for (n in 1:N) {
    
    Sigma             = matrix(object$posterior$Sigma[n,1:n,], nrow = n)
    out_g$Sigma[[n]]  = matrix(
      cbind(
        apply(Sigma, 1, mean),
        apply(Sigma, 1, sd),
        apply(Sigma, 1, quantile, probs = 0.05),
        apply(Sigma, 1, quantile, probs = 0.95)
      ),
      ncol = 4
    )
    colnames(out_g$Sigma[[n]]) = c("mean", "sd", "5% quantile", "95% quantile")
    rownames(out_g$Sigma[[n]]) = paste0("Sigma[", n, ",", 1:n, "]")  
    
    A      = object$posterior$A[,n,]  
    out_g$A[[n]] = cbind(
      apply(A, 1, mean),
      apply(A, 1, sd),
      apply(A, 1, quantile, probs = 0.05),
      apply(A, 1, quantile, probs = 0.95)
    )
    colnames(out_g$A[[n]]) = c("mean", "sd", "5% quantile", "95% quantile")
    rownames(out_g$A[[n]]) = Anames
    
  } # END n loop
  
  names(out_g$Sigma)  = paste0("equation", 1:N)
  names(out_g$A)      = paste0("equation", 1:N)
  
  
  for (k in 1:K) {
    
    V             = matrix(object$posterior$V[k,1:k,], nrow = k)
    out_g$V[[k]]  = matrix(
      cbind(
        apply(V, 1, mean),
        apply(V, 1, sd),
        apply(V, 1, quantile, probs = 0.05),
        apply(V, 1, quantile, probs = 0.95)
      ),
      ncol = 4
    )
    colnames(out_g$V[[k]]) = c("mean", "sd", "5% quantile", "95% quantile")
    rownames(out_g$V[[k]]) = paste0("V[", k, ",", 1:k, "]")  
    
  } # END k loop
  
  hyper         = t(cbind(
    object$posterior$nu,
    object$posterior$m,
    object$posterior$w,
    object$posterior$s
  ))
  out_g$hyper  = matrix(
    cbind(
      apply(hyper, 1, mean),
      apply(hyper, 1, sd),
      apply(hyper, 1, quantile, probs = 0.05),
      apply(hyper, 1, quantile, probs = 0.95)
    ),
    ncol = 4
  )
  colnames(out_g$hyper) = c("mean", "sd", "5% quantile", "95% quantile")
  rownames(out_g$hyper) = c("nu", "m", "w", "s")
  
  out$global = out_g
  
  return(out)
} # END summary.PosteriorBVARPANEL






#' @title Provides posterior summary of forecast error variance decompositions
#'
#' @description Provides posterior means of the forecast error variance 
#' decompositions of each variable at all horizons.
#' 
#' @param object an object of class \code{PosteriorFEVDPANEL} obtained using the
#' \code{compute_variance_decompositions()} function containing draws from the 
#' posterior distribution of the forecast error variance decompositions. 
#' @param which_c a positive integer or a character string specifying the country 
#' for which the forecast should be plotted.
#' @param ... additional arguments affecting the summary produced.
#' 
#' @return A list reporting the posterior mean of the forecast error variance 
#' decompositions of each variable at all horizons.
#' 
#' @method summary PosteriorFEVDPANEL
#' 
#' @seealso \code{\link{compute_variance_decompositions.PosteriorBVARPANEL}}, \code{\link{plot}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
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
#' summary(fevd, which_c = "POL")
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_cubic_panel |>
#'   specify_bvarPANEL$new(p = 1) |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20) |> 
#'   compute_variance_decompositions(horizon = 4) |> 
#'   summary(which_c = "global")
#' 
#' @export
summary.PosteriorFEVDPANEL = function(
    object,
    which_c,
    ...
) {
  
  if (is.numeric(which_c)) {
    stopifnot("Argument which_c must be a positive integer indicating one of the countries."
              = length(which_c) == 1 & which_c %% 1 == 0 & which_c > 0 & which_c <= length(object))
  } else if (is.character(which_c)) {
    stopifnot("Argument which_c must be a character string indicating one of the countries."
              = which_c %in% names(object))
  } else {
    stop("Argument which_c must be either a positive integer or a character string.")
  }
  
  summary(object[[which_c]], ...)
}





#' @title Provides posterior summary of country-specific Forecasts
#'
#' @description Provides posterior summary of the forecasts including their 
#' mean, standard deviations, as well as 5 and 95 percentiles.
#' 
#' @param object an object of class \code{ForecastsPANEL} obtained using the
#' \code{forecast()} function containing draws the predictive density. 
#' @param which_c a positive integer or a character string specifying the country 
#' for which the forecast should be plotted.
#' @param ... additional arguments affecting the summary produced.
#' 
#' @return A list reporting the posterior mean, standard deviations, as well as 
#' 5 and 95 percentiles of the forecasts for each of the variables and forecast 
#' horizons.
#' 
#' @method summary ForecastsPANEL
#' 
#' @seealso \code{\link{forecast.PosteriorBVARPANEL}}, \code{\link{plot}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' data(ilo_cubic_panel)                                   # load the data
#' data(ilo_exogenous_variables)                           # load the exogenous variables
#' data(ilo_exogenous_forecasts)                           # load the exogenous forecast
#' set.seed(123)
#' 
#' # specify the model
#' specification = specify_bvarPANEL$new(ilo_cubic_panel, exogenous = ilo_exogenous_variables)
#' burn_in       = estimate(specification, 10)             # run the burn-in
#' posterior     = estimate(burn_in, 10)                   # estimate the model
#' 
#' # forecast 6 years ahead
#' predictive    = forecast(posterior, 6, exogenous_forecast = ilo_exogenous_forecasts)
#' summary(predictive, which_c = "POL")
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_cubic_panel |>
#'   specify_bvarPANEL$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20) |> 
#'   forecast(horizon = 2) |> 
#'   summary(which_c = "POL")
#' 
#' # conditional forecasting 6 years ahead conditioning on 
#' #  provided future values for the Gross Domestic Product 
#' #  growth rate
#' ############################################################
#' data(ilo_conditional_forecasts)                        # load the conditional forecasts of dgdp
#' specification = specify_bvarPANEL$new(ilo_cubic_panel)    # specify the model
#' burn_in       = estimate(specification, 10)               # run the burn-in
#' posterior     = estimate(burn_in, 10)                     # estimate the model
#' # forecast 6 years ahead
#' predictive    = forecast(posterior, 6, conditional_forecast = ilo_conditional_forecasts)
#' summary(predictive, which_c = "POL")
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_cubic_panel |>
#'   specify_bvarPANEL$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20) |> 
#'   forecast(
#'     horizon = 6, 
#'     conditional_forecast = ilo_conditional_forecasts
#'   ) |> 
#'   summary(which_c = "POL")
#' 
#' @export
summary.ForecastsPANEL = function(
    object,
    which_c,
    ...
) {
  
  if (is.numeric(which_c)) {
    stopifnot("Argument which_c must be a positive integer indicating one of the countries."
              = length(which_c) == 1 & which_c %% 1 == 0 & which_c > 0 & which_c <= length(object))
  } else if (is.character(which_c)) {
    stopifnot("Argument which_c must be a character string indicating one of the countries."
              = which_c %in% names(object))
  } else {
    stop("Argument which_c must be either a positive integer or a character string.")
  }
  
  summary(object[[which_c]], ...)  
} # END summary.ForecastsPANEL
