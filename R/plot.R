
#' @title Plots fitted values of dependent variables
#'
#' @description Plots of fitted values of dependent variables including their 
#' median and percentiles.
#' 
#' @param x an object of class \code{ForecastsPANEL} obtained using the
#' \code{forecast()} function containing posterior draws of 
#' fitted values of dependent variables.
#' @param which_c a positive integer or a character string specifying the country 
#' for which the forecast should be plotted.
#' @param probability a parameter determining the interval to be plotted. The 
#' interval stretches from the \code{0.5 * (1 - probability)} to 
#' \code{1 - 0.5 * (1 - probability)} percentile of the posterior distribution.
#' @param data_in_plot a fraction value in the range (0, 1) determining how many
#' of the last observations in the data should be plotted with the forecasts.
#' @param col a colour of the plot line and the ribbon
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot ForecastsPANEL
#' 
#' @seealso \code{\link{forecast.PosteriorBVARPANEL}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' specification = specify_bvarPANEL$new(ilo_dynamic_panel)  # specify the model
#' burn_in       = estimate(specification, 10)               # run the burn-in
#' posterior     = estimate(burn_in, 10)                     # estimate the model
#'
#' # forecast 6 years ahead
#' predictive    = forecast(posterior, 6, conditional_forecast = ilo_conditional_forecasts)
#' plot(predictive, which_c = "POL")                                # plot forecasts
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_dynamic_panel |>
#'   specify_bvarPANEL$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 10) |> 
#'   forecast(horizon = 6, conditional_forecast = ilo_conditional_forecasts) |>
#'   plot(which_c = 135)
#' 
#' @export
plot.ForecastsPANEL = function(
    x,
    which_c,
    probability = 0.9,
    data_in_plot = 1,
    col = "#ff69b4",
    main,
    xlab,
    mar.multi = c(1, 4.6, 0, 2.1),
    oma.multi = c(6, 0, 5, 0),
    ...
) {
  
  if (is.numeric(which_c)) {
    stopifnot("Argument which_c must be a positive integer indicating one of the countries."
              = length(which_c) == 1 & which_c %% 1 == 0 & which_c > 0 & which_c <= length(x))
  } else if (is.character(which_c)) {
    stopifnot("Argument which_c must be a character string indicating one of the countries."
              = which_c %in% names(x))
  } else {
    stop("Argument which_c must be either a positive integer or a character string.")
  }
  
  plot(
    x[[which_c]],
    probability,
    data_in_plot,
    col,
    main,
    xlab,
    mar.multi,
    oma.multi,
    ...
  )
} # END plot.ForecastsPANEL








#' @title Plots forecast error variance decompositions
#'
#' @description Plots of the posterior means of the forecast error variance 
#' decompositions.
#' 
#' @param x an object of class \code{PosteriorFEVDPANEL} obtained using the
#' \code{compute_variance_decompositions()} function containing posterior draws of 
#' forecast error variance decompositions.
#' @param which_c a positive integer or a character string specifying the country 
#' for which the forecast should be plotted.
#' @param cols an \code{N}-vector with colours of the plot
#' @param main an alternative main title for the plot
#' @param xlab an alternative x-axis label for the plot
#' @param mar.multi the default \code{mar} argument setting in \code{graphics::par}. Modify with care!
#' @param oma.multi the default \code{oma} argument setting in \code{graphics::par}. Modify with care!
#' @param ... additional arguments affecting the summary produced.
#' 
#' @method plot PosteriorFEVDPANEL
#' 
#' @seealso \code{\link{compute_variance_decompositions.PosteriorBVARPANEL}}
#'
#' @author Tomasz Woźniak \email{wozniak.tom@pm.me}
#' 
#' @examples
#' set.seed(123)
#' specification  = specify_bvarPANEL$new(ilo_dynamic_panel)
#' 
#' # run the burn-in
#' burn_in        = estimate(specification, 10)
#' 
#' # estimate the model
#' posterior      = estimate(burn_in, 20)
#' 
#' # compute forecast error variance decomposition 4 years ahead
#' fevd           = compute_variance_decompositions(posterior, horizon = 4)
#' plot(fevd, which_c = "POL")
#' 
#' # workflow with the pipe |>
#' ############################################################
#' set.seed(123)
#' ilo_dynamic_panel |>
#'   specify_bvarPANEL$new() |>
#'   estimate(S = 10) |> 
#'   estimate(S = 20) |> 
#'   compute_variance_decompositions(horizon = 4) |>
#'   plot(which_c = "POL")
#' 
#' @export
plot.PosteriorFEVDPANEL = function(
    x,
    which_c,
    cols,
    main,
    xlab,
    mar.multi = c(1, 4.6, 0, 4.6),
    oma.multi = c(6, 0, 5, 0),
    ...
) {
 
  if (is.numeric(which_c)) {
    stopifnot("Argument which_c must be a positive integer indicating one of the countries."
              = length(which_c) == 1 & which_c %% 1 == 0 & which_c > 0 & which_c <= length(x))
  } else if (is.character(which_c)) {
    stopifnot("Argument which_c must be a character string indicating one of the countries."
              = which_c %in% names(x))
  } else {
    stop("Argument which_c must be either a positive integer or a character string.")
  }
  
  plot(
    x[[which_c]],
    cols,
    main,
    xlab,
    mar.multi,
    oma.multi,
    ...
  )
} # END plot.PosteriorFEVDPANEL
