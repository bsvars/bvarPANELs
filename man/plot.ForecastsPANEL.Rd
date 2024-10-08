% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.ForecastsPANEL}
\alias{plot.ForecastsPANEL}
\title{Plots fitted values of dependent variables}
\usage{
\method{plot}{ForecastsPANEL}(
  x,
  which_c,
  probability = 0.9,
  data_in_plot = 1,
  col = "#1614B1",
  main,
  xlab,
  mar.multi = c(1, 4.6, 0, 2.1),
  oma.multi = c(6, 0, 5, 0),
  ...
)
}
\arguments{
\item{x}{an object of class \code{ForecastsPANEL} obtained using the
\code{forecast()} function containing posterior draws of 
fitted values of dependent variables.}

\item{which_c}{a positive integer or a character string specifying the country 
for which the forecast should be plotted.}

\item{probability}{a parameter determining the interval to be plotted. The 
interval stretches from the \code{0.5 * (1 - probability)} to 
\code{1 - 0.5 * (1 - probability)} percentile of the posterior distribution.}

\item{data_in_plot}{a fraction value in the range (0, 1) determining how many
of the last observations in the data should be plotted with the forecasts.}

\item{col}{a colour of the plot line and the ribbon}

\item{main}{an alternative main title for the plot}

\item{xlab}{an alternative x-axis label for the plot}

\item{mar.multi}{the default \code{mar} argument setting in \code{graphics::par}. Modify with care!}

\item{oma.multi}{the default \code{oma} argument setting in \code{graphics::par}. Modify with care!}

\item{...}{additional arguments affecting the summary produced.}
}
\description{
Plots of fitted values of dependent variables including their 
median and percentiles.
}
\examples{
specification = specify_bvarPANEL$new(ilo_dynamic_panel)  # specify the model
burn_in       = estimate(specification, 10)               # run the burn-in
posterior     = estimate(burn_in, 10)                     # estimate the model

# forecast 6 years ahead
predictive    = forecast(posterior, 6, conditional_forecast = ilo_conditional_forecasts)
plot(predictive, which_c = "POL")                                # plot forecasts

# workflow with the pipe |>
############################################################
set.seed(123)
ilo_dynamic_panel |>
  specify_bvarPANEL$new() |>
  estimate(S = 10) |> 
  estimate(S = 10) |> 
  forecast(horizon = 6, conditional_forecast = ilo_conditional_forecasts) |>
  plot(which_c = 135)

}
\seealso{
\code{\link{forecast.PosteriorBVARPANEL}}
}
\author{
Tomasz Woźniak \email{wozniak.tom@pm.me}
}
