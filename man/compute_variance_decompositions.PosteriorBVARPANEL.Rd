% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/compute_variance_decompositions.PosteriorBVARPANEL.R
\name{compute_variance_decompositions.PosteriorBVARPANEL}
\alias{compute_variance_decompositions.PosteriorBVARPANEL}
\title{Computes posterior draws of the forecast error variance decomposition}
\usage{
\method{compute_variance_decompositions}{PosteriorBVARPANEL}(posterior, horizon)
}
\arguments{
\item{posterior}{posterior estimation outcome - an object of class 
\code{PosteriorBVARPANEL} obtained by running the \code{estimate} function.}

\item{horizon}{a positive integer number denoting the forecast horizon for 
the forecast error variance decompositions.}
}
\value{
An object of class \code{PosteriorFEVDPANEL}, that is, a list with 
\code{C} elements containing \code{NxNx(horizon+1)xS} arrays of class 
\code{PosteriorFEVD} with \code{S} draws of country-specific forecast error 
variance decompositions.
}
\description{
For each country, each of the draws from the posterior estimation 
of the model is transformed into a draw from the posterior distribution of the forecast 
error variance decomposition.
}
\examples{
# upload data
data(ilo_dynamic_panel)

# specify the model and set seed
set.seed(123)
specification  = specify_bvarPANEL$new(ilo_dynamic_panel, p = 1)

# run the burn-in
burn_in        = estimate(specification, 10)

# estimate the model
posterior      = estimate(burn_in, 20)

# compute forecast error variance decomposition 4 years ahead
fevd           = compute_variance_decompositions(posterior, horizon = 4)

# workflow with the pipe |>
############################################################
set.seed(123)
ilo_dynamic_panel |>
  specify_bvarPANEL$new(p = 1) |>
  estimate(S = 10) |> 
  estimate(S = 20) |> 
  compute_variance_decompositions(horizon = 4) -> fevd

}
\references{
Lütkepohl, H. (2017). Structural VAR Tools, Chapter 4, In: Structural vector autoregressive analysis. Cambridge University Press.
}
\seealso{
\code{\link{estimate.PosteriorBVARPANEL}}, 
\code{\link{summary.PosteriorFEVDPANEL}}, 
\code{\link{plot.PosteriorFEVDPANEL}}
}
\author{
Tomasz Woźniak \email{wozniak.tom@pm.me}
}
