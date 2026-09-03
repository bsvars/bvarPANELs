# Provides posterior summary of country-specific Forecasts

Provides posterior summary of the forecasts including their mean,
standard deviations, as well as 5 and 95 percentiles.

## Usage

``` r
# S3 method for class 'ForecastsPANEL'
summary(object, which_c, ...)
```

## Arguments

- object:

  an object of class `ForecastsPANEL` obtained using the
  [`forecast()`](https://generics.r-lib.org/reference/forecast.html)
  function containing draws the predictive density.

- which_c:

  a positive integer or a character string specifying the country for
  which the forecast should be plotted.

- ...:

  additional arguments affecting the summary produced.

## Value

A list reporting the posterior mean, standard deviations, as well as 5
and 95 percentiles of the forecasts for each of the variables and
forecast horizons.

## See also

[`forecast.PosteriorBVARPANEL`](http://bsvars.org/bpvars/dev/reference/forecast.PosteriorBVARPANEL.md),
[`plot`](https://rdrr.io/r/graphics/plot.default.html)

## Author

Tomasz Woźniak <wozniak.tom@pm.me>

## Examples

``` r
# specify the model
specification = specify_bvarPANEL$new(
      ilo_dynamic_panel[1:5], 
      exogenous = ilo_exogenous_variables[1:5])
burn_in       = estimate(specification, 5)             # run the burn-in
#> **************************************************|
#> bpvars: Forecasting with Bayesian Panel VARs      |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
posterior     = estimate(burn_in, 5)                   # estimate the model
#> **************************************************|
#> bpvars: Forecasting with Bayesian Panel VARs      |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|

# forecast 3 years ahead
predictive    = forecast(
      posterior, 
      3, 
      exogenous_forecast = ilo_exogenous_forecasts[1:5])
summary(predictive, which_c = "ARG")
#> $`variable 1`
#>       mean         sd 5% quantile 95% quantile
#> 1 27.09490 0.07532608    27.01597     27.17225
#> 2 27.05525 0.13824947    26.89764     27.21759
#> 3 27.07207 0.14301723    26.96418     27.26212
#> 
#> $`variable 2`
#>       mean       sd 5% quantile 95% quantile
#> 1 7.661769 2.358624    4.452608     9.371222
#> 2 8.966183 3.164470    4.601118    10.985469
#> 3 8.727773 3.729852    3.817163    11.901094
#> 
#> $`variable 3`
#>       mean       sd 5% quantile 95% quantile
#> 1 56.60508 1.723899    54.96657     58.85585
#> 2 55.48046 2.271852    53.88060     58.50043
#> 3 55.28687 2.555210    53.05066     58.68573
#> 
#> $`variable 4`
#>       mean        sd 5% quantile 95% quantile
#> 1 61.32586 0.5545387    60.61799     61.86675
#> 2 60.95622 0.5845405    60.50316     61.74591
#> 3 60.56579 0.5043224    60.13459     61.19236
#> 
```
