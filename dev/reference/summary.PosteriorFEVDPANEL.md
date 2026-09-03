# Provides posterior summary of forecast error variance decompositions

Provides posterior means of the forecast error variance decompositions
of each variable at all horizons.

## Usage

``` r
# S3 method for class 'PosteriorFEVDPANEL'
summary(object, which_c, ...)
```

## Arguments

- object:

  an object of class `PosteriorFEVDPANEL` obtained using the
  [`compute_variance_decompositions()`](https://bsvars.org/bsvars/reference/compute_variance_decompositions.html)
  function containing draws from the posterior distribution of the
  forecast error variance decompositions.

- which_c:

  a positive integer or a character string specifying the country for
  which the forecast should be plotted.

- ...:

  additional arguments affecting the summary produced.

## Value

A list reporting the posterior mean of the forecast error variance
decompositions of each variable at all horizons.

## See also

[`compute_variance_decompositions.PosteriorBVARPANEL`](http://bsvars.org/bpvars/dev/reference/compute_variance_decompositions.PosteriorBVARPANEL.md),
[`plot`](https://rdrr.io/r/graphics/plot.default.html)

## Author

Tomasz Woźniak <wozniak.tom@pm.me>

## Examples

``` r
# specify the model and set seed
specification  = specify_bvarPANEL$new(ilo_dynamic_panel[1:5], p = 1)

# run the burn-in
burn_in        = estimate(specification, 5)
#> **************************************************|
#> bpvars: Forecasting with Bayesian Panel VARs      |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|

# estimate the model
posterior      = estimate(burn_in, 5)
#> **************************************************|
#> bpvars: Forecasting with Bayesian Panel VARs      |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|

# compute forecast error variance decomposition 4 years ahead
fevd           = compute_variance_decompositions(posterior, horizon = 4)
summary(fevd, which_c = "ARG")
#> $gdp
#>      shock1   shock2   shock3   shock4
#> 0 100.00000 0.000000  0.00000  0.00000
#> 1  62.60012 7.336113 17.61172 12.45205
#> 2  58.38758 8.146548 19.60396 13.86191
#> 3  56.41115 8.522853 20.53992 14.52607
#> 4  55.15127 8.761092 21.13676 14.95087
#> 
#> $UR
#>     shock1     shock2   shock3   shock4
#> 0 99.17236  0.8276363  0.00000  0.00000
#> 1 26.26367 15.1733689 34.61127 23.95169
#> 2 26.52568 14.0200824 34.54052 24.91372
#> 3 26.08783 13.9633933 34.77841 25.17037
#> 4 26.35455 13.8985996 34.66400 25.08285
#> 
#> $EPR
#>     shock1    shock2   shock3   shock4
#> 0 79.31049  8.466595 12.22291  0.00000
#> 1 18.25950 16.074478 38.43187 27.23415
#> 2 24.61275 14.666349 35.50098 25.21992
#> 3 27.41521 14.035865 34.20968 24.33924
#> 4 28.55232 13.780424 33.68156 23.98569
#> 
#> $LFPR
#>      shock1   shock2   shock3   shock4
#> 0  9.236772 15.98483 42.01035 32.76804
#> 1 12.959706 16.30191 40.77315 29.96524
#> 2 19.744691 15.11950 37.65608 27.47973
#> 3 22.835913 14.59352 36.23612 26.33445
#> 4 24.308039 14.34704 35.56105 25.78387
#> 
```
