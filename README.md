
# bpvars <a href="https://bsvars.org/bpvars/"><img src="man/figures/logo.png" align="right" height="139" alt="bpvars website" /></a>

An **R** package for Forecasting with Bayesian Panel Vector
Autoregressions

<!-- badges: start -->

[![R-CMD-check](https://github.com/bsvars/bpvars/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bpvars/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Provides Bayesian estimation and forecasting of dynamic panel data using
Bayesian Panel Vector Autoregressions with hierarchical prior
distributions following the specification by [Sanchez-Martinez & Woźniak
(2026)](https://doi.org/10.48550/arXiv.2606.14143). The models include
country-specific Vector Autoregressions (VARs) that share a global prior
distribution that extend the model by [Jarociński
(2010)](https://doi.org/10.1002/jae.1082). Under this prior expected
value, each country’s system follows a global VAR with country-invariant
parameters. Further flexibility is provided by the hierarchical prior
structure that retains the Minnesota prior interpretation for the global
VAR and features estimated prior covariance matrices, shrinkage, and
persistence levels. Bayesian forecasting is developed for models
including exogenous variables, allowing conditional forecasts given the
future trajectories of some variables and restricted forecasts assuring
that rates are forecasted to stay positive and less than 100. The
package implements the model specification, estimation, and forecasting
routines, facilitating coherent workflows and reproducibility. Beautiful
plots, informative summary functions, and extensive documentation
complement all this. Extraordinary computational speed is achieved
thanks to employing frontier econometric and numerical techniques and
algorithms written in **C++**. The **bpvars** package is aligned
regarding objects, workflows, and code structure with the **R** packages
**bsvars** by [Woźniak
(2024)](https://cran.r-project.org/package=bsvars), **bsvarSIGNs** by
[Wang & Woźniak (2024)](https://cran.r-project.org/package=bsvarSIGNs),
and **bvars** by [Liu, Ramirez Hassan, & Woźniak
(2026)](https://doi.org/10.32614/CRAN.package.bvars), and they
constitute an integrated toolset.

<a href="https://bsvars.org">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/solid/house.svg" width="40" height="40"/>
</a> <a href="mailto:contact@bsvars.org">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/solid/envelope.svg" width="40" height="40"/>
</a> <a href="https://github.com/bsvars/bpvars">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/brands/github.svg" width="40" height="40"/>
</a> <a href="https://bsky.app/profile/bsvars.org">
<img src="https://upload.wikimedia.org/wikipedia/commons/7/7a/Bluesky_Logo.svg" width="40" height="40"/>
</a> <a href="https://fosstodon.org/@bsvars">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/brands/mastodon.svg" width="40" height="40"/>
</a>

<a href="https://bsvars.org/"><img src="https://raw.githubusercontent.com/bsvars/hex/refs/heads/main/bsvars.org/bsvars.org.png" width="120" alt="bsvars.org website" /></a>
<a href="https://bsvars.org/bsvars/"><img src="https://raw.githubusercontent.com/bsvars/hex/refs/heads/main/bsvars/bsvars.png" width="120" alt="bsvars website" /></a>
<a href="https://bsvars.org/bsvarSIGNs/"><img src="https://raw.githubusercontent.com/bsvars/hex/refs/heads/main/bsvarSIGNs/bsvarSIGNs.png" width="120" alt="bsvarSIGNs website" /></a>
<a href="https://bsvars.org/bpvars/"><img 
src="https://raw.githubusercontent.com/bsvars/hex/refs/heads/main/bpvars/bpvars.png" width="120" alt="bpvars website" /></a>
<a href="https://bsvars.org/bvars/"><img 
src="https://raw.githubusercontent.com/bsvars/hex/refs/heads/main/bvars/bvars.png" width="120" alt="bvars website" /></a>
<a href="https://bsvars.org/StealLikeBayes/"><img src="https://raw.githubusercontent.com/bsvars/hex/refs/heads/main/StealLikeBayes/StealLikeBayes.png" width="120" alt="StealLikeBayes website" /></a>

## Features

#### Forecasting with Bayesian Panel Vector Autoregressions

- The model in the **bpvars** package features a country-specific Vector
  Autoregressive equation for the country-specific:
  - dependent variables `Yc`,
  - lagged dependent variables `Xc`,
  - error terms `Ec`,
  - autoregressive parameters `Ac`, and
  - error term covariance matrix `Sc`.

  <!-- -->

           Yc = Xc Ac + Ec         (VAR equation)
      Ec | Xc ~ MN(0, I, Sc)       (error term normality)
- The error terms feature a zero-mean matrix-variate normal distribution
  with row-specific covariance matrix `Sc` and column-specific
  covariance equal to the identity matrix of order `T`, denoted by `I`.
- The Hierarchical Panel VAR model features a sophisticated hierarchical
  prior structure that grants the model flexibility, interpretability,
  and improved forecasting performance.
- The country-specific parameters follow a prior distribution that, at
  its mean value, represents a global VAR model with a global
  autoregressive parameter matrix `A` and a global error term covariance
  matrix `S`:

<!-- -->

           Yc = Xc A + Ec          (VAR equation)
      Ec | Xc ~ MN(0, I, S)        (error term normality)

- See more details in the
  [vignette](https://doi.org/10.48550/arXiv.2606.14143) and this
  [presentation](https://bsvars.org/2026-03-bpvars-ilo/)

#### Simple workflows

- Specify the models using the `specify_bvarPANEL$new()` function
- Estimate the models using the `estimate()` method
- Predict the future using the `forecast()` method
- Compute forecast error variance decompositions using function
  `compute_variance_decompositions()`
- Use `plot()` and `summary()` methods to gain the insights into the
  core of the empirical problem.

#### Fast and efficient computations

- Extraordinary computational speed is obtained by combining
  - the application of frontier econometric and numerical techniques,
    and
  - the algorithms written in **cpp**
- It combines the best of two worlds: the ease of data analysis with
  **R** and fast **cpp** algorithms
- The algorithms used here are very fast. But still, Bayesian estimation
  might take a little time. Look at our beautiful **progress bar** in
  the meantime:

<!-- -->

    **************************************************|
    bpvars: Forecasting with Bayesian Hierarchical    |
                Panel Vector Autoregressions          |
    **************************************************|
     Progress of the MCMC simulation for 10000 draws
        Every draw is saved via MCMC thinning
     Press Esc to interrupt the computations
    **************************************************|
    0%   10   20   30   40   50   60   70   80   90   100%
    [----|----|----|----|----|----|----|----|----|----|
    *************************************

#### The hexagonal logo

This beautiful logo can be reproduced in R using [this
file](https://github.com/bsvars/hex/blob/main/bpvars/bpvars.R).

<p>

</p>

<a href="https://bsvars.org/bpvars/"><img src="man/figures/logo.png" height="400" alt="bpvars website" /></a>
<p>

</p>

## Resources

> [\[CRAN\]](https://cran.r-project.org/package=bpvars)
> [\[website\]](https://bsvars.org/bpvars/)
> [\[vignette\]](https://doi.org/10.48550/arXiv.2606.14143)
> [\[manual\]](https://cran.r-project.org/package=bpvars/refman/bpvars.html)
> [\[repo\]](https://github.com/bsvars/bpvars)
> [\[resources\]](https://bsvars.org/bpvars/#resources)
> [\[bsvars.org\]](https://bsvars.org/)
>
> **Youtube Recordings**\
> \[[International Labour Organization
> 2026-03](https://www.ilo.org/meetings-and-events/forecasting-labour-market-outcomes-using-bayesian-hierarchical-panel-vars)
> [youtube recording](https://www.youtube.com/watch?v=ef3eXbqNbr8)\]\
>
> **Presentations**\
> \[[Workshops for
> Ukraine](https://sites.google.com/view/dariia-mykhailyshyna/main/r-workshops-for-ukraine)
> [2026-09 featuring **bpvars**
> 2.0](http://bsvars.org/2026-09-bpvars-w4UKR/)\]\
> \[[International Labour Organization](https://www.ilo.org/) [2026-03
> featuring **bpvars** 1.0](https://bsvars.org/2026-03-bpvars-ilo/)\]\
> \[[International Labour Organization](https://www.ilo.org/) [2025-03
> featuring **bvarPANELs**
> 0.2](https://bsvars.org/2025-03-bvarPANELs-ilo/)\]\

## Start your Bayesian analysis of data

The beginnings are as easy as ABC:

``` r
library(bpvars)                                     # load the package

spec = specify_bvarPANEL$new(                           # specify the model
  ilo_dynamic_panel                                     # data
)

burn = estimate(spec, S = 10000)                        # run the burn-in
post = estimate(burn, S = 10000)                        # estimate the model

fore = forecast(                                        # forecast the model 
  post,                                                 # estimation output
  horizon = 6                                           # forecast horizon
)

plot(fore, "COL", main = "Forecasts for Colombia")
summary(fore, "COL")$variable2

post |>                                                 # estimation output
  compute_variance_decompositions(horizon = 6) |>       # compute variance decompositions
  plot(which_c = "COL")                                 # plot variance decompositions
```

The **bpvars** package supports a simplified workflow using the `|>`
pipe:

``` r
ilo_dynamic_panel |>                                    # data
  specify_bvarPANEL$new() |>                            # specify the model
  estimate(S = 10000) |>                                # run the burn-in
  estimate(S = 10000) -> post                           # estimate the model

post |> forecast(                                       # forecast the model 
  horizon = 6                                           # forecast horizon
) |> plot("COL", main = "Forecasts for Colombia")
```

Now, you’re ready to analyse your model and forecasts!

## Installation

#### The first time you install the package

You must have a **cpp** compiler. Follow the instructions from [Section
1.3. by Eddelbuettel & François
(2023)](https://cran.r-project.org/package=Rcpp/vignettes/Rcpp-FAQ.pdf).
In short, for **Windows:** install
[RTools](https://CRAN.R-project.org/bin/windows/Rtools/), for **macOS:**
install [Xcode Command Line
Tools](https://www.freecodecamp.org/news/install-xcode-command-line-tools/),
and for **Linux:** install the standard development packages.

#### Once that’s done:

Just open your **R** and install the package by running:

    instal.packages("bpvars")

To install the package from its developer’s repository just type:

    devtools::install_github("bsvars/bpvars")

## Development

The package is under intensive development. Your help is most welcome!
Please, have a look at our
[issues](https://github.com/bsvars/bsvars/issues) to learn what we’re
working on. Thank you!

## About the author

**Tomasz** is a Bayesian econometrician and a Senior Lecturer at the
University of Melbourne. He develops methodology for empirical
macroeconomic analyses and programs in **R** and **C++** using **Rcpp**.

<a href="mailto:twozniak@unimelb.edu.au">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/solid/envelope.svg" width="40" height="40"/>
</a> <a href="https://github.com/donotdespair">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/brands/github.svg" width="40" height="40"/>
</a> <a href="https://orcid.org/0000-0003-2212-2378">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/brands/orcid.svg" width="40" height="40"/>
</a> <a href="https://www.linkedin.com/in/tomaszwwozniak">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/brands/linkedin.svg" width="40" height="40"/>
</a>
<a href="https://scholar.google.com/citations?user=2uWpFrYAAAAJ&hl">
<img src="https://raw.githubusercontent.com/jpswalsh/academicons/refs/heads/master/svg/google-scholar-square.svg" width="40" height="40"/>
</a> <a href="https://arxiv.org/a/wozniak_t_1">
<img src="https://raw.githubusercontent.com/jpswalsh/academicons/refs/heads/master/svg/arxiv-square.svg" width="40" height="40"/>
</a> <a href="https://www.researchgate.net/profile/Tomasz-Wozniak-2">
<img src="https://raw.githubusercontent.com/jpswalsh/academicons/refs/heads/master/svg/researchgate-square.svg" width="40" height="40"/>
</a> <a href="https://fosstodon.org/@tomaszwozniak">
<img src="https://raw.githubusercontent.com/FortAwesome/Font-Awesome/6.x/svgs/brands/mastodon.svg" width="40" height="40"/>
</a> <a href="https://bsky.app/profile/tomaszwozniak.bsky.social">
<img src="https://upload.wikimedia.org/wikipedia/commons/7/7a/Bluesky_Logo.svg" width="40" height="40"/>
</a>
