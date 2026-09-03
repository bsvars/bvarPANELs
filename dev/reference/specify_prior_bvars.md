# R6 Class Representing PriorBVARs

The class PriorBVARs presents a prior specification for the Bayesian VAR
model for each country.

## Public fields

- `M`:

  an `KxN` matrix, the mean of the MNIW prior distribution for the
  autoregressive matrices \\\mathbf{A}\_c\\

- `W`:

  a `KxK` column-specific covariance matrix of the MNIW prior
  distribution for the autoregressive matrices \\\mathbf{A}\_c\\

- `S_inv`:

  an `NxN` row-specific precision matrix of the MNIW prior distribution
  for the covariance matrices \\\mathbf{\Sigma}\_c\\

- `lambda`:

  a positive shape of the exponential prior distribution for the shape
  parameter \\\nu\\.

- `mu_m`:

  a scalar mean of the normal prior distribution for the average
  persistence parameter \\m\\.

- `sigma2_m`:

  a positive scalar variance of the normal prior distribution for the
  average persistence parameter \\m\\.

- `s_w`:

  a positive scalar scale of the inverse-gamma 2 prior distribution for
  parameter \\w\\.

- `nu_w`:

  a positive scalar shape of the inverse-gamma 2 prior distribution for
  parameter \\w\\.

- `s_s`:

  a positive scalar scale parameter of the gamma prior distribution for
  parameter \\s\\.

- `a_s`:

  a positive scalar shape parameter of the gamma prior distribution for
  parameter \\s\\.

## Methods

### Public methods

- [`PriorBVARs$new()`](#method-PriorBVARs-initialize)

- [`PriorBVARs$get_prior()`](#method-PriorBVARs-get_prior)

- [`PriorBVARs$clone()`](#method-PriorBVARs-clone)

------------------------------------------------------------------------

### `PriorBVARs$new()`

Create a new prior specification PriorBVARs.

#### Usage

    PriorBVARs$new(C, N, p, d = 0, stationary = rep(FALSE, N))

#### Arguments

- `C`:

  a positive integer - the number of countries in the data.

- `N`:

  a positive integer - the number of dependent variables in the model.

- `p`:

  a positive integer - the autoregressive lag order of the SVAR model.

- `d`:

  a positive integer - the number of `exogenous` variables in the model.

- `stationary`:

  an `N` logical vector - its element set to `FALSE` sets the prior mean
  for the autoregressive parameters of the `N`th equation to the random
  walk process, otherwise to white noise.

#### Returns

A new prior specification PriorBVARs.

#### Examples

    # a prior for 2-country, 3-variable example with one lag and stationary data
    prior = specify_prior_bvars$new(C = 2, N = 3, p = 1)
    prior$M

------------------------------------------------------------------------

### `PriorBVARs$get_prior()`

Returns the elements of the prior specification PriorBVARs as a `list`.

#### Usage

    PriorBVARs$get_prior()

#### Examples

    # a prior for 2-country, 3-variable example with four lags
    prior = specify_prior_bvars$new(C = 2, N = 3, p = 4)
    prior$get_prior() # show the prior as list

------------------------------------------------------------------------

### `PriorBVARs$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PriorBVARs$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
prior = specify_prior_bvars$new(C = 2, N = 3, p = 1)
prior$M
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> [4,]    0    0    0


## ------------------------------------------------
## Method `PriorBVARs$new()`
## ------------------------------------------------

# a prior for 2-country, 3-variable example with one lag and stationary data
prior = specify_prior_bvars$new(C = 2, N = 3, p = 1)
prior$M
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> [4,]    0    0    0


## ------------------------------------------------
## Method `PriorBVARs$get_prior()`
## ------------------------------------------------

# a prior for 2-country, 3-variable example with four lags
prior = specify_prior_bvars$new(C = 2, N = 3, p = 4)
prior$get_prior() # show the prior as list
#> $M
#>       [,1] [,2] [,3]
#>  [1,]    1    0    0
#>  [2,]    0    1    0
#>  [3,]    0    0    1
#>  [4,]    0    0    0
#>  [5,]    0    0    0
#>  [6,]    0    0    0
#>  [7,]    0    0    0
#>  [8,]    0    0    0
#>  [9,]    0    0    0
#> [10,]    0    0    0
#> [11,]    0    0    0
#> [12,]    0    0    0
#> [13,]    0    0    0
#> 
#> $W
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#>  [1,]    1    0    0    0    0    0    0    0    0     0     0     0     0
#>  [2,]    0    1    0    0    0    0    0    0    0     0     0     0     0
#>  [3,]    0    0    1    0    0    0    0    0    0     0     0     0     0
#>  [4,]    0    0    0    4    0    0    0    0    0     0     0     0     0
#>  [5,]    0    0    0    0    4    0    0    0    0     0     0     0     0
#>  [6,]    0    0    0    0    0    4    0    0    0     0     0     0     0
#>  [7,]    0    0    0    0    0    0    9    0    0     0     0     0     0
#>  [8,]    0    0    0    0    0    0    0    9    0     0     0     0     0
#>  [9,]    0    0    0    0    0    0    0    0    9     0     0     0     0
#> [10,]    0    0    0    0    0    0    0    0    0    16     0     0     0
#> [11,]    0    0    0    0    0    0    0    0    0     0    16     0     0
#> [12,]    0    0    0    0    0    0    0    0    0     0     0    16     0
#> [13,]    0    0    0    0    0    0    0    0    0     0     0     0    10
#> 
#> $S_inv
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $lambda
#> [1] 30
#> 
#> $mu_m
#> [1] 1
#> 
#> $sigma2_m
#> [1] 0.1
#> 
#> $s_w
#> [1] 1
#> 
#> $nu_w
#> [1] 3
#> 
#> $s_s
#> [1] 1
#> 
#> $a_s
#> [1] 1
#> 
```
