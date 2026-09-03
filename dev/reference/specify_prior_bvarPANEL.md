# R6 Class Representing PriorBVARPANEL

The class PriorBVARPANEL presents a prior specification for the Bayesian
hierarchical panel VAR model.

## Public fields

- `M`:

  an `KxN` matrix, the mean of the second-level MNIW prior distribution
  for the global parameter matrices \\\mathbf{A}\\ and \\\mathbf{V}\\

- `W`:

  a `KxK` column-specific covariance matrix of the second-level MNIW
  prior distribution for the global parameter matrices \\\mathbf{A}\\
  and \\\mathbf{V}\\

- `S_inv`:

  an `NxN` row-specific precision matrix of the second-level MNIW prior
  distribution for the global parameter matrices \\\mathbf{A}\\ and
  \\\mathbf{V}\\

- `S_Sigma_inv`:

  an `NxN` precision matrix of the second-level Wishart prior
  distribution for the global parameter matrix \\\mathbf{\Sigma}\\.

- `eta`:

  a positive shape parameter of the second-level MNIW prior distribution
  for the global parameter matrices \\\mathbf{A}\\ and \\\mathbf{V}\\

- `mu_Sigma`:

  a positive shape parameter of the second-level Wishart prior
  distribution for the global parameter matrix \\\mathbf{\Sigma}\\.

- `lambda`:

  a positive shape of the second-level exp prior distribution for the
  shape parameter \\\nu\\.

- `mu_m`:

  a scalar mean of the third-level normal prior distribution for the
  global average persistence parameter \\m\\.

- `sigma2_m`:

  a positive scalar variance of the third-level normal prior
  distribution for the global average persistence parameter \\m\\.

- `s_w`:

  a positive scalar scale of the third-level gamma prior distribution
  for parameter \\w\\.

- `a_w`:

  a positive scalar shape of the third-level gamma prior distribution
  for parameter \\w\\.

- `s_s`:

  a positive scalar scale parameter of the third-level inverted-gamma 2
  prior distribution for parameter \\s\\.

- `nu_s`:

  a positive scalar shape parameter of the third-level inverted-gamma 2
  prior distribution for parameter \\s\\.

## Methods

### Public methods

- [`PriorBVARPANEL$new()`](#method-PriorBVARPANEL-initialize)

- [`PriorBVARPANEL$get_prior()`](#method-PriorBVARPANEL-get_prior)

- [`PriorBVARPANEL$clone()`](#method-PriorBVARPANEL-clone)

------------------------------------------------------------------------

### `PriorBVARPANEL$new()`

Create a new prior specification PriorBVARPANEL.

#### Usage

    PriorBVARPANEL$new(C, N, p, d = 0, stationary = rep(FALSE, N))

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

A new prior specification PriorBVARPANEL.

#### Examples

    # a prior for 2-country, 3-variable example with one lag and stationary data
    prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 1)
    prior$M

------------------------------------------------------------------------

### `PriorBVARPANEL$get_prior()`

Returns the elements of the prior specification PriorBSVAR as a `list`.

#### Usage

    PriorBVARPANEL$get_prior()

#### Examples

    # a prior for 2-coutnry, 3-variable example with four lags
    prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 4)
    prior$get_prior() # show the prior as list

------------------------------------------------------------------------

### `PriorBVARPANEL$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PriorBVARPANEL$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 1)
prior$M
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> [4,]    0    0    0


## ------------------------------------------------
## Method `PriorBVARPANEL$new()`
## ------------------------------------------------

# a prior for 2-country, 3-variable example with one lag and stationary data
prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 1)
prior$M
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> [4,]    0    0    0


## ------------------------------------------------
## Method `PriorBVARPANEL$get_prior()`
## ------------------------------------------------

# a prior for 2-coutnry, 3-variable example with four lags
prior = specify_prior_bvarPANEL$new(C = 2, N = 3, p = 4)
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
#> $S_Sigma_inv
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
#> 
#> $eta
#> [1] 4
#> 
#> $mu_Sigma
#> [1] 4
#> 
#> $lambda
#> [1] 72
#> 
#> $mu_m
#> [1] 1
#> 
#> $sigma2_m
#> [1] 1
#> 
#> $s_w
#> [1] 1
#> 
#> $a_w
#> [1] 1
#> 
#> $s_s
#> [1] 1
#> 
#> $nu_s
#> [1] 3
#> 
```
