# R6 Class Representing StartingValuesBVARGROUPPRIORPANEL

The class StartingValuesBVARGROUPPRIORPANEL presents starting values for
the Bayesian hierarchical panel VAR model with country grouping

## Public fields

- `group_allocation`:

  a numeric vector with integer numbers denoting group allocations

- `A_c`:

  an `KxNxC` array of starting values for the local parameter
  \\\mathbf{A}\_c\\.

- `Sigma_c`:

  an `NxNxC` array of starting values for the local parameter
  \\\mathbf{\Sigma}\_c\\.

- `A_g`:

  an `KxNxG` array of starting values for the group parameter
  \\\mathbf{A}\_g\\.

- `Sigma_g`:

  an `NxNxG` array of starting values for the group parameter
  \\\mathbf{\Sigma}\_g\\.

- `V`:

  an `KxK` matrix of starting values for the global parameter
  \\\mathbf{V}\\.

- `nu`:

  a positive scalar with starting values for the global parameter
  \\\nu\\.

- `m`:

  a positive scalar with starting values for the global hyper-parameter
  \\m\\.

- `w`:

  a positive scalar with starting values for the global hyper-parameter
  \\w\\.

- `s`:

  a positive scalar with starting values for the global hyper-parameter
  \\s\\.

## Methods

### Public methods

- [`StartingValuesBVARGROUPPRIORPANEL$new()`](#method-StartingValuesBVARGROUPPRIORPANEL-initialize)

- [`StartingValuesBVARGROUPPRIORPANEL$get_starting_values()`](#method-StartingValuesBVARGROUPPRIORPANEL-get_starting_values)

- [`StartingValuesBVARGROUPPRIORPANEL$set_starting_values()`](#method-StartingValuesBVARGROUPPRIORPANEL-set_starting_values)

- [`StartingValuesBVARGROUPPRIORPANEL$clone()`](#method-StartingValuesBVARGROUPPRIORPANEL-clone)

------------------------------------------------------------------------

### `StartingValuesBVARGROUPPRIORPANEL$new()`

Create new starting values StartingValuesBVARGROUPPRIORPANEL

#### Usage

    StartingValuesBVARGROUPPRIORPANEL$new(
      group_allocation = 1:C,
      C,
      G = C,
      N,
      p,
      d = 0
    )

#### Arguments

- `group_allocation`:

  a numeric vector with integer numbers denoting group allocations

- `C`:

  a positive integer - the number of countries in the data.

- `G`:

  a positive integer specifying the number of country groups.

- `N`:

  a positive integer - the number of dependent variables in the model.

- `p`:

  a positive integer - the autoregressive lag order of the SVAR model.

- `d`:

  a positive integer - the number of `exogenous` variables in the model.

#### Returns

Starting values StartingValuesBVARGROUPPRIORPANEL

#### Examples

    # starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
    sv = specify_starting_values_bvarGroupPriorPANEL$new(C = 2, N = 3, p = 1)

------------------------------------------------------------------------

### `StartingValuesBVARGROUPPRIORPANEL$get_starting_values()`

Returns the elements of the starting values
StartingValuesBVARGROUPPRIORPANEL as a `list`.

#### Usage

    StartingValuesBVARGROUPPRIORPANEL$get_starting_values()

#### Examples

    # starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
    sv = specify_starting_values_bvarGroupPriorPANEL$new(rep(1,2), C = 2, N = 3, p = 1)
    sv$get_starting_values()   # show starting values as list

------------------------------------------------------------------------

### `StartingValuesBVARGROUPPRIORPANEL$set_starting_values()`

Returns the elements of the starting values
StartingValuesBVARGROUPPRIORPANEL as a `list`.

#### Usage

    StartingValuesBVARGROUPPRIORPANEL$set_starting_values(last_draw)

#### Arguments

- `last_draw`:

  a list containing the same elements as object
  StartingValuesBVARGROUPPRIORPANEL

#### Returns

An object of class StartingValuesBVARGROUPPRIORPANEL including the last
draw of the current MCMC as the starting value to be passed to the
continuation of the MCMC estimation.

#### Examples

    sv = specify_starting_values_bvarGroupPriorPANEL$new(rep(1,2), C = 2, N = 3, p = 1)

    # Modify the starting values by:
    sv_list = sv$get_starting_values()   # getting them as list
    sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    sv$set_starting_values(sv_list)      # providing to the class object

------------------------------------------------------------------------

### `StartingValuesBVARGROUPPRIORPANEL$clone()`

The objects of this class are cloneable with this method.

#### Usage

    StartingValuesBVARGROUPPRIORPANEL$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# starting values for a Bayesian Panel VAR
sv = specify_starting_values_bvarGroupPriorPANEL$new(rep(1,2), C = 2, G = 1, N = 3, p = 1)


## ------------------------------------------------
## Method `StartingValuesBVARGROUPPRIORPANEL$new()`
## ------------------------------------------------

# starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
sv = specify_starting_values_bvarGroupPriorPANEL$new(C = 2, N = 3, p = 1)


## ------------------------------------------------
## Method `StartingValuesBVARGROUPPRIORPANEL$get_starting_values()`
## ------------------------------------------------

# starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
sv = specify_starting_values_bvarGroupPriorPANEL$new(rep(1,2), C = 2, N = 3, p = 1)
sv$get_starting_values()   # show starting values as list
#> $group_allocation
#> [1] 1 1
#> 
#> $A_c
#> , , 1
#> 
#>               [,1]          [,2]          [,3]
#> [1,]  0.0010854390 -0.0001952577 -1.189237e-03
#> [2,] -0.0011230851 -0.0006162470 -1.620542e-03
#> [3,]  0.0004833151 -0.0016868432 -4.490966e-04
#> [4,]  0.0006956435 -0.0025254535  2.626057e-05
#> 
#> , , 2
#> 
#>               [,1]          [,2]          [,3]
#> [1,] -0.0002268467 -0.0003871743  0.0020317907
#> [2,]  0.0006195275 -0.0009275109 -0.0002915996
#> [3,]  0.0007735930 -0.0026059300 -0.0000645518
#> [4,] -0.0007990488  0.0005208404 -0.0004474033
#> 
#> 
#> $Sigma_c
#> , , 1
#> 
#>           [,1]       [,2]       [,3]
#> [1,] 5.3658937  1.4905428  0.1601149
#> [2,] 1.4905428  0.4961124 -0.2118962
#> [3,] 0.1601149 -0.2118962  2.0230488
#> 
#> , , 2
#> 
#>           [,1]       [,2]       [,3]
#> [1,] 8.1158809  0.1763141  4.6298083
#> [2,] 0.1763141  3.3207061 -0.1427474
#> [3,] 4.6298083 -0.1427474  3.9067713
#> 
#> 
#> $A_g
#> , , 1
#> 
#>               [,1]         [,2]          [,3]
#> [1,] -0.0007856546 8.447181e-04 -0.0001079591
#> [2,] -0.0009051365 4.785831e-04 -0.0006402019
#> [3,]  0.0004795506 7.835054e-04  0.0005378910
#> [4,]  0.0016374050 3.705442e-05 -0.0013234740
#> 
#> , , 2
#> 
#>               [,1]         [,2]          [,3]
#> [1,] -1.842816e-03  0.000386922 -8.296462e-04
#> [2,] -1.404784e-04  0.001538865 -7.219126e-04
#> [3,] -1.357421e-03 -0.000471580 -4.668189e-04
#> [4,] -9.566352e-05  0.001099753  2.357448e-06
#> 
#> 
#> $Sigma_g
#> , , 1
#> 
#>            [,1]       [,2]       [,3]
#> [1,]  3.5512817 -0.7573283 -1.4083869
#> [2,] -0.7573283  0.7362663 -0.1521234
#> [3,] -1.4083869 -0.1521234  1.5235443
#> 
#> , , 2
#> 
#>           [,1]     [,2]     [,3]
#> [1,] 0.8544013 1.437194 1.096880
#> [2,] 1.4371943 5.944971 3.330098
#> [3,] 1.0968797 3.330098 4.911840
#> 
#> 
#> $V
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 15.6255709 -0.2971182  0.9258805 -0.8899933
#> [2,] -0.2971182  2.0571223 -1.8757152  0.6671329
#> [3,]  0.9258805 -1.8757152  2.2618081 -0.1161461
#> [4,] -0.8899933  0.6671329 -0.1161461  4.1104864
#> 
#> $nu
#> [1] 4.1
#> 
#> $m
#> [1] -0.001187267
#> 
#> $w
#> [1] 0.199078
#> 
#> $s
#> [1] 0.4179427
#> 


## ------------------------------------------------
## Method `StartingValuesBVARGROUPPRIORPANEL$set_starting_values()`
## ------------------------------------------------

sv = specify_starting_values_bvarGroupPriorPANEL$new(rep(1,2), C = 2, N = 3, p = 1)

# Modify the starting values by:
sv_list = sv$get_starting_values()   # getting them as list
sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
sv$set_starting_values(sv_list)      # providing to the class object
```
