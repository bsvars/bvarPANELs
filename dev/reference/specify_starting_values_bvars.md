# R6 Class Representing StartingValuesBVARs

The class StartingValuesBVARs presents starting values for the Bayesian
hierarchical panel VAR model.

## Public fields

- `A_c`:

  an `KxNxC` array of starting values for the local parameter
  \\\mathbf{A}\_c\\.

- `Sigma_c`:

  an `NxNxC` array of starting values for the local parameter
  \\\mathbf{\Sigma}\_c\\.

- `nu`:

  a `C`-vector of positive starting values for the parameter \\\nu\\.

- `m`:

  a `C`-vector of starting values for the parameter \\m\\.

- `w`:

  a `C`-vector of positive starting values for the parameter \\w\\.

- `s`:

  a `C`-vector of positive starting values for the parameter \\s\\.

## Methods

### Public methods

- [`StartingValuesBVARs$new()`](#method-StartingValuesBVARs-initialize)

- [`StartingValuesBVARs$get_starting_values()`](#method-StartingValuesBVARs-get_starting_values)

- [`StartingValuesBVARs$set_starting_values()`](#method-StartingValuesBVARs-set_starting_values)

- [`StartingValuesBVARs$clone()`](#method-StartingValuesBVARs-clone)

------------------------------------------------------------------------

### `StartingValuesBVARs$new()`

Create new starting values StartingValuesBVARs

#### Usage

    StartingValuesBVARs$new(C, N, p, d = 0)

#### Arguments

- `C`:

  a positive integer - the number of countries in the data.

- `N`:

  a positive integer - the number of dependent variables in the model.

- `p`:

  a positive integer - the autoregressive lag order of the SVAR model.

- `d`:

  a positive integer - the number of `exogenous` variables in the model.

#### Returns

Starting values StartingValuesBVARs

#### Examples

    # starting values for Bayesian VARs 2-country model with 4 lags for a 3-variable system.
    sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 4)

------------------------------------------------------------------------

### `StartingValuesBVARs$get_starting_values()`

Returns the elements of the starting values StartingValuesBVARs as a
`list`.

#### Usage

    StartingValuesBVARs$get_starting_values()

#### Examples

    # starting values for bvars with 1 lag for a 3-variable system
    sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 1)
    sv$get_starting_values()   # show starting values as list

------------------------------------------------------------------------

### `StartingValuesBVARs$set_starting_values()`

Returns the elements of the starting values StartingValuesBVARs as a
`list`.

#### Usage

    StartingValuesBVARs$set_starting_values(last_draw)

#### Arguments

- `last_draw`:

  a list containing the same elements as object StartingValuesBVARs.

#### Returns

An object of class StartingValuesBVARs including the last draw of the
current MCMC as the starting value to be passed to the continuation of
the MCMC estimation.

#### Examples

    sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 1)

    # Modify the starting values by:
    sv_list = sv$get_starting_values()   # getting them as list
    sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    sv$set_starting_values(sv_list)      # providing to the class object

------------------------------------------------------------------------

### `StartingValuesBVARs$clone()`

The objects of this class are cloneable with this method.

#### Usage

    StartingValuesBVARs$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# starting values for a Bayesian Panel VAR
sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 1)


## ------------------------------------------------
## Method `StartingValuesBVARs$new()`
## ------------------------------------------------

# starting values for Bayesian VARs 2-country model with 4 lags for a 3-variable system.
sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 4)


## ------------------------------------------------
## Method `StartingValuesBVARs$get_starting_values()`
## ------------------------------------------------

# starting values for bvars with 1 lag for a 3-variable system
sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 1)
sv$get_starting_values()   # show starting values as list
#> $A_c
#> , , 1
#> 
#>               [,1]          [,2]          [,3]
#> [1,] -0.0011085675 -0.0001978905 -0.0013404320
#> [2,] -0.0004023415 -0.0021867963  0.0003385493
#> [3,]  0.0008131702 -0.0007578514 -0.0015615494
#> [4,]  0.0000595284  0.0013349899  0.0017041468
#> 
#> , , 2
#> 
#>               [,1]          [,2]          [,3]
#> [1,]  5.413216e-05  0.0001977299 -0.0001768038
#> [2,]  3.775696e-04 -0.0001645073 -0.0005625921
#> [3,] -9.796387e-04  0.0002826896  0.0010619837
#> [4,] -1.067624e-03 -0.0003604902  0.0007249878
#> 
#> 
#> $Sigma_c
#> , , 1
#> 
#>           [,1]      [,2]      [,3]
#> [1,]  8.528677  8.324654 -1.665982
#> [2,]  8.324654 10.461182 -1.722317
#> [3,] -1.665982 -1.722317  1.821717
#> 
#> , , 2
#> 
#>           [,1]      [,2]       [,3]
#> [1,] 5.9911612  2.126774  0.8801813
#> [2,] 2.1267736  4.206568 -2.3698911
#> [3,] 0.8801813 -2.369891  2.8129706
#> 
#> 
#> $nu
#> [1] 4.1 4.1
#> 
#> $m
#> [1] -0.0006002325  0.0016053817
#> 
#> $w
#> [1] 0.9798671 3.4109203
#> 
#> $s
#> [1] 0.8229502 0.7573296
#> 


## ------------------------------------------------
## Method `StartingValuesBVARs$set_starting_values()`
## ------------------------------------------------

sv = specify_starting_values_bvars$new(C = 2, N = 3, p = 1)

# Modify the starting values by:
sv_list = sv$get_starting_values()   # getting them as list
sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
sv$set_starting_values(sv_list)      # providing to the class object
```
