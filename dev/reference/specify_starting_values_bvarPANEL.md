# R6 Class Representing StartingValuesBVARPANEL

The class StartingValuesBVARPANEL presents starting values for the
Bayesian hierarchical panel VAR model.

## Public fields

- `A_c`:

  an `KxNxC` array of starting values for the local parameter
  \\\mathbf{A}\_c\\.

- `Sigma_c`:

  an `NxNxC` array of starting values for the local parameter
  \\\mathbf{\Sigma}\_c\\.

- `A`:

  an `KxN` matrix of starting values for the global parameter
  \\\mathbf{A}\\.

- `V`:

  an `KxK` matrix of starting values for the global parameter
  \\\mathbf{V}\\.

- `Sigma`:

  an `NxN` matrix of starting values for the global parameter
  \\\mathbf{\Sigma}\\.

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

- [`StartingValuesBVARPANEL$new()`](#method-StartingValuesBVARPANEL-initialize)

- [`StartingValuesBVARPANEL$get_starting_values()`](#method-StartingValuesBVARPANEL-get_starting_values)

- [`StartingValuesBVARPANEL$set_starting_values()`](#method-StartingValuesBVARPANEL-set_starting_values)

- [`StartingValuesBVARPANEL$clone()`](#method-StartingValuesBVARPANEL-clone)

------------------------------------------------------------------------

### `StartingValuesBVARPANEL$new()`

Create new starting values StartingValuesBVARPANEL

#### Usage

    StartingValuesBVARPANEL$new(C, N, p, d = 0)

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

Starting values StartingValuesBVARPANEL

#### Examples

    # starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
    sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 4)

------------------------------------------------------------------------

### `StartingValuesBVARPANEL$get_starting_values()`

Returns the elements of the starting values StartingValuesBVARPANEL as a
`list`.

#### Usage

    StartingValuesBVARPANEL$get_starting_values()

#### Examples

    # starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
    sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)
    sv$get_starting_values()   # show starting values as list

------------------------------------------------------------------------

### `StartingValuesBVARPANEL$set_starting_values()`

Returns the elements of the starting values StartingValuesBVARPANEL as a
`list`.

#### Usage

    StartingValuesBVARPANEL$set_starting_values(last_draw)

#### Arguments

- `last_draw`:

  a list containing the same elements as object StartingValuesBVARPANEL.

#### Returns

An object of class StartingValuesBVARPANEL including the last draw of
the current MCMC as the starting value to be passed to the continuation
of the MCMC estimation.

#### Examples

    sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)

    # Modify the starting values by:
    sv_list = sv$get_starting_values()   # getting them as list
    sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    sv$set_starting_values(sv_list)      # providing to the class object

------------------------------------------------------------------------

### `StartingValuesBVARPANEL$clone()`

The objects of this class are cloneable with this method.

#### Usage

    StartingValuesBVARPANEL$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# starting values for a Bayesian Panel VAR
sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)


## ------------------------------------------------
## Method `StartingValuesBVARPANEL$new()`
## ------------------------------------------------

# starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 4)


## ------------------------------------------------
## Method `StartingValuesBVARPANEL$get_starting_values()`
## ------------------------------------------------

# starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)
sv$get_starting_values()   # show starting values as list
#> $A_c
#> , , 1
#> 
#>               [,1]         [,2]         [,3]
#> [1,]  0.0001852815 0.0008190492 0.0003883010
#> [2,]  0.0008873240 0.0005189019 0.0004066908
#> [3,]  0.0002059307 0.0003624973 0.0009074711
#> [4,] -0.0022305519 0.0008953130 0.0004045349
#> 
#> , , 2
#> 
#>              [,1]          [,2]          [,3]
#> [1,] 1.549831e-03 -0.0004536336  9.979931e-05
#> [2,] 4.397918e-04  0.0007811692  5.313644e-04
#> [3,] 1.272904e-03  0.0006758488 -9.064482e-04
#> [4,] 2.125325e-05  0.0007196550 -1.020071e-03
#> 
#> 
#> $Sigma_c
#> , , 1
#> 
#>             [,1]       [,2]        [,3]
#> [1,]  0.63697928 -0.1340802  0.02598222
#> [2,] -0.13408016  4.1581714 -0.13626879
#> [3,]  0.02598222 -0.1362688  1.13343815
#> 
#> , , 2
#> 
#>           [,1]        [,2]        [,3]
#> [1,]  6.962050 -1.19908343 -1.35700920
#> [2,] -1.199083  6.58234263 -0.08728492
#> [3,] -1.357009 -0.08728492  1.18854308
#> 
#> 
#> $A
#>               [,1]          [,2]          [,3]
#> [1,]  9.999070e-01  0.0011537442 -0.0007686428
#> [2,]  1.090224e-03  1.0003949745 -0.0003847204
#> [3,] -6.428193e-05 -0.0005097384  1.0009751815
#> [4,] -5.323123e-05 -0.0005282103  0.0001009873
#> 
#> $V
#>             [,1]       [,2]       [,3]        [,4]
#> [1,]  2.59563991  0.2130839  0.7007901 -0.06621208
#> [2,]  0.21308392  4.4276015 -0.2114750 -3.05028808
#> [3,]  0.70079014 -0.2114750  2.4657567 -0.21524988
#> [4,] -0.06621208 -3.0502881 -0.2152499  2.37413718
#> 
#> $Sigma
#>           [,1]     [,2]      [,3]
#> [1,]  4.773482 2.410161 -2.829625
#> [2,]  2.410161 8.139780  1.150940
#> [3,] -2.829625 1.150940  6.214514
#> 
#> $nu
#> [1] 4.1
#> 
#> $m
#> [1] 0.001868389
#> 
#> $w
#> [1] 2.103943
#> 
#> $s
#> [1] 1.403311
#> 


## ------------------------------------------------
## Method `StartingValuesBVARPANEL$set_starting_values()`
## ------------------------------------------------

sv = specify_starting_values_bvarPANEL$new(C = 2, N = 3, p = 1)

# Modify the starting values by:
sv_list = sv$get_starting_values()   # getting them as list
sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
sv$set_starting_values(sv_list)      # providing to the class object
```
