# R6 Class Representing StartingValuesBVARGROUPPANEL

The class StartingValuesBVARGROUPPANEL presents starting values for the
Bayesian hierarchical panel VAR model with country grouping

## Super class

`StartingValuesBVARPANEL` -\> `StartingValuesBVARGROUPPANEL`

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

- [`StartingValuesBVARGROUPPANEL$new()`](#method-StartingValuesBVARGROUPPANEL-initialize)

- [`StartingValuesBVARGROUPPANEL$get_starting_values()`](#method-StartingValuesBVARGROUPPANEL-get_starting_values)

- [`StartingValuesBVARGROUPPANEL$set_starting_values()`](#method-StartingValuesBVARGROUPPANEL-set_starting_values)

- [`StartingValuesBVARGROUPPANEL$clone()`](#method-StartingValuesBVARGROUPPANEL-clone)

------------------------------------------------------------------------

### `StartingValuesBVARGROUPPANEL$new()`

Create new starting values StartingValuesBVARGROUPPANEL

#### Usage

    StartingValuesBVARGROUPPANEL$new(group_allocation = 1:C, C, G = C, N, p, d = 0)

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

Starting values StartingValuesBVARGROUPPANEL

#### Examples

    # starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
    sv = specify_starting_values_bvarGroupPANEL$new(C = 2, N = 3, p = 1)

------------------------------------------------------------------------

### `StartingValuesBVARGROUPPANEL$get_starting_values()`

Returns the elements of the starting values StartingValuesBVARGROUPPANEL
as a `list`.

#### Usage

    StartingValuesBVARGROUPPANEL$get_starting_values()

#### Examples

    # starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
    sv = specify_starting_values_bvarGroupPANEL$new(rep(1,2), C = 2, N = 3, p = 1)
    sv$get_starting_values()   # show starting values as list

------------------------------------------------------------------------

### `StartingValuesBVARGROUPPANEL$set_starting_values()`

Returns the elements of the starting values StartingValuesBVARGROUPPANEL
as a `list`.

#### Usage

    StartingValuesBVARGROUPPANEL$set_starting_values(last_draw)

#### Arguments

- `last_draw`:

  a list containing the same elements as object
  StartingValuesBVARGROUPPANEL

#### Returns

An object of class StartingValuesBVARGROUPPANEL including the last draw
of the current MCMC as the starting value to be passed to the
continuation of the MCMC estimation.

#### Examples

    sv = specify_starting_values_bvarGroupPANEL$new(rep(1,2), C = 2, N = 3, p = 1)

    # Modify the starting values by:
    sv_list = sv$get_starting_values()   # getting them as list
    sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
    sv$set_starting_values(sv_list)      # providing to the class object

------------------------------------------------------------------------

### `StartingValuesBVARGROUPPANEL$clone()`

The objects of this class are cloneable with this method.

#### Usage

    StartingValuesBVARGROUPPANEL$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# starting values for a Bayesian Panel VAR
sv = specify_starting_values_bvarGroupPANEL$new(rep(1,2), C = 2, G = 1, N = 3, p = 1)


## ------------------------------------------------
## Method `StartingValuesBVARGROUPPANEL$new()`
## ------------------------------------------------

# starting values for Bayesian Panel VAR 2-country model with 4 lags for a 3-variable system.
sv = specify_starting_values_bvarGroupPANEL$new(C = 2, N = 3, p = 1)


## ------------------------------------------------
## Method `StartingValuesBVARGROUPPANEL$get_starting_values()`
## ------------------------------------------------

# starting values for a homoskedastic bsvar with 1 lag for a 3-variable system
sv = specify_starting_values_bvarGroupPANEL$new(rep(1,2), C = 2, N = 3, p = 1)
sv$get_starting_values()   # show starting values as list
#> $group_allocation
#> [1] 1 1
#> 
#> $A_c
#> , , 1
#> 
#>               [,1]          [,2]          [,3]
#> [1,]  0.0006047524 -0.0008940336 -0.0010916108
#> [2,] -0.0010688077  0.0011750599  0.0008125073
#> [3,]  0.0007085073  0.0001561209  0.0013138399
#> [4,]  0.0003461859  0.0013191522 -0.0016021002
#> 
#> , , 2
#> 
#>               [,1]          [,2]          [,3]
#> [1,]  0.0006047524 -0.0008940336 -0.0010916108
#> [2,] -0.0010688077  0.0011750599  0.0008125073
#> [3,]  0.0007085073  0.0001561209  0.0013138399
#> [4,]  0.0003461859  0.0013191522 -0.0016021002
#> 
#> 
#> $Sigma_c
#> , , 1
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 4.1962278 1.0119765 0.9181699
#> [2,] 1.0119765 0.3227215 0.1101306
#> [3,] 0.9181699 0.1101306 2.6560153
#> 
#> , , 2
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 4.1962278 1.0119765 0.9181699
#> [2,] 1.0119765 0.3227215 0.1101306
#> [3,] 0.9181699 0.1101306 2.6560153
#> 
#> 
#> $A_g
#> , , 1
#> 
#>               [,1]          [,2]          [,3]
#> [1,]  0.0006047524 -0.0008940336 -0.0010916108
#> [2,] -0.0010688077  0.0011750599  0.0008125073
#> [3,]  0.0007085073  0.0001561209  0.0013138399
#> [4,]  0.0003461859  0.0013191522 -0.0016021002
#> 
#> , , 2
#> 
#>              [,1]          [,2]          [,3]
#> [1,] 9.608396e-04 -0.0012649629 -0.0010374876
#> [2,] 3.590305e-04 -0.0013843557 -0.0004001225
#> [3,] 9.434266e-05  0.0003603438 -0.0004708866
#> [4,] 1.320436e-03  0.0022819910 -0.0013311328
#> 
#> 
#> $Sigma_g
#> , , 1
#> 
#>           [,1]      [,2]      [,3]
#> [1,] 4.1962278 1.0119765 0.9181699
#> [2,] 1.0119765 0.3227215 0.1101306
#> [3,] 0.9181699 0.1101306 2.6560153
#> 
#> , , 2
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  3.1383742 -2.932472  0.7522226
#> [2,] -2.9324717  7.475599 -1.3346634
#> [3,]  0.7522226 -1.334663  6.7560085
#> 
#> 
#> $A
#>               [,1]          [,2]          [,3]
#> [1,]  1.0004495506 -0.0009355935  0.0006422114
#> [2,]  0.0001504719  1.0001053060 -0.0026264663
#> [3,]  0.0007026927  0.0024083851  0.9997600898
#> [4,] -0.0005661534  0.0000210839  0.0006022699
#> 
#> $V
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  0.9029576 -0.3895099 -0.5340549  0.4344648
#> [2,] -0.3895099  3.2317527  1.9467499 -1.3870081
#> [3,] -0.5340549  1.9467499  2.6951423 -1.2111609
#> [4,]  0.4344648 -1.3870081 -1.2111609  2.9850043
#> 
#> $Sigma
#>            [,1]      [,2]       [,3]
#> [1,] 0.78430070 0.2084214 0.05653299
#> [2,] 0.20842138 1.2416036 1.12074847
#> [3,] 0.05653299 1.1207485 2.17284821
#> 
#> $nu
#> [1] 4.1
#> 
#> $m
#> [1] 0.0005246072
#> 
#> $w
#> [1] 0.4979873
#> 
#> $s
#> [1] 0.9800526
#> 


## ------------------------------------------------
## Method `StartingValuesBVARGROUPPANEL$set_starting_values()`
## ------------------------------------------------

sv = specify_starting_values_bvarGroupPANEL$new(rep(1,2), C = 2, N = 3, p = 1)

# Modify the starting values by:
sv_list = sv$get_starting_values()   # getting them as list
sv_list$A <- matrix(rnorm(12), 3, 4) # modifying the entry
sv$set_starting_values(sv_list)      # providing to the class object
```
