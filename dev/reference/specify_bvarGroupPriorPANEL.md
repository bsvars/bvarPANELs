# R6 Class representing the specification of the BVARGROUPPRIORPANEL model

The class BVARGROUPPRIORPANEL presents complete specification for the
Bayesian Panel Vector Autoregression with county grouping for global
prior parameters. The groups can be pre-specified, which requires the
argument `group_allocation` to be provided, or estimated, which requires
the argument `G` for the number of groups to be provided and the
argument `group_allocation` to be left empty.

## References

Zellner, Hong (1989). Forecasting international growth rates using
Bayesian shrinkage and other procedures. *Journal of Econometrics*,
**40**(1), 183–202,
[doi:10.1016/0304-4076(89)90036-5](https://doi.org/10.1016/0304-4076%2889%2990036-5)
.

## Public fields

- `p`:

  a non-negative integer specifying the autoregressive lag order of the
  model.

- `G`:

  a non-negative integer specifying the number of country groupings.

- `estimate_groups`:

  a logical value denoting whether the groups are to be estimated.

- `prior`:

  an object PriorBSVAR with the prior specification.

- `data_matrices`:

  an object DataMatricesBVARPANEL with the data matrices.

- `starting_values`:

  an object StartingValuesBVARGROUPPRIORPANEL with the starting values.

- `adaptiveMH`:

  a vector of four values setting the adaptive MH sampler for nu:
  adaptive rate, target acceptance rate, the iteration at which to start
  adapting, the initial scaling rate

## Methods

### Public methods

- [`BVARGROUPPRIORPANEL$new()`](#method-BVARGROUPPRIORPANEL-initialize)

- [`BVARGROUPPRIORPANEL$get_data_matrices()`](#method-BVARGROUPPRIORPANEL-get_data_matrices)

- [`BVARGROUPPRIORPANEL$get_prior()`](#method-BVARGROUPPRIORPANEL-get_prior)

- [`BVARGROUPPRIORPANEL$get_starting_values()`](#method-BVARGROUPPRIORPANEL-get_starting_values)

- [`BVARGROUPPRIORPANEL$get_type()`](#method-BVARGROUPPRIORPANEL-get_type)

- [`BVARGROUPPRIORPANEL$set_global2pooled()`](#method-BVARGROUPPRIORPANEL-set_global2pooled)

- [`BVARGROUPPRIORPANEL$set_adaptiveMH()`](#method-BVARGROUPPRIORPANEL-set_adaptiveMH)

- [`BVARGROUPPRIORPANEL$clone()`](#method-BVARGROUPPRIORPANEL-clone)

------------------------------------------------------------------------

### `BVARGROUPPRIORPANEL$new()`

Create a new specification of the Bayesian Panel VAR model with country
grouping for global prior parameters BVARGROUPPRIORPANEL. The groups can
be pre-specified, which requires the argument `group_allocation` to be
provided, or estimated, which requires the argument `G` for the number
of groups to be provided and the argument `group_allocation` to be left
empty.

#### Usage

    BVARGROUPPRIORPANEL$new(
      data,
      p = 1L,
      exogenous = NULL,
      stationary = rep(FALSE, ncol(data[[1]])),
      type = rep("real", ncol(data[[1]])),
      G = NULL,
      group_allocation = NULL
    )

#### Arguments

- `data`:

  a list with `C` elements of `(T_c+p)xN` matrices with time series
  data.

- `p`:

  a positive integer providing model's autoregressive lag order.

- `exogenous`:

  a `(T+p)xd` matrix of exogenous variables.

- `stationary`:

  an `N` logical vector - its element set to `FALSE` sets the prior mean
  for the autoregressive parameters of the `N`th equation to the random
  walk process, otherwise to white noise.

- `type`:

  an `N` character vector with elements set to "rate" or "real"
  determining the truncation of the predictive density to `[0, 100]` and
  `(-Inf, Inf)` (no truncation) for each of the variables.

- `G`:

  a positive integer specifying the number of country groups. Its
  specification is required if `group_allocation` is not provided and
  the country groups to be estimated.

- `group_allocation`:

  an argument that can be provided as a numeric vector with integer
  numbers denoting group allocations to pre-specify the the country
  groups, in which case they are not estimated, or left empty if the
  country groups are to be estimated.

#### Returns

A new complete specification for the Bayesian Panel VAR model BVARPANEL.

------------------------------------------------------------------------

### `BVARGROUPPRIORPANEL$get_data_matrices()`

Returns the data matrices as the DataMatricesBVARPANEL object.

#### Usage

    BVARGROUPPRIORPANEL$get_data_matrices()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$get_data_matrices()

------------------------------------------------------------------------

### `BVARGROUPPRIORPANEL$get_prior()`

Returns the prior specification as the PriorBVARPANEL object.

#### Usage

    BVARGROUPPRIORPANEL$get_prior()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$get_prior()

------------------------------------------------------------------------

### `BVARGROUPPRIORPANEL$get_starting_values()`

Returns the starting values as the StartingValuesBVARPANEL object.

#### Usage

    BVARGROUPPRIORPANEL$get_starting_values()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$get_starting_values()

------------------------------------------------------------------------

### `BVARGROUPPRIORPANEL$get_type()`

Returns the type of the model.

#### Usage

    BVARGROUPPRIORPANEL$get_type()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$get_type()

------------------------------------------------------------------------

### `BVARGROUPPRIORPANEL$set_global2pooled()`

Sets the prior mean of the global autoregressive parameters to the OLS
pooled panel estimator following Zellner, Hong (1989).

#### Usage

    BVARGROUPPRIORPANEL$set_global2pooled(x)

#### Arguments

- `x`:

  a vector of four values setting the adaptive MH sampler for nu:
  adaptive rate, target acceptance rate, the iteration at which to start
  adapting, the initial scaling rate

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$set_global2pooled()

------------------------------------------------------------------------

### `BVARGROUPPRIORPANEL$set_adaptiveMH()`

Sets the parameters of adaptive Metropolis-Hastings sampler for the
parameter nu.

#### Usage

    BVARGROUPPRIORPANEL$set_adaptiveMH(x)

#### Arguments

- `x`:

  a vector of four values setting the adaptive MH sampler for nu:
  adaptive rate, target acceptance rate, the iteration at which to start
  adapting, the initial scaling rate

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel[1:5]
    )
    spec$set_adaptiveMH(c(0.6, 0.4, 10, 0.1))

------------------------------------------------------------------------

### `BVARGROUPPRIORPANEL$clone()`

The objects of this class are cloneable with this method.

#### Usage

    BVARGROUPPRIORPANEL$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
spec = specify_bvarGroupPriorPANEL$new(
   data = ilo_dynamic_panel,
   G = 2
)
#> Country groupings will be estimated.


## ------------------------------------------------
## Method `BVARGROUPPRIORPANEL$get_data_matrices()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$get_data_matrices()
#> <DataMatricesBVARPANEL>
#>   Public:
#>     Y: list
#>     clone: function (deep = FALSE) 
#>     exogenous: list
#>     get_data_matrices: function () 
#>     initialize: function (data, p = 1L, exogenous = NULL, type = rep("real", 
#>     missing: list
#>     type: real real real real


## ------------------------------------------------
## Method `BVARGROUPPRIORPANEL$get_prior()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$get_prior()
#> <PriorBVARPANEL>
#>   Public:
#>     M: 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0
#>     S_Sigma_inv: 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
#>     S_inv: 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1
#>     W: 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 10
#>     a_w: 1
#>     clone: function (deep = FALSE) 
#>     eta: 5
#>     get_prior: function () 
#>     initialize: function (C, N, p, d = 0, stationary = rep(FALSE, N)) 
#>     lambda: 72
#>     mu_Sigma: 5
#>     mu_m: 1
#>     nu_s: 3
#>     s_s: 1
#>     s_w: 1
#>     sigma2_m: 1


## ------------------------------------------------
## Method `BVARGROUPPRIORPANEL$get_starting_values()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$get_starting_values()
#> <StartingValuesBVARPANEL>
#>   Public:
#>     A: 1.00116165770763 0.000278610576251996 0.0004603944288530 ...
#>     A_c: 0.000961575795264991 0.000465729787011991 0.000901217167 ...
#>     Sigma: 0.957294128849276 -0.894005328637003 0.611378256564174 0 ...
#>     Sigma_c: 2.73878583772208 -3.18613427139964 2.12132858657354 -3.1 ...
#>     V: 7.64473455219002 0.20870502349516 -0.709831326477997 -0. ...
#>     clone: function (deep = FALSE) 
#>     get_starting_values: function () 
#>     initialize: function (C, N, p, d = 0) 
#>     m: -0.000806155390666568
#>     nu: 5.1
#>     s: 0.0493133345701443
#>     set_starting_values: function (last_draw) 
#>     w: 1.97055002989595


## ------------------------------------------------
## Method `BVARGROUPPRIORPANEL$get_type()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$get_type()
#> [1] "wozniak"


## ------------------------------------------------
## Method `BVARGROUPPRIORPANEL$set_global2pooled()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$set_global2pooled()


## ------------------------------------------------
## Method `BVARGROUPPRIORPANEL$set_adaptiveMH()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel[1:5]
)
spec$set_adaptiveMH(c(0.6, 0.4, 10, 0.1))
```
