# R6 Class representing the specification of the BVARGROUPPANEL model

The class BVARGROUPPANEL presents complete specification for the
Bayesian Panel Vector Autoregression with county groups. The groups can
be pre-specified, which requires the argument `group_allocation` to be
provided, or estimated, which requires the argument `G` for the number
of groups to be provided and the argument `group_allocation` to be left
empty.

## References

Zellner, Hong (1989). Forecasting international growth rates using
Bayesian shrinkage and other procedures. *Journal of Econometrics*,
**40**(1), 183–202,
[doi:10.1016/0304-4076(89)90036-5](https://doi.org/10.1016/0304-4076%2889%2990036-5)
.

## Super class

`BVARPANEL` -\> `BVARGROUPPANEL`

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

  an object StartingValuesBVARGROUPPANEL with the starting values.

- `adaptiveMH`:

  a vector of four values setting the adaptive MH sampler for nu:
  adaptive rate, target acceptance rate, the iteration at which to start
  adapting, the initial scaling rate

## Methods

### Public methods

- [`BVARGROUPPANEL$new()`](#method-BVARGROUPPANEL-initialize)

- [`BVARGROUPPANEL$set_global2pooled()`](#method-BVARGROUPPANEL-set_global2pooled)

- [`BVARGROUPPANEL$clone()`](#method-BVARGROUPPANEL-clone)

Inherited methods

- `BVARPANEL$get_data_matrices()`
- `BVARPANEL$get_prior()`
- `BVARPANEL$get_starting_values()`
- `BVARPANEL$get_type()`
- `BVARPANEL$set_adaptiveMH()`
- `BVARPANEL$set_to_Jarocinski()`

------------------------------------------------------------------------

### `BVARGROUPPANEL$new()`

Create a new specification of the Bayesian Panel VAR model with country
grouping BVARGROUPPANEL. The groups can be pre-specified, which requires
the argument `group_allocation` to be provided, or estimated, which
requires the argument `G` for the number of groups to be provided and
the argument `group_allocation` to be left empty.

#### Usage

    BVARGROUPPANEL$new(
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

### `BVARGROUPPANEL$set_global2pooled()`

Sets the prior mean of the global autoregressive parameters to the OLS
pooled panel estimator following Zellner, Hong (1989).

#### Usage

    BVARGROUPPANEL$set_global2pooled(x)

#### Arguments

- `x`:

  a vector of four values setting the adaptive MH sampler for nu:
  adaptive rate, target acceptance rate, the iteration at which to start
  adapting, the initial scaling rate

#### Examples

    spec = specify_bvarGroupPANEL$new(
       data = ilo_dynamic_panel[1:5],
       G = 2
    )
    spec$set_global2pooled()

------------------------------------------------------------------------

### `BVARGROUPPANEL$clone()`

The objects of this class are cloneable with this method.

#### Usage

    BVARGROUPPANEL$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
spec = specify_bvarGroupPANEL$new(
   data = ilo_dynamic_panel[1:5],
   G = 2
)
#> Country groupings will be estimated.


## ------------------------------------------------
## Method `BVARGROUPPANEL$set_global2pooled()`
## ------------------------------------------------

spec = specify_bvarGroupPANEL$new(
   data = ilo_dynamic_panel[1:5],
   G = 2
)
#> Country groupings will be estimated.
spec$set_global2pooled()
```
