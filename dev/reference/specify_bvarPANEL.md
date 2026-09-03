# R6 Class representing the specification of the BVARPANEL model

The class BVARPANEL presents complete specification for the Bayesian
Panel Vector Autoregression.

## References

Jarocinski (2010). Responses to monetary policy shocks in the east and
the west of Europe: a comparison. *Journal of Applied Econometrics*,
**25**(5), 833-868,
[doi:10.1002/jae.1082](https://doi.org/10.1002/jae.1082) .

Dieppe, Legrand, Roye (2016). The BEAR toolbox, *ECB Working Papers*,
**1934**, [doi:10.2866/292952](https://doi.org/10.2866/292952) .

Zellner, Hong (1989). Forecasting international growth rates using
Bayesian shrinkage and other procedures. *Journal of Econometrics*,
**40**(1), 183–202,
[doi:10.1016/0304-4076(89)90036-5](https://doi.org/10.1016/0304-4076%2889%2990036-5)
.

## Public fields

- `p`:

  a non-negative integer specifying the autoregressive lag order of the
  model.

- `prior`:

  an object PriorBSVAR with the prior specification.

- `data_matrices`:

  an object DataMatricesBVARPANEL with the data matrices.

- `starting_values`:

  an object StartingValuesBVARPANEL with the starting values.

- `adaptiveMH`:

  a vector of four values setting the adaptive MH sampler for nu:
  adaptive rate, target acceptance rate, the iteration at which to start
  adapting, the initial scaling rate

## Methods

### Public methods

- [`BVARPANEL$new()`](#method-BVARPANEL-initialize)

- [`BVARPANEL$set_to_Jarocinski()`](#method-BVARPANEL-set_to_Jarocinski)

- [`BVARPANEL$get_data_matrices()`](#method-BVARPANEL-get_data_matrices)

- [`BVARPANEL$get_prior()`](#method-BVARPANEL-get_prior)

- [`BVARPANEL$get_starting_values()`](#method-BVARPANEL-get_starting_values)

- [`BVARPANEL$get_type()`](#method-BVARPANEL-get_type)

- [`BVARPANEL$set_global2pooled()`](#method-BVARPANEL-set_global2pooled)

- [`BVARPANEL$set_adaptiveMH()`](#method-BVARPANEL-set_adaptiveMH)

- [`BVARPANEL$clone()`](#method-BVARPANEL-clone)

------------------------------------------------------------------------

### `BVARPANEL$new()`

Create a new specification of the Bayesian Panel VAR model BVARPANEL.

#### Usage

    BVARPANEL$new(
      data,
      p = 1L,
      exogenous = NULL,
      stationary = rep(FALSE, ncol(data[[1]])),
      type = rep("real", ncol(data[[1]]))
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

#### Returns

A new complete specification for the Bayesian Panel VAR model BVARPANEL.

------------------------------------------------------------------------

### `BVARPANEL$set_to_Jarocinski()`

Sets the model in line with the specification by Jarocinski (2010) as
presented by Dieppe, Legrand, Roye (2016).

#### Usage

    BVARPANEL$set_to_Jarocinski()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$set_to_Jarocinski()

------------------------------------------------------------------------

### `BVARPANEL$get_data_matrices()`

Returns the data matrices as the DataMatricesBVARPANEL object.

#### Usage

    BVARPANEL$get_data_matrices()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$get_data_matrices()

------------------------------------------------------------------------

### `BVARPANEL$get_prior()`

Returns the prior specification as the PriorBVARPANEL object.

#### Usage

    BVARPANEL$get_prior()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$get_prior()

------------------------------------------------------------------------

### `BVARPANEL$get_starting_values()`

Returns the starting values as the StartingValuesBVARPANEL object.

#### Usage

    BVARPANEL$get_starting_values()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$get_starting_values()

------------------------------------------------------------------------

### `BVARPANEL$get_type()`

Returns the type of the model.

#### Usage

    BVARPANEL$get_type()

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$get_type()

------------------------------------------------------------------------

### `BVARPANEL$set_global2pooled()`

Sets the prior mean of the global autoregressive parameters to the OLS
pooled panel estimator following Zellner, Hong (1989).

#### Usage

    BVARPANEL$set_global2pooled(x)

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

### `BVARPANEL$set_adaptiveMH()`

Sets the parameters of adaptive Metropolis-Hastings sampler for the
parameter nu.

#### Usage

    BVARPANEL$set_adaptiveMH(x)

#### Arguments

- `x`:

  a vector of four values setting the adaptive MH sampler for nu:
  adaptive rate, target acceptance rate, the iteration at which to start
  adapting, the initial scaling rate

#### Examples

    spec = specify_bvarPANEL$new(
       data = ilo_dynamic_panel
    )
    spec$set_adaptiveMH(c(0.6, 0.4, 10, 0.1))

------------------------------------------------------------------------

### `BVARPANEL$clone()`

The objects of this class are cloneable with this method.

#### Usage

    BVARPANEL$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)


## ------------------------------------------------
## Method `BVARPANEL$set_to_Jarocinski()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$set_to_Jarocinski()
#> Setting the model in line with the specification by Jarocinski (2010).


## ------------------------------------------------
## Method `BVARPANEL$get_data_matrices()`
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
## Method `BVARPANEL$get_prior()`
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
## Method `BVARPANEL$get_starting_values()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$get_starting_values()
#> <StartingValuesBVARPANEL>
#>   Public:
#>     A: 0.997754037109657 -0.000423835907009905 -6.1932740730671 ...
#>     A_c: 0.00164229858684123 -0.000139758321123161 0.000729969440 ...
#>     Sigma: 6.20624295852711 0.220604762489741 1.99876110021032 -1.1 ...
#>     Sigma_c: 3.44528607420808 -2.51613085850537 0.399326418676065 0.4 ...
#>     V: 7.07913277311316 -3.91204706711624 -1.51783357946632 -0. ...
#>     clone: function (deep = FALSE) 
#>     get_starting_values: function () 
#>     initialize: function (C, N, p, d = 0) 
#>     m: -0.000314940687504978
#>     nu: 5.1
#>     s: 0.197289542155412
#>     set_starting_values: function (last_draw) 
#>     w: 0.490511207907454


## ------------------------------------------------
## Method `BVARPANEL$get_type()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$get_type()
#> [1] "wozniak"


## ------------------------------------------------
## Method `BVARPANEL$set_global2pooled()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$set_global2pooled()


## ------------------------------------------------
## Method `BVARPANEL$set_adaptiveMH()`
## ------------------------------------------------

spec = specify_bvarPANEL$new(
   data = ilo_dynamic_panel
)
spec$set_adaptiveMH(c(0.6, 0.4, 10, 0.1))
```
