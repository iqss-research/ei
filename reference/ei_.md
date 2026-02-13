# Run (tidy) Ecological Inference Estimation and Simulation

Run (tidy) Ecological Inference Estimation and Simulation

## Usage

``` r
ei_(
  data,
  x,
  t,
  n,
  Zb = NULL,
  Zw = NULL,
  id = NA,
  erho = c(0.5, 3, 5, 0.1, 10),
  esigma = 0.5,
  ebeta = 0.5,
  ealphab = NA,
  ealphaw = NA,
  truth = NA,
  simulate = TRUE,
  ndraws = 99,
  nsims = 100,
  covariate = NULL,
  lambda1 = 4,
  lambda2 = 2,
  covariate.prior.list = NULL,
  tune.list = NULL,
  start.list = NULL,
  sample = 1000,
  thin = 1,
  burnin = 1000,
  verbose = 0,
  ret.beta = "r",
  ret.mcmc = TRUE,
  usrfun = NULL
)
```

## Arguments

- data:

  data where `x`, `t`, `total`, `Zb`, `Zw` are found

- x:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  column of subgroup proportions in data

- t:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  column of turnout in data

- n:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  column of total in data

- Zb:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  columns of covariates in data

- Zw:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  columns of covariates in data

- id:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  column of unique ids in data

- erho:

  The standard deviation of the normal prior on \\\phi_5\\ for the
  correlation. Numeric vector, used one at a time, in order. Default
  `c(.5, 3, 5, .1, 10)`.

- esigma:

  The standard deviation of an underlying normal distribution, from
  which a half normal is constructed as a prior for both
  \\\breve{\sigma}\_b\\ and \\\breve{\sigma}\_w\\. Default \\= 0.5\\

- ebeta:

  Standard deviation of the "flat normal" prior on \\\breve{B}^b\\ and
  \\\breve{B}^w\\. The flat normal prior is uniform within the unit
  square and dropping outside the square according to the normal
  distribution. Set to zero for no prior. Setting to positive values
  probabilistically keeps the estimated mode within the unit square.
  Default\\=0.5\\

- ealphab:

  cols(Zb) x 2 matrix of means (in the first column) and standard
  deviations (in the second) of an independent normal prior distribution
  on elements of \\\alpha^b\\. If you specify Zb, you should probably
  specify a prior, at least with mean zero and some variance (default is
  no prior). (See Equation 9.2, page 170, to interpret \\\alpha^b\\).

- ealphaw:

  cols(Zw) x 2 matrix of means (in the first column) and standard
  deviations (in the second) of an independent normal prior distribution
  on elements of \\\alpha^w\\. If you specify Zw, you should probably
  specify a prior, at least with mean zero and some variance (default is
  no prior). (See Equation 9.2, page 170, to interpret \\\alpha^w\\).

- truth:

  A length(t) x 2 matrix of the true values of the quantities of
  interest.

- simulate:

  default = TRUE:see documentation in `eiPack` for options for RxC ei.

- ndraws:

  integer. The number of draws. Default is 99.

- nsims:

  integer. The number of simulations within each draw. Default is 100.

- covariate:

  see documentation in `eiPack` for options for RxC ei.

- lambda1:

  default = 4:see documentation in `eiPack` for options for RxC ei.

- lambda2:

  default = 2:see documentation in `eiPack` for options for RxC ei.

- covariate.prior.list:

  see documentation in `eiPack` for options for RxC ei.

- tune.list:

  see documentation in `eiPack` for options for RxC ei.

- start.list:

  see documentation in `eiPack` for options for RxC ei.

- sample:

  default = 1000

- thin:

  default = 1

- burnin:

  default = 1000

- verbose:

  default = 0:see documentation in `eiPack` for options for RxC ei.

- ret.beta:

  default = "r": see documentation in `eiPack` for options for RxC ei.

- ret.mcmc:

  default = TRUE: see documentation in `eiPack` for options for RxC ei.

- usrfun:

  see documentation in `eiPack` for options for RxC ei.

## Value

an `ei_tbl`

## Examples

``` r
data(sample_ei)
dbuf <- ei_(sample_ei, x, t, n)
#> ℹ Running 2x2 ei
#> ℹ Maximizing likelihood for `erho` = 0.5.
#> ℹ Running 2x2 ei
#> ✔ Running 2x2 ei [1ms]
#> 
#> ⠙ Beginning importance sampling.
#> ✔ Beginning importance sampling. [86ms]
#> 
```
