# Ecological Inference Estimation

`ei` is the main command in the package `EI`. It gives observation-level
estimates (and various related statistics) of \\\beta_i^b\\ and
\\\beta_i^w\\ given variables \\T_i\\ and \\X_i\\ (\\i=1,...,n\\) in
this accounting identity: \\T_i=\beta_i^b\*X_i + \beta_i^w\*(1-X_i)\\.
Results are stored in an `ei` object, that can be read with
[`summary()`](https://rdrr.io/r/base/summary.html) or
[`eiread()`](https://iqss-research.github.io/ei/reference/eiread.md) and
graphed in [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Usage

``` r
ei(
  formula,
  total = NULL,
  Zb = 1,
  Zw = 1,
  id = NA,
  data,
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

- formula:

  A formula of the form \\t ~x\\ in the \\2x2\\ case and
  \\cbind(col1,col2,...) ~ cbind(row1,row2,...)\\ in the RxC case.

- total:

  \`total' is the name of the variable in the dataset that contains the
  number of individuals in each unit

- Zb:

  \\p\\ x \\k^b\\ matrix of covariates or the name of covariates in the
  dataset

- Zw:

  \\p\\ x \\k^w\\ matrix of covariates or the name of covariates in the
  dataset

- id:

  `id' is the name of the variable in the dataset that identifies the precinct. Used for `movie'
  and \`movieD' plot functions.

- data:

  data frame that contains the variables that correspond to formula. If
  using covariates and data is specified, data should also contain `Zb`
  and `Zw`.

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

`ei` object

## Details

The `EI` algorithm is run using the `ei` command. A summary of the
results can be seen graphically using `plot(ei.object)` or numerically
using `summary(ei.object)`. Quantities of interest can be calculated
using `eiread(ei.object)`.

## References

Gary King (1997). A Solution to the Ecological Inference Problem.
Princeton: Princeton University Press.

## Author

Gary King \<\<email: king@harvard.edu\>\> and Molly Roberts \<\<email:
molly.e.roberts@gmail.com\>\>

## Examples

``` r
data(sample_ei)
form <- t ~ x
dbuf <- ei(form, total = "n", data = sample_ei)
#> ℹ Running 2x2 ei
#> ℹ Maximizing likelihood for `erho` = 0.5.
#> ℹ Running 2x2 ei
#> ✔ Running 2x2 ei [2ms]
#> 
#> ⠙ Beginning importance sampling.
#> ✔ Beginning importance sampling. [324ms]
#> 
summary(dbuf)
#> 
#> ── Summary ─────────────────────────────────────────────────────────────────────
#> 
#> ── Erho ──
#> 
#> 0.5
#> 
#> 
#> ── Esigma ──
#> 
#> 0.5
#> 
#> 
#> ── Ebeta ──
#> 
#> 0.5
#> 
#> 
#> ── N ──
#> 
#> 75
#> 
#> 
#> ── Resamp ──
#> 
#> 23
#> 
#> 
#> ── Maximum likelihood results in scale of estimation (and se's) ──
#> 
#>         Bb0       Bw0       sigB       sigW       rho Zb0 Zw0
#>  -1.1533845 0.8214511 -2.2860120 -1.8816353 0.4128456   0   0
#>   0.1010396 0.1017911  0.2489234  0.1575911 0.3654937   0   0
#> 
#> 
#> ── Untruncated psi's ──
#> 
#>         BB       BW        SB        SW       RHO
#>  0.1992554 0.724253 0.1053115 0.1549569 0.3513778
#> 
#> 
#> ── Truncated psi's (ultimate scale) ──
#> 
#>         BB        BW       SB        SW       RHO
#>  0.2049919 0.7125456 0.097689 0.1426185 0.2853098
#> 
#> 
#> ── Aggregate Bounds ──
#> 
#>            betab     betaw
#> lower 0.09549504 0.3979843
#> upper 0.60893461 0.8031228
#> 
#> 
#> ── Estimates of Aggregate Quantities of Interest ──
#> 
#>         mean         sd
#> Bb 0.2046378 0.01172334
#> Bw 0.7170018 0.00925051
#> 
#> 
#> ── Precision ──
#> 
#> 4
#> 
```
