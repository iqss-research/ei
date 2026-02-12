# Run (tidy) Ecological Inference Estimation

Run (tidy) Ecological Inference Estimation

## Usage

``` r
ei_est(
  data,
  t,
  x,
  n,
  id = seq_len(nrow(data)),
  Zb = NULL,
  Zw = NULL,
  erho = 0.5,
  esigma = 0.5,
  ebeta = 0.5,
  ealphab = NA,
  ealphaw = NA,
  truth = NA
)
```

## Arguments

- data:

  data where `x`, `t`, `total`, `Zb`, `Zw` are found

- t:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  column of turnout in data

- x:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  column of subgroup proportions in data

- n:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  column of total in data

- id:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_data_masking.html)\>
  column of unique ids in data

- Zb:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  columns of covariates in data

- Zw:

  \<[`data-masking`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)\>
  columns of covariates in data

- erho:

  The standard deviation of the normal prior on \\\phi_5\\ for the
  correlation. Numeric vector, used one at a time, in order. Default
  `c(.5, 3, 5)`.

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

## Value

ei_tbl

## Examples

``` r
data(sample_ei)
dbuf <- ei_est(sample_ei, x, t, n)
#> ℹ Maximizing likelihood
```
