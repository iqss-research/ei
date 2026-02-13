# Run Ecological Inference Simulation

Run Ecological Inference Simulation

## Usage

``` r
ei_sim(data, ndraws = 99, nsims = 100)
```

## Arguments

- data:

  an `ei_tbl` object from
  [`ei_est()`](https://iqss-research.github.io/ei/reference/ei_est.md)

- ndraws:

  integer, default 99. The number of draws.

- nsims:

  integer, default 10. The number of simulations with each draw.

## Value

ei_tbl

## Examples

``` r
data(sample_ei)
dbuf <- ei_est(sample_ei, x, t, n) %>% ei_sim()
#> ℹ Maximizing likelihood
#> ⠙ Beginning importance sampling.
#> ✔ Beginning importance sampling. [879ms]
#> 
```
