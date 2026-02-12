# Simulate EI Solution via Importance Sampling

Simulate EI Solution via Importance Sampling

## Usage

``` r
ei.sim(ei.object, ndraws = 99, nsims = 100)
```

## Arguments

- ei.object:

  `ei` object

- ndraws:

  integer. The number of draws. Default is 99.

- nsims:

  integer. The number of simulations within each draw. Default is 100.

## Value

`ei.sim` object

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
ei_obj <- ei(form, total = "n", data = sample_ei, simulate = FALSE)
#> ℹ Running 2x2 ei
#> ℹ Maximizing likelihood for `erho` = 0.5.
#> ℹ Running 2x2 ei
#> ✔ Running 2x2 ei [1ms]
#> 
sims <- ei.sim(ei_obj)
#> ⠙ Beginning importance sampling.
#> ✔ Beginning importance sampling. [107ms]
#> 
```
