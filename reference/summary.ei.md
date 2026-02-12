# Summarize Ecological Inference Estimates

`summary' method for the class `ei'.

## Usage

``` r
# S3 method for class 'ei'
summary(object, ...)
```

## Arguments

- object:

  An `ei` object from the function `ei`.

- ...:

  A list of options to return in graphs. See values below.

## Value

formatted summary object

## References

Gary King (1997). A Solution to the Ecological Inference Problem.
Princeton: Princeton University Press.

## Author

Gary King \<\<email: king@harvard.edu\>\> and Molly Roberts \<\<email:
molly.e.roberts@gmail.com\>\>

## Examples

``` r
data(sample_ei)
formula <- t ~ x
dbuf <- ei(formula = formula, total = "n", data = sample_ei)
#> ℹ Running 2x2 ei
#> ℹ Maximizing likelihood for `erho` = 0.5.
#> ℹ Running 2x2 ei
#> ✔ Running 2x2 ei [1ms]
#> 
#> ⠙ Beginning importance sampling.
#> ✔ Beginning importance sampling. [64ms]
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
#> 19
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
#>         BB        BW        SB        SW       RHO
#>  0.1948911 0.7254949 0.1076914 0.1534774 0.3826016
#> 
#> 
#> ── Truncated psi's (ultimate scale) ──
#> 
#>        BB        BW         SB        SW       RHO
#>  0.203037 0.7147618 0.09824348 0.1394448 0.3389257
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
#> Bb 0.2030510 0.01126210
#> Bw 0.7182539 0.00888656
#> 
#> 
#> ── Precision ──
#> 
#> 4
#> 
```
