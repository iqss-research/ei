# Convert to `ei_tbl` objects

Convert to `ei_tbl` objects

## Usage

``` r
as_ei_tbl(x)

ei_as_ei_tbl(ei.object)
```

## Arguments

- x:

  an object to be coerced

- ei.object:

  list-based ei object to convert to tibble-based object

## Value

ei_tbl object

## Examples

``` r
data(sample_ei)
form <- t ~ x
dbuf <- ei(form, total = "n", data = sample_ei)
#> ℹ Running 2x2 ei
#> ℹ Maximizing likelihood for `erho` = 0.5.
#> ℹ Running 2x2 ei
#> ✔ Running 2x2 ei [1ms]
#> 
#> ⠙ Beginning importance sampling.
#> ✔ Beginning importance sampling. [110ms]
#> 
dbuf <- ei_as_ei_tbl(dbuf)
```
