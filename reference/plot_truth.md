# Visualizing EI (with truth)

Compares truth to estimates at the district and precinct-level. Requires
the `truth` argument in the `ei` object.

## Usage

``` r
plot_truth(ei.object)
```

## Arguments

- ei.object:

  The output of
  [`ei()`](https://iqss-research.github.io/ei/reference/ei.md)

## Value

a ggplot object

## Examples

``` r
data(matproii)
truth <- cbind(matproii$tb, matproii$tw)
suppressMessages({
  ei_res <- ei(formula = t ~ x, total = "n", truth = truth, data = matproii)
})
plot_truth(ei_res)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the ei package.
#>   Please report the issue at <https://github.com/iqss-research/ei/issues>.
```
