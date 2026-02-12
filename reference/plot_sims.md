# Visualizing EI (simulation)

Visualizing EI (simulation)

## Usage

``` r
plot_sims(ei.object)
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
suppressMessages({
  ei_res <- ei(formula = t ~ x, total = "n", data = matproii)
})
plot_sims(ei_res)
```
