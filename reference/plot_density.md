# Visualizing EI (density)

Visualizing EI (density)

## Usage

``` r
plot_density(ei.object, options = list(parameter = "betab"))
```

## Arguments

- ei.object:

  The output of
  [`ei()`](https://iqss-research.github.io/ei/reference/ei.md)

- options:

  The list of options

  - **parameter**: A parameter to plot. It takes either `betab` or
    `betaw`.

## Value

a ggplot object

## Examples

``` r
data(matproii)
suppressMessages({
  ei_res <- ei(formula = t ~ x, total = "n", data = matproii)
})
plot_density(ei_res, options = list(parameter = "betab"))

plot_density(ei_res, options = list(parameter = "betaw"))
```
