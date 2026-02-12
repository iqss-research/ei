# Visualizing EI (bound)

Visualizing EI (bound)

## Usage

``` r
plot_bound(ei.object, options = list())
```

## Arguments

- ei.object:

  The output of
  [`ei()`](https://iqss-research.github.io/ei/reference/ei.md) (it
  should be used with the `truth` argument)

- options:

  The list of options

  - **parameter**: A parameter to plot. It takes either `betab` or
    `betaw`. This option is only for the 2x2 case.

## Value

a ggplot object

## Examples

``` r
# 2x2
data(matproii)
truth <- cbind(matproii$tb, matproii$tw)
suppressMessages({
  ei_res <- ei(formula = t ~ x, total = "n", truth = truth, data = matproii)
})
plot_bound(ei_res, options = list(parameter = "betab"))

plot_bound(ei_res, options = list(parameter = "betaw"))

# RxC
data(RxCdata)
formula <- cbind(turnout, noturnout) ~ cbind(white, black, hisp)
suppressMessages({
  ei_resRxC <- ei(formula, data = RxCdata)
})
plot_bound(ei_resRxC)
```
