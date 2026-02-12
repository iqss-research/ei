# Visualizing EI (tomography plot for RxC)

A tomography plot for an estimated Ecological Inference model in RxC
data. This function supports the 2x3 case.

## Usage

``` r
plot_tomogRxC(formula, data, total = NULL)
```

## Arguments

- formula:

  A formula of the form `cbind(col1, col2,...)~cbind(row1,row2,...)`

- data:

  data that contains the data that corresponds to the formula

- total:

  \`total' is the name of the variable in the dataset that contains the
  number of individuals in each unit

## Value

a ggplot object

## Examples

``` r
data(RxCdata)
formula <- cbind(turnout, noturnout) ~ cbind(white, black, hisp)
plot_tomogRxC(formula, RxCdata)
```
