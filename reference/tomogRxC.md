# Plotting Ecological Inference Estimates with eiRxC information

A tomography plot for an estimated Ecological Inference model in RxC
data.

## Usage

``` r
tomogRxC(formula, data, total = NULL, refine = 100)
```

## Arguments

- formula:

  A formula of the form `cbind(col1, col2,...)~cbind(row1,row2,...)`

- data:

  data that contains the data that corresponds to the formula

- total:

  \`total' is the name of the variable in the dataset that contains the
  number of individuals in each unit

- refine:

  specifies the amount of refinement for the image. Higher numbers mean
  better resolution.

## References

Gary King (1997). A Solution to the Ecological Inference Problem.
Princeton: Princeton University Press.

## Author

Gary King \<\<email: king@harvard.edu\>\> and Molly Roberts \<\<email:
molly.e.roberts@gmail.com\>\>

## Examples

``` r
data(RxCdata)
formula <- cbind(turnout, noturnout) ~ cbind(white, black, hisp)
tomogRxC(formula, data = RxCdata)
#> Warning: `tomogRxC()` was deprecated in ei 2.0.0.
#> ℹ Please use `plot_tomogRxC()` instead.
#> There are 0 tomography polygons with no information
```
