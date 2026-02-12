# Plotting 2x3 Ecological Inference Estimates in 3 dimensions

A tomography plot in 3 dimensions for RxC Ecological Inference data and
an estimated Ecological Inference model in RxC data.

## Usage

``` r
tomogRxC3d(
  formula,
  data,
  total = NULL,
  lci = TRUE,
  estimates = FALSE,
  ci = FALSE,
  level = 0.95,
  seed = 1234,
  color = hcl(h = 30, c = 100, l = 60),
  transparency = 0.75,
  light = FALSE,
  rotate = TRUE
)
```

## Arguments

- formula:

  A formula of the form `cbind(col1, col2,...)~cbind(row1,row2,...)`

- data:

  data that contains the data that corresponds to the formula

- total:

  \`total' is the name of the variable in the dataset that contains the
  number of individuals in each unit

- lci:

  logical value specifying the use of the Law of Conservation of Ink,
  where the implicit information in the data is represented through
  color gradients, i.e. the color of the plane is a function of the area
  of the tomography plane.

- estimates:

  logical value specifying whether the point estimates of \\\beta\\'s
  are included for each observation on the tomography plot.

- ci:

  logical value specifying whether the estimated confidence ellipse is
  included on the tomography plot.

- level:

  numeric value from 0 to 1 specifying the significance level of the
  confidence ellipse; eg. .95 refers to 95% confidence ellipse.

- seed:

  seed value for model estimation.

- color:

  color of tomography planes if lci=F.

- transparency:

  numeric value from 0 to 1 specifying transparency of tomography
  planes; 0 is entirely transparent.

- light:

  logical value specifying whether lights should be included in the rgl
  interface. The inclusion of lights will create shadows in the plot
  that may distort colors.

- rotate:

  logical value specifying whether the plot will rotate for 20 seconds.

## Value

a base plot

## Details

Requires rgl package and rgl viewer.

## References

Gary King (1997). A Solution to the Ecological Inference Problem.
Princeton: Princeton University Press.

## Author

Gary King \<\<email: king@harvard.edu\>\>; Molly Roberts \<\<email:
molly.e.roberts@gmail.com\>\>; Soledad Prillaman \<\<email:
soledadartiz@fas.harvard.edu..
