# Plotting Ecological Inference Estimates

`plot' method for the class `ei'.

## Usage

``` r
# S3 method for class 'ei'
plot(x, ...)
```

## Arguments

- x:

  An `ei` object from the function `ei`.

- ...:

  A list of options to return in graphs. See values below.

## Value

- tomogD:

  Tomography plot with the data only. See Figure 5.1, page 81.

- tomog:

  Tomography plot with ML contours. See Figure 10.2, page 204.

- tomogCI:

  Tomography plot with \\80\\\\ confidence intervals. Confidence
  intervals appear on the screen in red with the remainder of the
  tomography line in yellow. The confidence interval portion is also
  printed thicker than the rest of the line. See Figure 9.5, page 179.

- tomogCI95:

  Tomography plot with \\95\\\\ confidence intervals. Confidence
  intervals appear on the screen in red with the remainder of the
  tomography line in yellow. The confidence interval portion is also
  printed thicker than the rest of the line. See Figure 9.5, page 179.

- tomogE:

  Tomography plot with estimated mean posterior \\\beta_i^b\\ and
  \\\beta_i^w\\ points.

- tomogP:

  Tomography plot with mean posterior contours.

- betab:

  Density estimate (i.e., a smooth version of a histogram) of point
  estimates of \\\beta_i^b\\'s with whiskers.

- betaw:

  Density estimate (i.e., a smooth version of a histogram) of point
  estimates of \\\beta_i^w\\'s with whiskers.

- xt:

  Basic \\X_i\\ by \\T_i\\ scatterplot.

- xtc:

  Basic \\X_i\\ by \\T_i\\ scatterplot with circles sized proportional
  to \\N_i\\.

- xtfit:

  \\X_i\\ by \\T_i\\ plot with estimated \\E(T_i\|X_i)\\ and conditional
  \\80\\\\ confidence intervals. See Figure 10.3, page 206.

- xtfitg:

  `xtfit` with Goodman's regression line superimposed.

- estsims:

  All the simulated \\\beta_i^b\\'s by all the simulated
  \\\beta_i^w\\'s. The simulations should take roughly the same shape of
  the mean posterior contours, except for those sampled from outlier
  tomography lines.

- boundXb:

  \\X_i\\ by the bounds on \\\beta_i^b\\ (each precinct appears as one
  vertical line), see the lines in the left graph in Figure 13.2, page
  238.

- boundXw:

  \\X_i\\ by the bounds on \\\beta_i^w\\ (each precinct appears as one
  vertical line), see the lines in the right graph in Figure 13.2, page
  238.

- truth:

  Compares truth to estimates at the district and precinct-level.
  Requires `truth` in the `ei` object. See Figures 10.4 (page 208) and
  10.5 (page 210).

- movieD:

  For each observation, one tomography plot appears with the line for
  the particular observation darkened. After the graph for each
  observation appears, the user can choose to view the next observation
  (hit return), jump to a specific observation number (type in the
  number and hit return), or stop (hit "s" and return).

- movie:

  For each observation, one page of graphics appears with the posterior
  distribution of \\\beta_i^b\\ and \\\beta_i^w\\ and a plot of the
  simulated values of \\\beta_i^b\\ and \\\beta_i^w\\ from the
  tomography line. The user can choose to view the next observation (hit
  return), jump to a specific observation number (type in the number and
  hit return), or stop (hit “s" and return).

a base plot

## Details

Returns any of a set of possible graphical objects, mirroring those in
the examples in King (1997). Graphical option `lci` is a logical value
specifying the use of the Law of Conservation of Ink, where the implicit
information in the data is represented through color gradients, i.e. the
color of the line is a function of the length of the tomography line.
This can be passed as an argument and is used for `tomogD'' and `tomog”
plots.

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
#> ✔ Beginning importance sampling. [56ms]
#> 
plot(dbuf, "tomog")
#> Warning: `plot()` was deprecated in ei 2.0.0.
#> ℹ Please check our reference for new functions:
#>   https://iqss-research.github.io/ei/reference/index.html
plot(dbuf, "tomog", "betab", "betaw", "xtfit")
```
