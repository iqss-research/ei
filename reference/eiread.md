# Quantities of Interest from Ecological Inference Estimation

`eiread` is the command that pulls quantities of interest from the `ei`
object. The command returns a list of quantities of interest requested
by the user.

## Usage

``` r
eiread(ei.object, ...)
```

## Arguments

- ei.object:

  An `ei` object from the function `ei`.

- ...:

  A list of quantities of interest for `eiread()` to return. See values
  below.

## Value

- betab:

  \\p\\ x \\1\\ point estimate of \\\beta_i^b\\ based on its mean
  posterior. See section 8.2

- betaw:

  \\p\\ x \\1\\ point estimate of \\\beta_i^w\\ based on its mean
  posterior. See section 8.2

- sbetab:

  \\p\\ x \\1\\ standard error for the estimate of \\\beta_i^b\\, based
  on the standard deviation of its posterior. See section 8.2

- sbetaw:

  \\p\\ x \\1\\ standard error for the estimate of \\\beta_i^w\\, based
  on the standard deviation of its posterior. See section 8.2

- phi:

  Maximum posterior estimates of the CML

- psisims:

  Matrix of random simulations of \\\psi\\. See section 8.2

- bounds:

  \\p\\ x \\4\\: bounds on \\\beta_i^b\\ and \\\beta_i^w\\, lowerB ~
  upperB ~ lowerW ~ upperW. See Chapter 5.

- abounds:

  \\2\\ x \\2\\: aggregate bounds rows:lower, upper; columns: betab,
  betaw. See Chapter 5.

- aggs:

  Simulations of district-level quantities of interest \\\hat{B^b}\\ and
  \\\hat{B^w}\\. See Section 8.3.

- maggs:

  Point estimate of 2 district-level parameters, \\\hat{B^b}\\ and
  \\\hat{B^w}\\ based on the mean of aggs. See Section 8.3.

- VCaggs:

  Variance matrix of 2 district-level parameters, \\\hat{B^b}\\ and
  \\\hat{B^w}\\. See Section 8.3.

- CI80b:

  \\p\\ x \\2\\: lower~upper \\80\\\\ confidence intervals for
  \\\beta_i^b\\. See section 8.2.

- CI80w:

  \\p\\ x \\2\\: lower~upper \\80\\\\ confidence intervals for
  \\\beta_i^w\\. See section 8.2.

- eaggbias:

  Regressions of estimated \\\beta_i^b\\ and \\\beta_i^w\\ on a constant
  term and \\X_i\\.

- goodman:

  Goodman's Regression. See Section 3.1

numeric values

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
#> ✔ Beginning importance sampling. [106ms]
#> 
eiread(dbuf, "phi")
#> [1] -1.1534  0.8214 -2.2861 -1.8817  0.4128  0.0000  0.0000
eiread(dbuf, "betab", "betaw")
#> $betab
#>  [1] 0.1911 0.1363 0.1832 0.2871 0.1431 0.1233 0.1053 0.1023 0.2180 0.2627
#> [11] 0.1720 0.2519 0.1738 0.2002 0.2070 0.1968 0.2827 0.1726 0.1528 0.1040
#> [21] 0.3008 0.2336 0.2374 0.2020 0.4004 0.1581 0.1619 0.4107 0.2098 0.2109
#> [31] 0.2145 0.1875 0.1769 0.2718 0.2213 0.1591 0.1964 0.3695 0.2268 0.2679
#> [41] 0.1775 0.1932 0.2113 0.1094 0.1211 0.1157 0.2139 0.2224 0.1336 0.1303
#> [51] 0.2737 0.1775 0.2248 0.2084 0.2441 0.1766 0.2347 0.2507 0.2837 0.1844
#> [61] 0.2125 0.1453 0.2329 0.1628 0.1307 0.1632 0.1288 0.3607 0.2255 0.0372
#> [71] 0.2745 0.1582 0.2014 0.1459 0.3110
#> 
#> $betaw
#>  [1] 0.6914 0.6705 0.6246 0.8805 0.5565 0.3792 0.5204 0.5558 0.7853 0.8980
#> [11] 0.5835 0.9260 0.5962 0.7383 0.7636 0.7106 0.7597 0.6939 0.6730 0.6555
#> [21] 0.9008 0.8741 0.7960 0.7177 0.8375 0.6332 0.5682 0.9458 0.7491 0.7636
#> [31] 0.7722 0.6966 0.7047 0.9280 0.7730 0.5919 0.6778 0.8430 0.7418 0.8159
#> [41] 0.6548 0.7047 0.7314 0.6435 0.5402 0.5015 0.7452 0.7621 0.6622 0.6371
#> [51] 0.8685 0.6601 0.7511 0.7474 0.8836 0.6783 0.8544 0.7672 0.8354 0.7116
#> [61] 0.7259 0.6592 0.7503 0.5929 0.5955 0.6691 0.5537 0.8982 0.7110 0.5022
#> [71] 0.7642 0.6659 0.6798 0.4695 0.9731
#> 
```
