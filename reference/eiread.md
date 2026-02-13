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
#> ✔ Beginning importance sampling. [73ms]
#> 
eiread(dbuf, "phi")
#> [1] -1.1534  0.8214 -2.2861 -1.8817  0.4128  0.0000  0.0000
eiread(dbuf, "betab", "betaw")
#> $betab
#>  [1] 0.2052 0.1498 0.1818 0.2915 0.1388 0.1392 0.0982 0.0971 0.2238 0.2450
#> [11] 0.1804 0.2491 0.1743 0.2066 0.2175 0.1986 0.2807 0.1768 0.1500 0.1053
#> [21] 0.3121 0.2355 0.2352 0.2095 0.3998 0.1740 0.1520 0.4011 0.2159 0.2096
#> [31] 0.2203 0.1804 0.1802 0.2421 0.2293 0.1389 0.1921 0.3706 0.2279 0.2632
#> [41] 0.1851 0.1956 0.2117 0.1059 0.1132 0.1204 0.2223 0.2369 0.1383 0.1444
#> [51] 0.2609 0.1764 0.2208 0.2082 0.2382 0.1810 0.2565 0.2469 0.2882 0.1971
#> [61] 0.2136 0.1419 0.2356 0.1717 0.1165 0.1731 0.1409 0.3581 0.2250 0.0452
#> [71] 0.2701 0.1583 0.2015 0.1440 0.3250
#> 
#> $betaw
#>  [1] 0.6908 0.6368 0.6248 0.8766 0.5595 0.3770 0.5301 0.5633 0.7842 0.9020
#> [11] 0.5833 0.9262 0.5961 0.7348 0.7636 0.7098 0.7792 0.6840 0.6917 0.6392
#> [21] 0.8881 0.8741 0.7977 0.7177 0.8398 0.6232 0.5703 0.9530 0.7354 0.7637
#> [31] 0.7695 0.7058 0.6910 0.9281 0.7680 0.6056 0.6778 0.8384 0.7300 0.8214
#> [41] 0.6523 0.7034 0.7313 0.6556 0.5489 0.4983 0.7249 0.7562 0.6497 0.6237
#> [51] 0.8731 0.6612 0.7577 0.7475 0.8843 0.6759 0.8514 0.7777 0.8296 0.7112
#> [61] 0.7171 0.6668 0.7360 0.5914 0.6231 0.6576 0.5494 0.9019 0.7223 0.4757
#> [71] 0.7921 0.6613 0.6798 0.4698 0.9697
#> 
```
