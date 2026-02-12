# Computes Analytical Bounds from Accounting Identity

Returns analytical bounds from accounting identity on unknown table
relationships beta_b, beta_w, from known, observed, table marginals, x,
t (and sample size n).

## Usage

``` r
ei.bounds(x, t, n)
```

## Arguments

- x:

  vector of characteristics, e.g. percentage of blacks in each district

- t:

  vector of characteristics, e.g. percentage of people that voted in
  each district

- n:

  size of each observation, e.g. number of voters in each district

## Value

a numeric matrix

## References

Gary King (1997). A Solution to the Ecological Inference Problem.
Princeton: Princeton University Press.

## Author

Gary King \<\<email: king@harvard.edu\>\> and Molly Roberts \<\<email:
molly.e.roberts@gmail.com\>\>

## Examples

``` r
data(census1910)
output <- ei.bounds(x = census1910$x, t = census1910$t, n = census1910$n)
```
