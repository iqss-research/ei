#' @importFrom grDevices chull hcl heat.colors rainbow rgb
#' @importFrom graphics abline image legend lines mtext par points polygon text
#' @importFrom stats coef cor cov density dnorm lm na.omit pnorm quantile runif sd terms.formula
#' @importFrom stats var weighted.mean
#' @importFrom utils head
#' @importFrom eiPack bounds ei.MD.bayes
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvnorm rmvnorm pmvnorm
#' @importFrom cubature adaptIntegrate
#' @importFrom mnormt sadmvn
#' @importFrom sp point.in.polygon
#' @importFrom ucminf ucminf
#' @importFrom plotrix draw.circle
#' @importFrom msm rtnorm
#' @importFrom ggplot2 aes
#' @importFrom dplyr %>% .data
NULL

## usethis namespace: start
#' @importFrom tibble tibble
## usethis namespace: end
NULL

utils::globalVariables(c('b_bounds', 'w_bounds'))
