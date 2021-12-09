#' @importFrom grDevices chull hcl heat.colors rainbow rgb
#' @importFrom graphics abline image legend lines mtext par points polygon text
#' @importFrom stats coef cor cov density dnorm lm na.omit pnorm quantile runif sd terms.formula
#' @importFrom stats var weighted.mean optim
#' @importFrom utils head
#' @importFrom eiPack bounds ei.MD.bayes
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvnorm rmvnorm pmvnorm
#' @importFrom cubature adaptIntegrate
#' @importFrom sp point.in.polygon
#' @importFrom plotrix draw.circle
#' @importFrom msm rtnorm
#' @importFrom ggplot2 aes
#' @importFrom dplyr %>% .data
#' @importFrom rlang enquo eval_tidy
NULL

## usethis namespace: start
#' @importFrom tibble tibble
#' @useDynLib ei, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

utils::globalVariables(c("b_bounds", "w_bounds"))
