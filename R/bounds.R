#' Computes Analytical Bounds from Accounting Identity
#'
#' Returns analytical bounds from accounting identity on unknown table
#' relationships beta_b, beta_w, from known, observed, table marginals, x, t
#' (and sample size n).
#'
#'
#' @param x vector of characteristics, e.g. percentage of blacks in each
#' district
#' @param t vector of characteristics, e.g. percentage of people that voted in
#' each district
#' @param n size of each observation, e.g. number of voters in each district
#'
#' @export
#' @return TODO
#'
#' @author Gary King <<email: king@@harvard.edu>> and Molly Roberts <<email:
#' molly.e.roberts@@gmail.com>>
#' @references Gary King (1997). A Solution to the Ecological Inference
#' Problem.  Princeton: Princeton University Press.
#'
#'
#' @examples
#' data(census1910)
#' output <- bounds1(x = census1910$x, t = census1910$t, n = census1910$n)
#'
bounds1 <- function(x, t, n) {
  # set basic values
  homindx <- NULL
  tx <- NULL
  tomx <- NULL
  LbetaB <- NULL
  UbetaB <- NULL
  LbetaW <- NULL
  UbetaW <- NULL
  omx <- 1 - x
  Nb <- x * n
  Nw <- omx * n
  p <- length(x)
  homoindx <- ifelse(x == 0, 1, 0)
  homoindx <- ifelse(x == 1, 2, homoindx)


  # Heterogenous precincts
  tx <- as.matrix(t / x)
  tomx <- as.matrix(t / omx)
  tomxx <- as.matrix(tx - (omx / x))
  txx <- as.matrix(tomx - x / (1 - x))
  LbetaB <- apply(tomxx, 1, function(x) max(0, x))
  UbetaB <- apply(tx, 1, function(x) min(x, 1))
  LbetaW <- apply(txx, 1, function(x) max(0, x))
  UbetaW <- apply(tomx, 1, function(x) min(x, 1))

  # Homogenously black
  bl <- homoindx == 2
  LbetaB[bl] <- t[bl]
  UbetaB[bl] <- t[bl]
  LbetaW[bl] <- NA
  UbetaW[bl] <- NA


  # Homogenously white
  wh <- homoindx == 1
  LbetaB[wh] <- NA
  UbetaB[wh] <- NA
  LbetaW[wh] <- t[wh]
  UbetaW[wh] <- t[wh]

  return(cbind(LbetaB, UbetaB, LbetaW, UbetaW))
}
