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
#' @examples
#' data(census1910)
#' output <- bounds1(x = census1910$x, t = census1910$t, n = census1910$n)
ei.bounds <- function(x, t, n) {
  bounds1(x, t, n)
}

#' @noRd
bounds1 <- function(x, t, n) {
  # TODO add our checks here in the future!
  bounds_cpp(x, t, n)
}
