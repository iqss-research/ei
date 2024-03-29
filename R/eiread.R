#' Quantities of Interest from Ecological Inference Estimation
#'
#' \code{eiread} is the command that pulls quantities of interest from the
#' \code{ei} object.  The command returns a list of quantities of interest
#' requested by the user.
#'
#'
#' @param ei.object An \code{ei} object from the function \code{ei}.
#' @param \dots A list of quantities of interest for \code{eiread()} to return.
#' See values below.
#' @return \item{betab}{\eqn{p} x \eqn{1} point estimate of \eqn{\beta_i^b}
#' based on its mean posterior.  See section 8.2} \item{betaw}{\eqn{p} x
#' \eqn{1} point estimate of \eqn{\beta_i^w} based on its mean posterior.  See
#' section 8.2} \item{sbetab}{\eqn{p} x \eqn{1} standard error for the estimate
#' of \eqn{\beta_i^b}, based on the standard deviation of its posterior.  See
#' section 8.2} \item{sbetaw}{\eqn{p} x \eqn{1} standard error for the estimate
#' of \eqn{\beta_i^w}, based on the standard deviation of its posterior.  See
#' section 8.2} \item{phi}{Maximum posterior estimates of the CML}
#' \item{psisims}{Matrix of random simulations of \eqn{\psi}.  See section 8.2}
#' \item{bounds}{\eqn{p} x \eqn{4}: bounds on \eqn{\beta_i^b} and
#' \eqn{\beta_i^w}, lowerB ~ upperB ~ lowerW ~ upperW.  See Chapter 5.}
#' \item{abounds}{\eqn{2} x \eqn{2}: aggregate bounds rows:lower, upper;
#' columns: betab, betaw.  See Chapter 5.} \item{aggs}{Simulations of
#' district-level quantities of interest \eqn{\hat{B^b}} and \eqn{\hat{B^w}}.
#' See Section 8.3.} \item{maggs}{Point estimate of 2 district-level
#' parameters, \eqn{\hat{B^b}} and \eqn{\hat{B^w}} based on the mean of aggs.
#' See Section 8.3.} \item{VCaggs}{Variance matrix of 2 district-level
#' parameters, \eqn{\hat{B^b}} and \eqn{\hat{B^w}}.  See Section 8.3.}
#' \item{CI80b}{\eqn{p} x \eqn{2}: lower~upper \eqn{80\%} confidence intervals
#' for \eqn{\beta_i^b}.  See section 8.2.} \item{CI80w}{\eqn{p} x \eqn{2}:
#' lower~upper \eqn{80\%} confidence intervals for \eqn{\beta_i^w}.  See
#' section 8.2.} \item{eaggbias}{Regressions of estimated \eqn{\beta_i^b} and
#' \eqn{\beta_i^w} on a constant term and \eqn{X_i}.} \item{goodman}{Goodman's
#' Regression.  See Section 3.1}
#'
#' @export
#' @return numeric values
#'
#' @author Gary King <<email: king@@harvard.edu>> and Molly Roberts <<email:
#' molly.e.roberts@@gmail.com>>
#' @references Gary King (1997). A Solution to the Ecological Inference
#' Problem.  Princeton: Princeton University Press.
#'
#' @examples
#' data(sample_ei)
#' formula <- t ~ x
#' dbuf <- ei(formula = formula, total = "n", data = sample_ei)
#' eiread(dbuf, "phi")
#' eiread(dbuf, "betab", "betaw")
eiread <- function(ei.object, ...) {
  function.list <- list(
    "betab" = .betaB, "betaw" = .betaW,
    "phi" = .phi, "sbetab" = .sbetab,
    "sbetaw" = .sbetaw, "psisims" = .psisims,
    "bounds" = .bounds, "CI80b" = .CI80b,
    "CI80w" = .CI80w, "abounds" = .abounds,
    "aggs" = .aggs, "maggs" = .maggs,
    "VCaggs" = .VCaggs, "eaggbias" = .eaggbias,
    "goodman" = .goodman
  )
  dec <- ei.object$precision
  arguments <- list(...)
  results <- list()
  for (arg in arguments) {
    if (arg %in% names(function.list)) {
      results[[arg]] <-
        floor(function.list[[arg]](ei.object) * 10^dec) / 10^dec
    } else {
      results[[arg]] <- NA
    }
  }
  if (length(results) == 1) {
    results <- results[[1]]
  }
  if (length(results) < 1) {
    cli::cli_warn("qi results object is empty")
  }
  return(results)
}
