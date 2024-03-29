% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy.R
\name{ei_est}
\alias{ei_est}
\title{Run (tidy) Ecological Inference Estimation}
\usage{
ei_est(
  data,
  t,
  x,
  n,
  id = seq_len(nrow(data)),
  Zb = NULL,
  Zw = NULL,
  erho = 0.5,
  esigma = 0.5,
  ebeta = 0.5,
  ealphab = NA,
  ealphaw = NA,
  truth = NA
)
}
\arguments{
\item{data}{data where \code{x}, \code{t}, \code{total}, \code{Zb}, \code{Zw} are found}

\item{t}{<\code{\link[=dplyr_data_masking]{data-masking}}> column of turnout in data}

\item{x}{<\code{\link[=dplyr_data_masking]{data-masking}}> column of subgroup proportions in data}

\item{n}{<\code{\link[=dplyr_data_masking]{data-masking}}> column of total in data}

\item{id}{<\code{\link[=dplyr_data_masking]{data-masking}}> column of unique ids in data}

\item{Zb}{<\code{\link[=dplyr_tidy_select]{data-masking}}> columns of covariates in data}

\item{Zw}{<\code{\link[=dplyr_tidy_select]{data-masking}}> columns of covariates in data}

\item{erho}{The standard deviation of the normal prior on \eqn{\phi_5} for
the correlation. Numeric vector, used one at a time, in order. Default \code{c(.5, 3, 5)}.}

\item{esigma}{The standard deviation of an underlying normal distribution,
from which a half normal is constructed as a prior for both
\eqn{\breve{\sigma}_b} and \eqn{\breve{\sigma}_w}. Default \eqn{= 0.5}}

\item{ebeta}{Standard deviation of the "flat normal" prior on
\eqn{\breve{B}^b} and \eqn{\breve{B}^w}.  The flat normal prior is uniform
within the unit square and dropping outside the square according to the
normal distribution.  Set to zero for no prior. Setting to positive values
probabilistically keeps the estimated mode within the unit square.
Default\eqn{=0.5}}

\item{ealphab}{cols(Zb) x 2 matrix of means (in the first column) and
standard deviations (in the second) of an independent normal prior
distribution on elements of \eqn{\alpha^b}.  If you specify Zb, you should
probably specify a prior, at least with mean zero and some variance (default
is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^b}).}

\item{ealphaw}{cols(Zw) x 2 matrix of means (in the first column) and
standard deviations (in the second) of an independent normal prior
distribution on elements of \eqn{\alpha^w}.  If you specify Zw, you should
probably specify a prior, at least with mean zero and some variance (default
is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^w}).}

\item{truth}{A length(t) x 2 matrix of the true values of the quantities of
interest.}
}
\value{
ei_tbl
}
\description{
Run (tidy) Ecological Inference Estimation
}
\examples{
data(sample_ei)
dbuf <- ei_est(sample_ei, x, t, n)
}
\concept{tidy}
