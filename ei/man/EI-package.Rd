\name{ei}
\alias{EI}
\alias{ei}
\docType{function}
\title{
Ecological Inference Estimation
}
\description{
'ei' is the main command in the package _ei_.  It gives observation-level estimates (and various related statistics) of \eqn{\beta_i^b} and \eqn{\beta_i^w} given variables \eqn{T_i} and \eqn{X_i} (\eqn{i=1,...,n}) in this accounting identity: \eqn{T_i=\beta_i^b*X_i + \beta_i^2*(1-X_i)}.  Results are stored in and ei object, that can be read with summary() or get() and graphed in plot().
}
\usage{
ei(t, x, n, Zb, Zw, erho=.5, esigma=.5, ebeta=0, ealphab=NA, ealphaw=NA, truth=NA)
}

\arguments{
	\item{t}{px1 vector}
	\item{x}{px1 vector}
	\item{n}{px1 vector}
	\item{Zb}{scalar 1 for no covariates or a 
\eqn{p x k^b} matrix of covariates}
	\item{Zw}{scalar 1 for no covariates or a \eqn{p x k^w} matrix of covariates}
	\item{erho}{The standard deviation of the normal prior on \eqn{\phi_5} for the correlation. Default=.05.}
	\item{esigma}{The standard deviation of an underlying normal distribution, from which a half normal is constructed as a prior for both \eqn{\breve{\sigma_b}} and \eqn{\breve{\sigma_w}}.}
	\item{ebeta}{Standard deviation of the "flat normal" prior on \eqn{\breve{\beta}^b} and \eqn{\breve{\beta}^w}.  The flat normal prior is uniform within the unit square and dropping outside the square according to the normal distribution.  Set to zero for no prior (default). Setting to positive values probabilitically keeps the estimated mode within the unit square. 0.25 is a resonable value to experiment with first.}
	\item{ealphab}{cols(Zb)x2 matrix of means (in the first column) and standard deviations (in the second)of an independent normal prior distribution on elements of \eqn{\alpha^b}.  If you specify Zb, you shoudl probably specify a prior, at least with mean zero and some variance (default is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^b}).}
	\item{ealphaw}{cols(Zw)x2 matrix of means (in the first column) and standard deviations (in the second)of an independent normal prior distribution on elements of \eqn{\alpha^w}.  If you specify Zw, you shoudl probably specify a prior, at least with mean zero and some variance (default is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^w}).}
}

\details{
	The ei algorithm is run using the ei command.  A summary of the results can be seen graphically using plot(eiobject) or numerically using summary(eiobject).  Quantities of interest can be calculated using get(eiobject).
}

\author{

}
\references{
Gary King (1997). A Solution to the Ecological Inference Problem.  Princeton: Princton University Press.
}

\examples{

}
