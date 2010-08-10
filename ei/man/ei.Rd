\name{ei}
\alias{EI}
\alias{ei}
\docType{function}
\title{
Ecological Inference Estimation
}
\description{
\code{ei} is the main command in the package _ei_.  It gives observation-level estimates (and various related statistics) of \eqn{\beta_i^b} and \eqn{\beta_i^w} given variables \eqn{T_i} and \eqn{X_i} (\eqn{i=1,...,n}) in this accounting identity: \eqn{T_i=\beta_i^b*X_i + \beta_i^w*(1-X_i)}.  Results are stored in and \code{ei} object, that can be read with \code{summary()} or \code{eiread()} and graphed in \code{plot()}.
}
\usage{
ei(t, x, n, Zb, Zw, erho=.5, esigma=.5, ebeta=0, ealphab=NA,
 ealphaw=NA, truth=NA)
}

\arguments{
	\item{t}{\eqn{p} x \eqn{1} vector}
	\item{x}{\eqn{p} x \eqn{1} vector}
	\item{n}{\eqn{p} x \eqn{1} vector}
	\item{Zb}{vector of 1's for no covariates or a 
\eqn{p} x \eqn{k^b} matrix of covariates}
	\item{Zw}{vector of 1's for no covariates or a 
\eqn{p} x \eqn{k^w} matrix of covariates}
	\item{erho}{The standard deviation of the normal prior on \eqn{\phi_5} for the correlation. Default=.05.}
	\item{esigma}{The standard deviation of an underlying normal distribution, from which a half normal is constructed as a prior for both \eqn{\breve{\sigma}_b} and \eqn{\breve{\sigma}_w}.}
	\item{ebeta}{Standard deviation of the "flat normal" prior on \eqn{\breve{B}^b} and \eqn{\breve{B}^w}.  The flat normal prior is uniform within the unit square and dropping outside the square according to the normal distribution.  Set to zero for no prior (default). Setting to positive values probabilitically keeps the estimated mode within the unit square. 0.25 is a resonable value to experiment with first.}
	\item{ealphab}{cols(Zb) x 2 matrix of means (in the first column) and standard deviations (in the second) of an independent normal prior distribution on elements of \eqn{\alpha^b}.  If you specify Zb, you should probably specify a prior, at least with mean zero and some variance (default is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^b}).}
	\item{ealphaw}{cols(Zw) x 2 matrix of means (in the first column) and standard deviations (in the second) of an independent normal prior distribution on elements of \eqn{\alpha^w}.  If you specify Zw, you should probably specify a prior, at least with mean zero and some variance (default is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^w}).}
}

\details{
	The EI algorithm is run using the \code{ei} command.  A summary of the results can be seen graphically using \code{plot(eiobject)} or numerically using \code{summary(eiobject)}.  Quantities of interest can be calculated using \code{eiread(eiobject)}.
}

\author{

}
\references{
Gary King (1997). A Solution to the Ecological Inference Problem.  Princeton: Princton University Press.
}

\examples{
print("!!!!!!!")
}
