\name{ei}
\alias{EI}
\alias{ei}
\title{
Ecological Inference Estimation
}
\description{
\code{ei} is the main command in the package \code{EI}.  It gives observation-level estimates (and various related statistics) of \eqn{\beta_i^b} and \eqn{\beta_i^w} given variables \eqn{T_i} and \eqn{X_i} (\eqn{i=1,...,n}) in this accounting identity: \eqn{T_i=\beta_i^b*X_i + \beta_i^w*(1-X_i)}.  Results are stored in an \code{ei} object, that can be read with \code{summary()} or \code{eiread()} and graphed in \code{plot()}.
}
\usage{
ei(formula, total, data Zb=1, Zw=1,data=NA, erho=.5, esigma=.5, ebeta=.5, ealphab=NA,
 ealphaw=NA, truth=NA, precision=4)
}

\arguments{
	\item{t}{\eqn{p} x \eqn{1} vector, or name of vector if data argument is included}
	\item{x}{\eqn{p} x \eqn{1} vector, or name of vector if data argument is included}
	\item{n}{\eqn{p} x \eqn{1} vector, or name of vector if data argument is included}
	\item{Zb}{
\eqn{p} x \eqn{k^b} matrix of covariates}
	\item{Zw}{
\eqn{p} x \eqn{k^w} matrix of covariates}
	\item{data}{optional data frame that if specified should contain \code{t},\code{x},\code{n}s.  If using covariates and data is specified, data should also contain \code{Zb} and \code{Zw}.}
	\item{erho}{The standard deviation of the normal prior on \eqn{\phi_5} for the correlation. Default \eqn{=0.5}.}
	\item{esigma}{The standard deviation of an underlying normal distribution, from which a half normal is constructed as a prior for both \eqn{\breve{\sigma}_b} and \eqn{\breve{\sigma}_w}. Default \eqn{= 0.5}}
	\item{ebeta}{Standard deviation of the "flat normal" prior on \eqn{\breve{B}^b} and \eqn{\breve{B}^w}.  The flat normal prior is uniform within the unit square and dropping outside the square according to the normal distribution.  Set to zero for no prior. Setting to positive values probabilistically keeps the estimated mode within the unit square. Default\eqn{=0.5}}
	\item{ealphab}{cols(Zb) x 2 matrix of means (in the first column) and standard deviations (in the second) of an independent normal prior distribution on elements of \eqn{\alpha^b}.  If you specify Zb, you should probably specify a prior, at least with mean zero and some variance (default is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^b}).}
	\item{ealphaw}{cols(Zw) x 2 matrix of means (in the first column) and standard deviations (in the second) of an independent normal prior distribution on elements of \eqn{\alpha^w}.  If you specify Zw, you should probably specify a prior, at least with mean zero and some variance (default is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^w}).}
	\item{truth}{A length(t) x 2 matrix of the true values of the quantities of interest.}
	\item{precision}{The number of decimal places that are printed in \code{summary} and \code{eiread} functions.}
}

\details{
	The \code{EI} algorithm is run using the \code{ei} command.  A summary of the results can be seen graphically using \code{plot(ei.object)} or numerically using \code{summary(ei.object)}.  Quantities of interest can be calculated using \code{eiread(ei.object)}.
}

\author{

}
\references{
Gary King (1997). A Solution to the Ecological Inference Problem.  Princeton: Princeton University Press.
}

\examples{
data(sample)
form <- t ~ x
dbuf <- ei(form,total="n",data=sample)
summary(dbuf)
}
