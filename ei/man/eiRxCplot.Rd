\name{eiRxCplot.ei}
\alias{eiRxCplot.ei}

\title{
  Plotting Ecological Inference Estimates with eiRxC information
}
\description{
 
}
\usage{
eiRxCplot(ei.object, groups, data, groupnames=groups, threshold=.65, random=FALSE,prop=.15,title="EI RxC Plot",
xaxis="betab",yaxis="betaw",  estimates=TRUE, legendpos=c(.69,1))
}

\arguments{
  \item{ei.object}{An \code{ei} object from the function \code{ei} or \code{ei.sim}.}
  \item{groups}{A list of groups to be represented in the eiRxC plot}
  \item{data}{data where the groups are contained and from which the ei
    object was produced}
  \item{groupnames}{A list of group names corresponding to the groups
    to be printed in the legend}
  \item{threshold}{The threshold where precincts in which a group
    has more than this proportion will have tomography lines colored to represent that group.}
  \item{random}{TRUE to print only a random selection of tomography
    lines}
  \item{prop}{Specifies the proportion of tomography lines to be printed
    if random is TRUE}
  \item{title}{Title for the plot}
  \item{xaxis}{X-axis label for the plot}
  \item{yaxis}{Y-axis label for the plot}
  \item{estimates}{Specifies whether or not the quantities of interest
    should be included in the plot.  If estimates=FALSE, \code{ei.sim}
    is not needed to use \code{eiRxCplot}.
   \item{legendpos}{Specifies coordinate position for the legend.}
}
 
\details{
%%  ~~ If necessary, more details than the description above ~~
}


\references{
Gary King (1997). A Solution to the Ecological Inference Problem.  Princeton: Princeton University Press.
}
\author{
%%  ~~who you are~~
}
\examples{
data(eiRxCsample)
dbuf <- ei("t","x","n",data=eiRxCsample)
dbufsim <- ei.sim(dbuf)
eiRxCplot(dbuf, groups=c("black.ei", "whit.ei", "hisp.ei"), data=eiRxCsample,
groupnames=c("Black", "White", "Hispanic"), title="Non-white vs. White
Turnout", xaxis="Non-white", yaxis="White")
}
