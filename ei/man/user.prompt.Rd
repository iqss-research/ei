\name{user.prompt}
\alias{user.prompt}

\title{Scan user input}

\description{
  Gets user response to continue with the demos simulation, jump to
  specific precinct or to quit.
}

\usage{user.prompt(p)}

\arguments{
  \item{p}{String or scalar to add to the prompt-messages. Default is NA
    that will not change the prompt-messages}
  }

\details{Allows users to terminate demos early. 
The argument p indicates an optional string to be printed with the messages.}

\value{None.}

\references{\url{http://gking.harvard.edu/ei}}

\author{Elena Villalon \email{evillalon@iq.harvard.edu}}

\seealso{See demos; for example: demo(eip.default)}

