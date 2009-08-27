\name{eigraph}
\alias{eigraph}

\title{Graphic Parameters and Options for  ei model.}

\description{
 See the help for plot.ei as they both perform the same graphics and
 computations. The difference is that plot.ei is a S3method that
 recognizes the output of ei as an object of class ei. 
}

\usage{eigraph(dbuf,str,psiu=NA,...)}

\arguments{
  \item{dbuf}{See help file for plot.ei., argument x}
  \item{str}{See help file for plot.ei, argument \eqn{\ldots}}
  \item{psiu}{Vector of length 5. Reparameterized \eqn{$\phi$} into
    untruncated scale (average for bb bw).Default NA.}
  \item{...}{Optional argumnets are the globals specific
  to the graphics:
  \describe{
  \item{eigraph.bb}{string: xlabel for \eqn{$\beta_i^b$} by \eqn{$\beta_i^w$} plots. Default="betaB".}

\item{eigraph.bbhi}{high end for betaB plots. Default=1.}

\item{eigraph.bblo}{low end for betaB plots. Default=0.}

\item{eigraph.bvsmth}{smoothing parameter for nonparametric
  density estimation; used only for option nonpar.
  Default=0.08. (The same parameter is used in ei
  for nonparametric estimation).}

\item{eigraph.bw}{string: ylabel for $ \beta_i^b$ by $
  \beta_i^w$ plots. Default="betaW".}

\item{eigraph.bwhi}{high end for betaW plots. Default=1.}

\item{eigraph.bwlo}{low end for betaW plots. Default=0.}

\item{eigraph.dbuf}{if 1, create data buffer called eigraph.dbuf
  with px, py, pz inputs to the contour plot
  when option nonpar is chosen.}
\item{eigraph.eval}{number of equally spaced evaluations
  on each axis of the nonparametric density estimation. Default=31.}
\item{eigraph.loess}{1 if show simulated and fitted loess
  for xtfit; 0 if show fit only (default)}
\item{eigraph.pro}{nrow($ \phi$)$ \times 2$.
  The first column is the lower limit and the second
  is the upper limit of the range in which eigraph's
  profile command uses for plotting.
  The default, when set to missing, is to use eiread's ebounds option.}

\item{eigraph.smpl}{Set this in the range (0,1] to randomly
  select this fraction of observations to use in tomography plots.
  This is useful if $ p$ is so large that
  it is difficult to see patterns because of the
  large number of lines. Default=1.}

\item{eigraph.T}{string: ylabel for xt plots. Default="T"}

\item{eigraph.Thi}{High end of T graphs. Default=1.}

\item{eigraph.Tlo}{How end of T graphs. Default=0.}

\item{eigraph.thick}{
    value to add to line thickness parameter (default=1).}

\item{eigraph.X}{
    string: xlabel for xt plots. Default="X".}

\item{eigraph.Xhi}{
    High end of X graphs. Default=1.}

\item{eigraph.Xlo}{
    How end of X graphs. Default=0.}

  \item{eigraphC}{scalar: Multiply by circle size to change sizes.
    Default=1, which means do not change. set to (0,1)
    to reduce circle size and to (1,$ \infty$]
    to increase circle size. Must be greater than zero.}

  \item{EIMetaR}{If x=dbuf is a meta-data buffer
    (i.e., from the ei2.mta output global of ei2),
    this global denotes which of the imputed data buffers
    stored in dbuf should be accessed when running eigraph. Default=1.}

  \item{tomogClr}{Colors for each contour drawn on tomography plots,
    default={ 12, 9, 10, 11, 13, 5 }.}

  \item{tomogPct}{Vector of percentage values at which to draw contours.
    One contour is drawn for each element with this fraction
    of this distribution's volume falling within
    the drawn ellipse. Default={.5, .95}; nrow(tomogClr)
    must be $ \geq$ rows(tomogPct).
  }
}
}
}
\details{Allows users to graphic the output of ei in different
  ways. Computations of different quantities to plot are obtained with summary.}

\value{None.}

\references{\url{http://gking.harvard.edu/ei}}

\author{Gary King \email{gking@harvard.edu} }
\note{Elena Villalon \email{evillalon@iq.harvard.edu} translated the
  original Gauss code into R.
}

\seealso{See demos; for example: demo(eip.default)}

