\name{plot.ei}
\alias{plot.ei}

\title{Graphic Parameters and Options for  ei model.}

\description{
 Graphs the results taken from the ei output data buffer (i.e., for $
 2\times 2$ tables) and, when it get completed, for ei2 data buffers
 ($ 2\times C$ tables). See summary.ei to extracting the numerical
 information input to the graphics.eigraph(dbuf,"name") or
 plot(x=dbuf,"name").
 See book (EIBAD) \emph{A Solution to the Ecological Inference Problem:
 Reconstructing Individual Behavior from Aggregate Data}, Gary King,
 Princeton University Press, 1997.}


\usage{plot.ei(x,...)}

\arguments{
  \item{x}{A list output of ei}
  
  \item{...}
  {
    A string with the name of the corresponding graph.
    Choose any of the items below in place of name to see the
    corresponding graph. (Items that use true values from
    individual-level data will only work if you have input this
    information into the output buffer \sQuote{x}.
    These are useful for verifying the method if individual-level
    data are available, but they are not useful for most ecological
    inference applications, where this information is not available.
    If no name is provided the default is \sQuote{tomogs}. The options
    are: \newline
    
   \item{beta}{betaB and betaW; betaB vector of length p\eqn{\times 1} with p
	the number of precincts is the mean posterior over
	the simulations of \eqn{\beta_i^b},and similarly for betaW.}
   \item{betaB}{
	density estimate (i.e., a smooth version of a histogram) of point
   	estimates of $ \beta^b_i$'s.}
   \item{betaW}{density estimate (i.e., a smooth version of a histogram)
	of point estimates of $ \beta^w_i$'s.}
   \item{betabw}{lines,Tlines,bivar,Tbivar. See Figure 10.8, page 214.} 
   \item{bias}{combination of biasB, biasW, TbiasB, TbiasW.
	See the points in EIBAD Figure 13.2, page 238.}
  \item{biasB}{precinct $X_i$ by estimated $\beta_i^b$.}
  \item{biasW}{precinct $X_i$ by estimated $ \beta_i^w$.}
  \item{bivar}{estimated $ \beta_i^b$ by $ \beta_i^w$. See the right
    graph in EIBAD, Figure 13.6, page 243.}
  \item{boundX}{combination of boundXB, boundXW, boundXB overlaid
    on TbiasB, boundXW overlaid on TbiasW. See EIBAD Figure 13.2, page
    238.}
  \item{boundXB}{precinct $X_i$ by the bounds on $ \beta_i^b$ (
    each precinct appears as one vertical line), see the lines in the
    left graph in EIBAD Figure 13.2, page 238.}
  \item{boundXW}{precinct $X_i$ by the bounds on $ \beta_i^w$
    (each precinct appears as one vertical line),
    see the lines in the right graph in EIBAD Figure 13.2, page 238.}
  \item{estsims}{All the simulated $ \beta_i^b$'s by all the simulated
    $ \beta_i^w$'s. The simulations should take roughly the shape of the
    mean posterior contours, except for those sampled from outlier
    tomography lines. Homogeneous precincts are excluded.}
  \item{fit}{Combination of xtfit and tomogP.}
  \item{fitT}{Combination of xtfit and tomogT.}
  \item{goodman}{precinct $X_i$ by $T_i$ plot with Goodman's regression
    line plotted. See EIBAD,Figure 4.1, page 60.}
  \item{lines}{precinct $X_i$ by $ T_i$ plot with one estimated line
    per precinct. Homogenous precincts do not appear.
    See EIBAD Figure 10.8, page 214.}
  \item{movie}{For each observation ( $ i=1,\ldots,p$), one page of
    graphics appears with the posterior distribution of $ \beta_i^b$
    and $ \beta_i^w$, a plot of simulated values of
    $ \beta_i^b$ by $ \beta_i^w$ from the tomography line)
    and the numerical values of $ T_i$, $ X_i$, $ N_i$, and $ i$.
    A message is displayed affter the graph for each observation appears,
    which containes the total number of precincts, and that allows
    the user to either view the next observation (hit return),
    jump to a specific observation number
    (type in the number and hit return),see the last
    observation that terminates also the graph, or quit the graph at any
    observation (hit ctrl-c ctrl-c and return).
    See EIBAD, Figure 8.1, page 148.}
  
  \item{movieD}{For each observation ( $ i=1,\ldots,p$), one tomography
    plot appears (like tomogD) with the line for observation $ i$
    darkened.After the graph for each observation appears, the user can
    choose to view the next observation (hit return),jump to a specific
    observation number (type in the number and hit return) including the
    last observation which number is shown in the message, or stop (hit
    ctr-c ctr-c and return). See Figure EIBAD 5.1, page 81.}

  \item{nonpar}{combination of tomogD and a nonparametric density
    estimation with contour plot and surface
    plot representations. See EIBAD Figure 9.9, page 195.}

  \item{post}{combination of postB and postW. See EIBAD Figure 10.4, page 208.}

  \item{postB}{density estimate of (weighted) district quantity of
    interest, the posterior distribution for $ B^b$. See EIBAD Figure 10.4, page 208.}

  \item{postW}{density estimate of (weighted) district quantity of
    interest,the posterior distribution for $ B^w$. See EIBAD Figure 10.4, page 208.}

  \item{prectB}{plot of estimated $ \beta^b_i$ by true $ \beta^b_i$. See
    EIBAD Figure 10.5, page 210.}

  \item{prectW}{plot of estimated $ \beta^w_i$ by true $ \beta^w_i$. See
    EIBAD Figure 10.5, page 210.}

  \item{profile}{Profile plots -- that is plots of the posterior of each
    element of $ \phi$, holding constant
    all other values at their maxima. The maximum posterior point, as
    found by the R function \code{nlminb} is displayed as a vertical
    line from the max down to the bottom of the graph. After seeing
    these plots,you can zoom in by changing eigraph.pro
    (as is useful for seeing how curved the likelihood is near the
    maximum).(Kinks in the profile are due to imprecision in the
    \code{cdfbvn} function necessitating truncation by EcdfTol.)
    This option is very useful for verifying whether the maximization
    procedure did indeed find the maximum.}

  \item{profileR}{Profile plots of the cdfbvn function, (Equation 6.15,
    page 104).These are plots of R by each element of $ \phi$, holding
    constant all other values at their maxima. The maximum posterior
    point, as found by R function \code{nlminb} is displayed as a
    vertical line from the max down to the bottom of the graph.
    After seeing these plots, you can zoom in by changing eigraph.pro.
    Kinks in the profile are due to imprecision in the cdfbvn function
    necessitating truncation by EcdfTol. This option is very useful for
    verifying whether the posterior is numerical stable near the maximum.}

  \item{ptile}{combination of ptileB and ptileW. See EIBAD Figure 10.7, page 213.}

  \item{ptileB}{true percentile at which $ \beta_i^b$ falls by estimated $
    \beta_i^b$.See Figure 10.7, page 213.}

  \item{ptileW}{true percentile at which $ \beta_i^w$ falls by estimated
    $ \beta_i^w$. See Figure 10.7, page 213.}

\item{results}{combination of postB, postW, betaB, betaW.}

\item{sims}{combination of simsB and simsW. See Figure 10.6, page 212.}

\item{simsB}{simulations of $ \beta_i^b$ by true $ \beta_i^b$. See Figure 10.6, page 212.}

\item{simsW}{simulations of $ \beta_i^w$ by true $ \beta_i^w$. See Figure 10.6, page 212.}

\item{TbiasB}{precincts $ X_i$ by true $ \beta_i^b$. See the points in
  the left graph in Figure 13.2, page 238.}

\item{TbiasW}{precincts $X_i$ by true $ \beta_i^w$. See the lines in the
  right graph in Figure 13.2, page 238.}

\item{Tbivar}{true $ \beta_i^b$ by $ \beta_i^w$. See Figure 13.2, page 238.}

\item{Tlines}{precincts $ X_i$ by $ T_i$ plot with one true line
  per precinct. Homogenous precincts do not appear. See Figure 3.1, page 41.}

\item{tomog}{tomography plot with ML contours. See EIBAD Figure 10.2,
  page 204. Dashed lines are added for observations
  Eselect'd out of the estimation stage, and a note appears
  at the bottom left and top right corners if any observations
  are included for which $ T_i=0$ or $ T_i=1$
  (since as tomography lines they are represented as barely visible
  points). The global EselRnd is ignored.}

\item{tomogCI}{tomography plot with 80\% confidence intervals.
  Confidence intervals appear on the screen in red with the remainder
  of the tomography line in yellow. On printed output, which is often
  easier to see than the screen display, the confidence interval
  portion is printed thicker than the rest of the line. See Figure 9.5, page 179.}

\item{tomogCI95}{tomography plot with 95\% confidence intervals.
  Confidence intervals appear on the screen in red with the remainder
  of the tomography line in yellow. On printed output, which is often
  easier to see than the screen display, the confidence interval
  portion is printed thicker than the rest of the line. See Figure 9.5, page 179.}

\item{tomogD}{tomography plot with data only. See Figure EIBAD 5.1, page
  81. Dashed lines are added for observations Eselect'd out of the
  estimation stage,and a note appears at the bottom left and top
  right corners if any observations are included for which $ T_i=0$ or 
  $T_i=1$(since as tomography lines they are represented as barely
  visible points).The global EselRnd is ignored.}

\item{tomogE}{tomography plot with estimated mean posterior
  $ \beta_i^b$ and $ \beta_i^w$ points. Dashed lines are
  added for observations Eselect'd out of the estimation stage,
  and a note appears at the bottom left and top right corners
  if any observations are included for which $ T_i=0$ or $ T_i=1$
  (since as tomography lines they are represented as barely visible
  points).The global EselRnd is ignored.}

\item{tomogP}{tomography plot with mean posterior contours.
  Dashed lines are added for observations Eselect'd out of
  the estimation stage, and a note appears at the bottom
  left and top right corners if any observations are included
  for which $ T_i=0$ or $ T_i=1$ (since as tomography lines
  they are represented as barely visible points). The global EselRnd is ignored.}

\item{tomogS}{combination of tomog, tomogp, tomogCI,
  and Tbivar (or estsims if truth isn't available).
  See EIBAD Figures 10.2 (page 204), 9.5 (page 179), 7.1 (page 126).
  Dashed lines are added for observations Eselect'd out of
  the estimation stage, and a note appears at the bottom
  left and top right corners if any observations are
  included for which $ T_i=0$ or $ T_i=1$ (since as tomography
  lines they are represented as barely visible points). The global EselRnd is ignored.}

\item{tomogT}{tomography plot with true $ \beta_i^b$ and $ \beta_i^w$
  points.See EIBAD Figure 7.1 (page 126). Dashed lines
  are added for observations Eselect'd out of the
  estimation stage, and a note appears at the bottom
  left and top right corners if any observations are
  included for which $ T_i=0$ or $ T_i=1$
  (since as tomography lines they are represented as barely visible points). The global EselRnd is ignored.}

\item{truth}{combination of post, precB, precW
  (compare truth to estimates at district and precinct-level).
  See Figures 10.4 (page 208) and 10.5 (page 210).}

\item{xgraph}{a scattercross with data plotted. See Figure 12.1, page 227.}

\item{XgraphC}{a scattercross with data plotted, with size proportional to $ N_i$.}

\item{xt}{basic $ X_i$ by $ T_i$ scatterplot.}

\item{xtC}{basic $ X_i$ by $ T_i$ scatterplot with circles sized
  proportional to $ {N_i}$ or some other variable defined by the global eigraph.circ.}

\item{xtfit}{precincts $ X_i$ by $ T_i$ plot with estimated
  E($T_i\vert X_i$) and conditional 80\% confidence intervals.
  See Figure 10.3, page 206.}

\item{xtfitg}{xtfit with Goodman's regression line superimposed.}
For additional globals see help for eigraph.
}
    }
  
\details{S3method that wraps around the R function eigraph, the user can
  call it as plot(x=dbuf,"name") or plot(x=dbuf) as \sQuote{dbuf} is an
  object of class \sQuote{ei}. 
}

\value{None.}

\references{\url{http://gking.harvard.edu/ei}
 \newline
  \emph{A Solution to the Ecological Inference Problem: Reconstructing
  Individual Behavior from Aggregate Data}, Gary King, Princeton
  University Press, 1997 (EIBAD)}

\author{Gary King \email{king.harvard.edu}}
\note{Elena Villalon \email{evillalon@iq.harvard.edu} translated the
  original Gauss code into R.}

\seealso{See demos; for example: demo(eip.default)}

