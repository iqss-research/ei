% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tomogRxc3d.R
\name{tomogRxC3d}
\alias{tomogRxC3d}
\title{Plotting 2x3 Ecological Inference Estimates in 3 dimensions}
\usage{
tomogRxC3d(
  formula,
  data,
  total = NULL,
  lci = TRUE,
  estimates = FALSE,
  ci = FALSE,
  level = 0.95,
  seed = 1234,
  color = hcl(h = 30, c = 100, l = 60),
  transparency = 0.75,
  light = FALSE,
  rotate = TRUE
)
}
\arguments{
\item{formula}{A formula of the form \code{cbind(col1,
col2,...)~cbind(row1,row2,...)}}

\item{data}{data that contains the data that corresponds to the formula}

\item{total}{`total' is the name of the variable in the dataset that
contains the number of individuals in each unit}

\item{lci}{logical value specifying the use of the Law of Conservation of
Ink, where the implicit information in the data is represented through color
gradients, i.e. the color of the plane is a function of the area of the
tomography plane.}

\item{estimates}{logical value specifying whether the point estimates of
\eqn{\beta}'s are included for each observation on the tomography plot.}

\item{ci}{logical value specifying whether the estimated confidence ellipse
is included on the tomography plot.}

\item{level}{numeric value from 0 to 1 specifying the significance level of
the confidence ellipse; eg. .95 refers to 95\% confidence ellipse.}

\item{seed}{seed value for model estimation.}

\item{color}{color of tomography planes if lci=F.}

\item{transparency}{numeric value from 0 to 1 specifying transparency of
tomography planes; 0 is entirely transparent.}

\item{light}{logical value specifying whether lights should be included in
the rgl interface. The inclusion of lights will create shadows in the plot
that may distort colors.}

\item{rotate}{logical value specifying whether the plot will rotate for 20
seconds.}
}
\value{
a base plot
}
\description{
A tomography plot in 3 dimensions for RxC Ecological Inference data and an
estimated Ecological Inference model in RxC data.
}
\details{
Requires rgl package and rgl viewer.
}
\references{
Gary King (1997). A Solution to the Ecological Inference
Problem.  Princeton: Princeton University Press.
}
\author{
Gary King <<email: king@harvard.edu>>; Molly Roberts <<email:
molly.e.roberts@gmail.com>>; Soledad Prillaman <<email:
soledadartiz@fas.harvard.edu..
}
\concept{deprecated}
