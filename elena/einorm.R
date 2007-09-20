###/***********************
## UNIVARIATE NORMALS **
###************************

###----------------------------------------------------------------
## cumulative normal distribution (not standardized)
##
## usage:  p = cdfnorm(y,mu,sigma2);
##
## INPUTS:
## mu = mean
##sigma2 = variance
## y = normal variate
## OUTPUT:
## p = Prob(Y<y|mu,sigma2), where y is a realization of the r.v. Y
##

 cdfnorm <- function(y,mu,sigma2,lower.tail=T,log.p=F){
 ### local sigma;
  if (sigma2<=0){
    message("cdfnorm: variance must be positive")
    return(NULL)
  }
  sigma=sqrt(sigma2);
  cdf <- pnorm(y, mean=mu,sd=sigma,lower.tail=lower.tail, log.p=log.p)
  return(cdf);
}
## The inverse normal CDF (not standardized)
## usage:  y = invcdfnm(p,mu,sigma2);
##
## INPUTS:
## mu = mean
## sigma2 = variance
## p = Prob(Y<y|mu,sigma2), where y is a realization
## of the random variable Y
## OUTPUT:
## y = normal variate
##

invcdfnm <- function(p,mu,sigma2, lower.tail=T, log.p=F){
 
  sigma <- sqrt(sigma2);
  qq <- qnorm(p, mean=mu, sd=sigma, lower.tail=lower.tail, log.p=log.p)
  return(qq)
}

