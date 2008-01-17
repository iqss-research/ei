###DESCRIPTION: Basic probabilities functions in Gauss
###             and their R counterparts; x is a matrix
###
cdfni <- function(x){qnorm(x)} ###inverse cumulative density integrated   
cdfn <- function(x){pnorm(x)} ### cumulative density integrated 
cdfnc <- function(x) {1- pnorm(x)}
pdf <- function(x){return(dnorm(x))}
###random numbers drawn from normal distribution
rndn <- function(r, c){
  return(matrix(rnorm(r*c, mean=0, sd=1), nrow=r, ncol=c, byrow=T))}
###random numbers drawn from the uniform distributrion 
rndu <- function(r, c){matrix(runif(r*c), nrow=r, ncol=c, byrow=TRUE)}

###DESCRIPTION Computes the cdf of the standardized bivariate normal
###            with lower limits in -Inf, i.e. lower tail. 
###            x and t are the upper limits for the two variables
###            and rho is the correlation coefficients
###            Wraps pmvnorm of of package mvtnorm
###            
cdfbvn <- function(x,t,rho, maxpts=25000, abseps=0.001, releps=0){
  if(!require(mvtnorm))
    stop("ei:To compute bivariate normal you need to install package mvtnorm")
  
  
  v  <- c(as.vector(x), as.vector(t))
  ln <- length(v)
  low <- rep(-Inf, ln)
 
  rho <- diag(ln)*rho
  p00 <- pmvnorm(lower=low, upper=v,mean=rep(0, ln), sigma=rho, maxpts=maxpts,abseps=abseps, releps=releps);
  return(p00)
 }
###DESCRIPTION As in the Gauss function based on cdfbvn or bivariate normal.
###
cdfbvn2 <- function(h,dh,k,dk,r){
  
y <- cdfbvn(h+dh, k+dk,r)+cdfbvn(h,k,r) - cdfbvn(h,k+dk,r) - cdfbvn(h+dh, k, r)
return(y)
}

###Computes the incomplete gamma function: uses pgamma

cdfgam <- function(x,intlim){
  if(any(x <= 0) || any(intlim < 0))
    stop("cdfgam: requires x > 0 and intlim >= 0")
  fn <- function(x, intlim){
    res <- pgamma(q=intlim, shape=x, rate=1,scale=1, lower.tail= TRUE,log.p=FALSE)
    res <- res - pgamma(q=0, shape=x, rate=1,scale=1, lower.tail= TRUE,log.p=FALSE)
    return(as.matrix(res))
  }
  if(length(x) <= 1)
    return(fn(x,intlim))
 
  lst <- lapply(as.vector(x), function(v){ return(fn(v, intlim))})
  mat <- matrix(,nrow=length(intlim), ncol=length(x))
  for(n in 1:length(lst))
    mat[, n ] <- lst[[n]]

  return(mat)
}

test.cdfgam <- function( x=c(0.5,1,3,10), intlim=c(0,0.2,0.4,0.6,0.8,1)){
  res <- c( 0.0000000, 0.0000000, 0.000000000, 0.000000e+00,
            0.4729107, 0.1812692, 0.001148481, 2.353069e-14,
            0.6289066, 0.3296800, 0.007926332, 2.009810e-11, 
            0.7266783, 0.4511884, 0.023115288, 9.669718e-10, 
            0.7940968, 0.5506710, 0.047422596, 1.433100e-08, 
            0.8427008, 0.6321206, 0.080301397, 1.114255e-07)
  mat <- matrix(res, nrow = length(intlim), ncol=length(x), byrow=TRUE)
  r <- cdfgam(x,intlim) -mat
 if (any(abs(r) >= 1.e-6))
  message("problems with cdfgam")
 return(r)
  
}

##/* INVCDFG - Inverse of the Gamma Cumulative Distribution function
##**           with shape parameter A
##**
##** Usage: x = invcdfg(p,a)
##**
##** Input:  P  - matrix of percentage points ( 0 < p < 1 )
##**         A  - matrix of shape parameters (conformable with p)
##**
##** Output: X  - matrix of critical values st Pr(x < X) = P and x ~ Gamma(A)
##*/
invcdfg <- function(p,a){return(qgamma(p,a))}

invcdfgGary <- function(p,a) {
   ###  local converge,negative,tol,tol2,x0,x1,f0,df0,k ;
     if(! (p >= 0 &&  p <= 1))
        stop("ERROR: INVCDFG - P not in range (0,1)")
       
     if(!(a > 0))
        stop( "ERROR: INVCDFG - A is not positive")
     
     tol  <- 1e-8 ;
     tol2 <- 1e-20 ;
     x0 <- a%dot*% matrix(1,nrow=rows(p),ncol=cols(p))
     converge <- k <- 0
     count <- 0
     while( converge<= 50 || k <= 50 ) {
       f0  = cdfgam(a,x0) 
       df0 = pdfgam(x0,a) 
       if (!(df0 > tol2))
         df0 <- df0 + (df0 < tol2)%dot*%(tol2 - df0) 
        
       x1 <- x0 + (p - f0)%dot/%df0 
       negative <- !all(x1 > tol) 
       if (negative) 
         x1 <- x1 + (x1 <=tol)%dot*%(x0%dot*%(0.5 + 1.5%dot*%(p> f0)) - x1) 
      
        converge <- abs(x0 - x1) < tol &  !negative 
        x0 <- x1 
        k <- k + 1
       count <- count + 1
       if(count >= 100)
         break
     }
     if(!converge)
        stop("Warning: INVCDFG has not converged ")
     
    ### ndpclex###COMPUTE
     return(x0) 
   }

###/* RNDGAM - Random numbers from a General Gamma Distribution
##**
##** Usage:  x = rndgam(nr,nc,a,b,c)
##**
##** Input:  NR - scalar of row size returned matrix
##**         NC - scalar of column size returned matrix
##**         A  - matrix of shape parameters for Gamma distribution
##**         B  - matrix of scale parameters for Gamma distribution
##**         C  - matrix of location parameters for Gamma distribution
##**             (ExE conformable with RxC matrix)
##**
##** Output: X - R X C matrix of random numbers ~ Gamma(a,b,c)
##*/

rndgam <- function(nr,nc,a,b,c) {
     if( !(a > 0)) 
        stop( "ERROR: RNDG - Shape parameter, A, is not positive" )
    
     rndu <- matrix(runif(nr*nc), nrow=nr, ncol=nc)
     res  <- b%dot*%qgamma(rndu,a) + c 
}

### /* PDFGAM - Standardized Gamma Probability Density function
##**
##** Usage:  d = pdfgam(x,a)
##**
##** Input:  X - matrix of Gamma values
##**         A - matrix of shape parameter for Gamma distribution
##**
##** Output: D - matrix of density of Gamma(a) function at X
##*/
pdfgam <- function(x, a){return(dgamma(x,a))}
     
pdfgamGary <- function(x,a) {
     if(!all(a > 0) )
        stop("ERROR: PDFGAM - Shape parameter is not positive")
    
    if( x < 100 )
       return(res <- x^(a-1)%dot*% exp(-x) %dot/% gamma(a)) ;
    
    
     return(res <- exp((a-1)%dot*%log(x) - x - lfactorial(a-1))) ;
    
   }
