##/*
##**  This archive is part of the program EI
##**  (C) Copyright 1995-2001 Gary King
##**  All Rights Reserved.
##**
##**  various methods of computing areas (or logs of areas) or 
##**  pdf's or cdf's of the univariate and bivariate normal distributions.
##**  and some related distributions
##*/
#include ei.ext;
##/*******************
##** RANDOM NUMBERS **
##********************
##*/
##/*
##   y = rndmn(mu,vc,n);
##**
##** inputs: mu = kx1 means
##**         vc = kxk variance matrix
##**          n = number of simulations
##**
##** output:  y = nxk matrix of dependent Multivariate Normal Random Variables
##**              each row of y is one 1xk simulation
##**
##** example:
##**
##**      c=rndn(30,5);
##**      vc=c'c;                 @ some theoretical var-cov matrix,
##**                                c'c for the example to assure vc ispos def @
##**      mu=0|5|-10|130|3;
##**      y=rndmn(mu,vc,500);
##**      "the theoretical correlation matrix";;
##**          d=1/sqrt(diag(vc));   d.*vc.*d';?;
##**      "simulated correlation matrix";;  corrx(y);
##**      "theoretical mean matrix: ";;     mu';
##**      "simulated mean matrix:   ";;     meanc(y)';
##**
##**  History:  
##**  12/6/94  Added capability to do cholsky decomposition on matrices 
##**           that have symmetric columns and rows of all zeros -- e.g.,
##**           the vc matrix returned from a restricted regression.  Also,
##**           vectorized computation of result.  Curt Signorino 
##*/
rndmn <- function(mu,vc,n,eps=1.e-5){

  k <- rows(mu)
  c <- cols(mu)
  r <- rows(vc)
  if ((r!=k || cols(vc)!= k) && (r!=1))
    stop( "rndmn: mu must be kx1, and vc kxk or scalar") 
      
  if (n<1) 
    stop("rndmn: number of simulations must be >=1   ") 
   
  if (c!=1 && c!=n)
    stop( "rndmn: mu must be kxn or kx1")

  if( scalzero(vc))
    return(t(mu)%dot*%as.matrix(rep(1,n))) 
  tmp <- colSums(as.matrix(dotfeq(vc,0)))
  i <- tmp == r     ##   @ which columns are all zeros?  @      
 
  if (all(i==FALSE)){##   @ no all-zero columns/rows      @
             
    a <- try(chol(vc,pivot=FALSE),silent=TRUE)
    if(class(a)=="try-error"){
      message("rndmn: chol with pivot=F fails trying pivot=T and eichol")
      a <- eichol(vc,pivot=TRUE,mess="rndmn",tol=eps)
    }
    a <- t(a)
               ##   @ matrix square root function   @
  }else{ ###                      @ all-zero columns/rows exist   @
    t <- subset(diag(r), subset=!i)###  @ create transform matrix       @
    vcd <- (t%*%vc)%*%t(t) ###  @ create nonsingular submatrix  @
    ad <- try(chol(vcd,pivot=FALSE),silent=TRUE)
    if(class(ad)=="try-error") {
      message("rndmn: chol with pivot=F fails trying pivot=T and eichol")
      ad <- eichol(vcd,pivot=TRUE,mess="rndmn",tol=eps)
    }
###  @ cholsky decomp of submatrix   @
    a <- (t(t)%*%ad)%*%t ###              @ rebuild full square-root matrix @
  }
  mat <- matrix(rnorm(k*n, mean=0, sd=1), nrow=k, ncol=n, byrow=T)
 
  
  res <- t(mu%plus%(a%*%mat)) ###      @ dep ran normals with mean mu, var vc @
   return(res)
}

rndnm.test <- function(n){
   vcv <- c(28.6072, -2.1666, -6.9943,  8.3078, -0.3080, -2.1666, 36.0326, -8.1411,  8.3657,  1.0486, 
            -6.9943, -8.1411, 37.5934,  5.8055, -3.8859,8.3078,  8.3657,  5.8055, 41.9698, -0.3851,
             -0.3080,  1.0486, -3.8859, -0.3851, 19.6271)
   vc <-  matrix(vcv,ncol=5,byrow=TRUE)
   mu <- as.matrix(c(0,5,-10,130,3))
   yv <- c(-4.5975, -3.4228, -7.1720, 130.2136, 10.1580,
           1.6090,  0.9583, -1.2137, 137.6933, -0.1357, 
           5.1135, -5.6087, -10.2555, 124.3189,  6.6084,
           -0.7678,  3.8019, -7.2281, 133.4664, 10.1348, 
           -4.0171, 12.3191, -8.4431, 128.2330,  5.4298,
            7.2640, -2.3761, -16.1102, 130.8514,  7.3573,
            2.9412,  0.7858, -3.4992, 132.9169, -0.7575,
            0.8296, -2.0460, -9.9832, 125.0697,  7.6052,
           -2.2463, 10.5770, -19.5932, 132.8424, -0.1762,
           -2.8260,  4.8684, -17.2458, 120.3651, 12.7777)
   ytest <- matrix(yv,ncol=5,byrow=TRUE)
   y <- rndmn(mu,vc,n)
   if(n ==10)
     print(y - ytest )
 }

##/*
##   y = rndmt(mu,vc,df,n);
##**
##** inputs: mu = kx1 means
##**         vc = kxk variance matrix
##**         df = scalar degrees of freedom
##**          n = scalar number of simulations
##**
##** output:  y = nxk matrix of dependent Multivariate T Random Variables
##**              each row of y is one 1xk simulation
##**
##*/
rndmt <- function(mu,vc,df,n,eps=1.e-5){
###  local k,c,r,i,t,vcd,ad,a,res;
  k <- rows(mu)
  c <- cols(mu)
  r <- rows(vc)
  if (( r!=k || cols(vc)!=k) && (r!=1))
    stop( "rndmt: mu must be kx1, and vc kxk or scalar")
    
  if( n<1) 
    stop("rndmt: number of simulations must be >=1   ") 

  if( !(c %in% c(1,n)))
    stop("rndmt: mu must be kxn or kx1")
 
  if( vc==0) 
    return(t(as.matrix(mu))%dot*% matrix(1,nrow=n,ncol=1)) 
  
  i <- colSums(as.matrix(dotfeq(vc,0))) ==r ### @ which columns are all zeros?  @
  if( all(!i)){                              ###@ no all-zero columns/rows      @
    a <- try(chol(vc,pivot=FALSE),silent=TRUE)
    if(class(a) == "try-error"){
      message("rndmt: chol with pivot=F fails trying pivot=T and eichol")
      a <- eichol(vc,pivot=TRUE,mess="rndmt",tol=eps)
    }
    a <- t(a)
###@ matrix square root function   @
  }else{                     ###@ all-zero columns/rows exist   @
    t <- subset(diag(r),subset=!i);      ###@ create transform matrix       @
    if(!length(t)) t <- NA
    vcd <- t%*%vc%*%t(t)            ###@ create nonsingular submatrix  @
    ad <- try(chol(vcd,pivot=FALSE), silent=TRUE)
    if(class(ad)=="try-error"){
      message("rndmt: chol with pivot=F fails trying pivot=T and eichol")
      ad <- eichol(vcd,pivot=TRUE,mess="rndmt",tol=eps)
    }
###@ cholsky decomp of submatrix   @
    a <- (t(t) %*% ad) %*% t               ####@ rebuild full square-root matrix @
  }
  rndn <- matrix(rnorm(k*n, mean=0, sd=1), nrow=k, ncol=n, byrow=T)
  res <- t(mu+a%*%(rndn %dot*%sqrt(df%dot/%rndchi(1,n,df))));
  return(res)
}

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
  if (any(sigma2<=0)){
    message("cdfnorm: variance must be positive")
    return(NULL)
  }
  sigma <- sqrt(sigma2);
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

###/* ----------------------------------------------------------------
##
##  y = Lcdfnormi(bounds,mu,sig2);
##
##  y = ln( cdfnorm(bounds[.,2],mu,sig2) - cdfnorm(bounds[.,1],mu,sig2) );
##  but much more numerically precise than doing it this way.
##
##*/
lcdfnormi <- function(bounds,mu,sig2){
 ### local res,a,b,res2,T,c1,c2,c3,c4,ba,c0,r,la,lb,lma,lmb,rs,s,bs,as,m,rb;
  rs <- nrow(as.matrix(mu));
  rb <- nrow(as.matrix(bounds))

  if(rb!=1 && rb!=rs)
    message("lcdfnormi: input error")
  else if( rb==1)
    bounds <- matrix(1,nrow=rs, ncol=2) %dot*% bounds;

  a <- bounds[,1];
  b <- bounds[,2];

  if (any((b-a)<0))
    stop("lcdfnormi: bounds error");
    
  s <- sqrt(sig2);
  bs <- (b-mu)%dot/%s;
  as <- (a-mu)%dot/%s;
  
  m <- NA
  res <- infrv(lncdfn2(as,bs-as),m,m);

  return(res);
}
lncdfn2 <- function(a, b, eps=1e-25){
  if(length(eps))
    .Machine$double.eps <- eps
  res <- log(abs(pnorm(a) - pnorm(a+b)))
  return(res)
}
##
## y = lncdfbvnu(bb,bw,sb,sw,rho);
##
## ln of the cdfbvn on the unit square, with various options
## for different procs to do the computations

lncdfbvnu <- function(bb,bw,sb,sw,rho){
 
  Bb <- bb
  Bw <- bw
 
  evbase <- get("evbase", env=parent.frame())
  Ecdfbvn <- get("Ecdfbvn", env=evbase)
  o <- 1
  if (equalev(Bb,Bb[1]) && equalev(Bw, Bw[1])){
    o <- matrix(1,nrow=rows(Bb),ncol=1)
    Bb <- Bb[1]
    Bw <- Bw[1]
  }else{
    o <- matrix(1,nrow=rows(Bb),ncol=1)
    Bb <- mean(Bb)
    Bw <- mean(Bw)
  }
  
  if(all(rho==0)){
  
    R <- lcdfnormi(cbind(0,1),Bb,sb^2)+lcdfnormi(cbind(0,1),Bw,sw^2)
  
  
  }else{
  
    if (all(Ecdfbvn==1)){
 
      R <- log(abs(cdfbvnormi(Bb,Bw,sb,sw,rho)))
 
    }else if (all(Ecdfbvn==2)){
 
      R <- log(abs(cdfbvnorme(Bb,Bw,sb,sw,rho,evbase)))
 
    }else if(all(Ecdfbvn==3)){
 
      R <- lcdfbvnorma(Bb,Bw,sb,sw,rho)
 
    }else if(all(Ecdfbvn==4)){
 
      R <- log(abs(cdfbvnunit(Bb,Bw,sb,sw,rho,evbase)))
 
    }else if(all(Ecdfbvn==5)){
    
          
      R <- lncdfbvn2i(Bb,Bw,sb,sw,rho)
 
    }else if(all(Ecdfbvn==6)){
 
      R <- log(abs(cdfbvnormig(Bb,Bw,sb,sw,rho)))
 
      }
     
  }

  R <- R%dot*%o


  return(R);
}
###/* ----------------------------------------------------------------
##   y = lnpdfmn(y,mu,vc);
##**
##**  y = kx1
##** mu = kx1 vector of means
##** vc = kxk var covar matrix
##**
##** y = log of the Multivariate Normal Density (scalar)
##**
##** GLOBAL:  _Eivc = invpd(vc) to save computation time (kxk)
##**                  or missing to compute 
##**
##** NOTE:  detl (system global) will be used, so do not call invpd() except
##**        to define _Eivc before calling this proc.
##*/
lnpdfmn <- function(y,mu,vc,Eivc,detl=0){
###  local a,b,c,k,ivc,w,ymu,res;

  if (any(is.na(Eivc)))
    ivc <- invpd(vc)
  else
    ivc <- Eivc;
  
  ymu <- y-mu
  res <- -(rows(mu)*1.83787706640934548+log(detl)+(t(ymu)%*%ivc)%*%ymu)/2
  return(res)
}

###/* another version without scale factor for importance sampling */
lnpdfmn2 <- function(y,mu,vc,Eivc){
 ### local a,b,c,k,ivc,w,ymu,res;

  if(any(is.na(Eivc)))
    ivc <- invpd(vc)
  else
    ivc <- Eivc;
  
  ymu <- y-mu;
  res <- -(t(ymu)%*%ivc)%*%ymu
  res <- res/2
  return(res)
}


###/* ---------------------------------------------------------------
##   y = lnpdfmt(y,mu,vc,df);
##**
##**  y = kx1
##** mu = kx1 vector of means
##** vc = kxk var covar matrix
##** df = degrees of freedom
##**
##** y = log of Multivariate T pdf without normalizing constant (scalar)
##**
##** GLOBAL:  _Eivc = invpd(vc) to save computation time or missing to compute
##**          must be kxk.
##**
##** NOTE:  detl (system global) will be used, so do not call invpd() except
##**        to define _Eivc before calling this proc.
##*/
lnpdfmt <- function(y,mu,vc,df,Eivc,detl=0){
###  local a,b,c,k,ivc,w,ymu;
  k <- rows(mu);

  ymu <- y-mu
  if (any(is.na(Eivc)))
    ivc <- invpd(vc)
  else
    ivc <- Eivc
  
  
  w <- (t(ymu)%*%ivc)%*%ymu;
  b <- -0.5*log(detl)
  c <- -(df+k)%dot*%log(1+w%dot/%df)/2
  return(b+c)
}


##/* ----------------------------------------------------------------
## cdfbvnormig(Bb,Bw,sigb,sigw,rho)
## cdf of the bivariate truncated normal with bounds at 0,1, 0,1
## using Genz's CDFBVNG -- slow, but precise -- with impermissible values
## truncated.
##
## Bb = px1 vector of means of var 1
## Bw = means of var 2
## sigb = stan dev of var 1
## sigw = stan dev of var 2
## rho = correlation

cdfbvnormig <- function(Bb,Bw,sigb,sigw,rho){
  evbase <- get("evbase", env=parent.frame())
  EcdfTol <- as.vector(get("EcdfTol", env=evbase))
  bL <- 0; ##					@ bounds  @
  bU <- 1;
  wL <- 0;
  wU <- 1;  
  bL <- (bL-Bb)/sigb;
  bU <- (bU-Bb)/sigb;
  wL <- (wL-Bw)/sigw;
  wU <- (wU-Bw)/sigw;
  c1 <- cdfbvng(bU,wU,rho);
  c2 <- cdfbvng(bL,wL,rho);
  c3 <- cdfbvng(bL,wU,rho);
  c4 <- cdfbvng(bU,wL,rho);
  res <- c1+c2-c3-c4;
  res <- recode(res,cbind(res<EcdfTol,res>1),rbind(EcdfTol,1));
  return(res);
}

##
##  a = cdfbvng(dh,dk,r);
##
##  a more accurate version of Gauss' internal CDFBVN().
##  a = double integral from -infty to dk and -infty to dk for a bivariate
##  standard normal integral with correlation r
##
cdfbvng <- function(dh,dk,r){
  dh <- as.matrix(dh)
  dk <- as.matrix(dk)
  tmp <- matrix(0,nrow=rows(dh),ncol=ncol(dh));
  for (i in 1:rows(dh)){
    for(j in 1:ncol(dh))
      tmp[i,j] <- bvnu(dh[i,j],dk[i,j],r)+bvnuphi(dh[i,j])+bvnuphi(dk[i,j])-1;
  }
 
  return(tmp);
}

##
## a = bvnu(dh,dk,r);
##
## a = the probability that x > dh and y > dk under the standard 
##     bivariate normal with correlation r.
##
## To compute the probability that x < dh and y < dk, use 
## bvnu(dh,dk,r) + bvnuphi(dh) + bvnuphi(dk) - 1. 
##
## INPUTS:
##   dh 1st lower integration limit
##   dk 2nd lower integration limit
##   r   correlation coefficient
##
##  by   Alan Genz
##       Department of Mathematics
##       Washington State University
##       Pullman, Wa 99164-3113
##       Email : alangenz@wsu.edu
##
##  based on the method in
##     Drezner, Z and G.O. Wesolowsky, (1989),
##    On the computation of the bivariate normal integral,
##     Journal of Statist. Comput. Simul. 35, pp. 101-107,
##  with major modifications for double precision, for |r| close to 1,
##  and for matlab by Alan Genz - last modifications 7/98. Maximum
##  absolute error for bvnu is approximately 5 x 10^(-16).  Ported
##  to Gauss by Micah Altman and Gary King.

bvnu <- function( dh, dk, r ){
  zero <- 0 
  twopi <- 6.283185307179586 
  w <- matrix(0,nrow=10,ncol=3)
  x <- w
  w[1,1] <- 0.1713244923791705
  w[2,1] <- 0.3607615730481384 
  w[3,1] <- 0.4679139345726904
  x[1,1] <- -0.9324695142031522;
  x[2,1] <- -0.6612093864662647;
  x[3,1] <- -0.2386191860831970;
  w[1,2] <- 0.04717533638651177;
  w[2,2] <- 0.1069393259953183 ;
  w[3,2] <- 0.1600783285433464;
  x[1,2] <- -0.9815606342467191;
  x[2,2] <- -0.9041172563704750 ;
  x[3,2] <- -0.7699026741943050;
  w[4,2] <- 0.2031674267230659;
  w[5,2] <- 0.2334925365383547  ;
  w[6,2] <- 0.2491470458134029;
  x[4,2] <- -0.5873179542866171;
  x[5,2] <- -0.3678314989981802;
  x[6,2] <- -0.1252334085114692;
 ### /* Gauss Legendre points and weights, n = 20*/
  w[1,3] <- 0.01761400713915212;
  w[2,3] <- 0.04060142980038694 ;
  w[3,3] <- 0.06267204833410906;
  w[4,3] <- 0.08327674157670475;
  w[5,3] <- 0.1019301198172404 ;
  w[6,3] <- 0.1181945319615184 ;
  w[7,3] <- 0.1316886384491766;
  w[8,3] <- 0.1420961093183821 ;
  w[9,3] <- 0.1491729864726037 ;
  w[10,3] <- 0.1527533871307259;
  x[1,3] <- -0.9931285991850949;
  x[2,3] <- -0.9639719272779138 ;
  x[3,3] <- -0.9122344282513259;
  x[4,3] <- -0.8391169718222188;
  x[5,3] <- -0.7463319064601508 ;
  x[6,3] <- -0.6360536807265150;
  x[7,3] <- -0.5108670019508271;
  x[8,3] <- -0.3737060887154196;
  x[9,3] <- -0.2277858511416451;
  x[10,3] <- -0.07652652113349733;
  if ( abs(r) < 0.3 ){
    ng <- 1
    lg <- 3
  }else if( abs(r) < 0.75 ){
    ng <- 2
    lg <- 6
  }else{
    ng <- 3
    lg <- 10
  }
  h <- dh;
  k <- dk;
  hk <- h*k;
  bvn <- 0;
  if ( abs(r) < 0.925 ) {
    hs <- ( h*h + k*k )/2
    asr <- asin(r)
    for(i in 1:lg){
      for (i2  in c(-1 , 1)){
        sn <- sin( asr*(  i2*x[i,ng] + 1 )/2) 
        bvn <- bvn + w[i,ng]*exp( ( sn*hk - hs )/( 1 - sn*sn ) )
      }
    }
  
    bvn <- bvn*asr/( 2*twopi ) + bvnuphi(-h)*bvnuphi(-k) 
  }else{
    if ( r < 0 ){
      k <- -k 
      hk <- -hk
    
    }
    if( abs(r) < 1 ){ 
      as <- ( 1 - r )*( 1 + r )
      a <- sqrt(as)
      bs <- ( h - k )^2
      c <- ( 4 - hk )/8
      d <- ( 12 - hk )/16
      asr <- -( bs/as + hk )/2
      if ( asr > -100 )
        bvn <- a*exp(asr)*( 1 - c*(bs-as)*(1-d*bs/5)/3 + c*d*as*as/5 )
     
      if ( -hk < 100 ) {
        b <- sqrt(bs)
        sp <- sqrt(twopi)*bvnuphi(-b/a)
        bvn <- bvn - exp(-hk/2)*sp*b*( 1 - c*bs*(1 - d*bs/5 )/3 )
      }
      a <- a/2;
      for(i in 1:lg){
        for (i2 in c(-1 , 1)){
          xs <- ( a*(  i2*x[i,ng] + 1 ) )^2
          rs <- sqrt( 1 - xs )
          asr <- -( bs/xs + hk )/2;
          if ( asr > -100 ){
            sp <- ( 1 + c*xs*( 1 + d*xs ) )
            ep <- exp( -hk*( 1 - rs )/( 2*( 1 + rs ) ) )/rs
            bvn <- bvn + a*w[i,ng]*exp(asr)*( ep - sp )
          }
        }
      }
      bvn <- -bvn/twopi;
    }
    if (r > 0)  
      bvn <- bvn + bvnuphi( (-1* maxr( h, k )) ); 
    
    if (r < 0) 
      bvn <- -bvn + maxr( zero, bvnuphi(-h) - bvnuphi(-k) ); 
    
  }
  p <- bvn;
  return(p)
}
##  pr = bvnuphi(z)
##  More accurate version of Gauss' cdfn().
##  Normal distribution probabilities accurate to 1.e-15.
##  z = no. of standard deviations from the mean.
## 
##  based upon algorithm 5666 for the error function, from:
##  Hart, J.F. et al, 'Computer Approximations', Wiley 1968
## 
##  programmer: Alan Miller
##  latest revision - 30 march 1986
##  modified for matlab by Alan Genz - 7/98
##  and for Gauss by Micah Altman  
##     
 bvnuphi <- function(z){
 
  p0 <- 220.2068679123761; p1 <- 221.2135961699311; p2 <- 112.0792914978709;
  p3 <- 33.91286607838300; p4 <- 6.373962203531650; p5 <- 0.7003830644436881; 
  p6 <- 0.03526249659989109;
  q0 <- 440.4137358247522; q1 <- 793.8265125199484; q2 <- 637.3336333788311;
  q3 <- 296.5642487796737; q4 <- 86.78073220294608; q5 <- 16.06417757920695; 
  q6 <- 1.755667163182642; q7 <- .08838834764831844;
  rootpi <- 2.506628274631001; cutoff <- 7.071067811865475;
  zabs <- abs(z);
  if ( zabs > 37 )
    p <- 0
  else{
    expntl <- exp(-zabs^2/2);
   ## /*   |z| < cutoff = 10/sqrt(2);
      
  if ( zabs < cutoff ) {
      pt <- ((((((p6*zabs+p5)*zabs+p4)*zabs+p3)*zabs+p2)*zabs+p1)*zabs+p0)
      qt <- ((((((q7*zabs+q6)*zabs+q5)*zabs+q4)*zabs+q3)*zabs+q2)*zabs+q1)
      p <- expntl*pt/( qt*zabs + q0 )
  ##    /*  |z| >= cutoff
    }else{
      p <- expntl/(zabs+1/(zabs+2/(zabs+3/(zabs+4/(zabs+0.65)))))/rootpi;
    }
  }
  
  if (z > 0)
    p <- 1 - p; 
  
  return(p)
}
##/* ----------------------------------------------------------------
## cdfbvnormi(Bb,Bw,sigb,sigw,rho)
## cdf of the bivariate truncated normal with bounds at 0,1, 0,1
## using gauss's CDFBVN -- fast, but imprecise -- with impermissible values
## truncated
##
## Bb = px1 vector of means of var 1
## Bw = means of var 2
## sigb = stan dev of var 1
## sigw = stan dev of var 2
## rho = correlation
##
cdfbvnormi <- function(Bb,Bw,sigb,sigw,rho){
  evbase <- get("evbase", env=parent.frame())
  EcdfTol <- get("EcdfTol", env=evbase)
  bL <- 0; ##					@ bounds  @
  bU <- 1;
  wL <- 0;
  wU <- 1;
  bL <- (bL-Bb)/sigb;
  bU <- (bU-Bb)/sigb;
  wL <- (wL-Bw)/sigw;
  wU <- (wU-Bw)/sigw;
  c1 <- cdfbvn(bU,wU,rho);
  c2 <- cdfbvn(bL,wL,rho);
  c3 <- cdfbvn(bL,wU,rho);
  c4 <- cdfbvn(bU,wL,rho);
  res <- c1+c2-c3-c4;
  res <- recode(res,cbind(res<EcdfTol,res>1),rbind(EcdfTol, 1));
  return(res);
}
##/* ----------------------------------------------------------------
## cdfbvnorme(Bb,Bw,sigb,sigw,rho)
## cdf of the bivariate truncated normal with bounds at 0,1, 0,1
## using Martin van der Ende's bvn_mask()
##
## Bb = px1 vector of means of var 1
## Bw = means of var 2
## sigb = stan dev of var 1
## sigw = stan dev of var 2
## rho = correlation
##
cdfbvnorme <- function(Bb,Bw,sigb,sigw,rho,evbase){
 
  bL <- 0;##					@ bounds  @
  bU <- 1;
  wL <- 0;
  wU <- 1;
  bL <- (bL-Bb)/sigb;
  bU <- (bU-Bb)/sigb;
  wL <- (wL-Bw)/sigw;
  wU <- (wU-Bw)/sigw;
  c1 <- bvn.mask(bU,wU,rho,evbase);
  c2 <- bvn.mask(bL,wL,rho,evbase);
  c3 <- bvn.mask(bL,wU,rho,evbase);
  c4 <- bvn.mask(bU,wL,rho,evbase);
  res <- c1+c2-c3-c4;
  EcdfTol <- as.vector(get("EcdfTol", env=evbase))
  res <- recode(res,cbind(res<EcdfTol,res>1),rbind(EcdfTol,1));
  return(res)
}  


##/* *************************************************************
##   *                                                           *
##   * The procedure bvn_mask calculates the cumulative Normal   *
##   * distribution function cdfbvn(h,k,rho) in a different way. *
##   *                                                           *
##   * The proc bvn_pi2 is a help procedure based on D.R. Divgi: *
##   * "Calculation of the univariate and bivariate Normal       *
##   * integral", The Annals of Statistics, 1979, 903-910.       *
##   *                                                           *
##   * Additions are:                                            *
##   * - some simplifications,                                   *
##   * - more approximating coefficients                         *
##   * - an asymptotical expansion for small probabilities.      *
##   *                                                           *
##   * Accuracy should be over 10 digits absolute and 3-4 digits *
##   * relative. Speed is similar to that of the current Gauss   *
##   * cdfbvn.                                                   *
##   *                                                           *
##   * These procedures may be freely distributed, but not sold. *
##   * I provide it hoping that it may be useful. But I disclaim *
##   * any reliability for possible errors.                      *
##   *                                                           *
##   * Martin van der Ende                                       *
##   *                                                           *
##   * May 28, 1996.                                             *
##   *                                                           *
##   * University of Amsterdam                                   *
##   * vakgroep algemene economie                                *
##   * Roeterssstraat 11                                         *
##   * 1018 WB Amsterdam                                         *
##   * The Netherlands                                           *
##   * email: mvdende@butler.fee.uva.nl                          *
##   * fax: +31.20.5254254                                       *
##   *                                                           *
##   ************************************************************* */
##/* ------------------------------------------------------------------
##** bvn_mask(a,b,c), same inputs as Gauss' cdfbvn, using Divgi's method
##** as extended and programmed by Martin van der Ende
##**
##** By changing comments in bvn_pi2, you can choose 10,20,30 or 40 digits,
##** asymptotic expansions with 1,3,5,7 or 9 terms, and angle reduction to
##** pi/2 (proc bvn_pi2) or pi/4 (proc bvn_pi4). By default, bvn_mask() 
##** uses 30 digits, 9 terms asymptotic expansion, and bvn_pi2.
##*/
bvn.mask <- function(h,k,rho,evbase){
##  /* calculates cdfbvn in a different way */
  
  rho <- rho*(rho<1) + (1-1E-15)*(rho>=1);
  rho <- rho*(rho>-1) - (1-1E-15)*(rho<=-1);
  rho2 <- sqrt(1-rho^2);
  h2 <- (h-rho*k)/rho2;##      /* Pr(x=h conditional on y=k)=Pr(h2) */
  k2 <- (k-rho*h)/rho2;##      /* Pr(y=k conditional on x=h)=Pr(k2) */
  rSq <- h^2 + k2^2;##           /* square of distance of wedge vortex to 0 */
  rSq <- rSq*(rSq<=1400) + 1400*(rSq>1400); ###/* to prevent underflow message */
  rDens <- exp(-rSq/2)/(2*pi); ##/* rho2 times pdfn(h,k,rho) */
  return(bvn.pi2(-h,-k,rho,-h2,-k2,rho2,rSq,rDens))
}

##/*
##** a support proc for bvn_mask, written by Martin van der Ende
##*/
 bvn.pi2 <- function(h,k,rho,h2,k2,rho2,rSq,rDens){
##/*    bvn_pi2 returns cdfbvn(-h,-k,rho).
##      cdfbvn(h,k,rho) returns normpr2(-h,-k,rho,-h2,-k2,rho2,rSq,rDens).
##      Definitions:
##      rho2   = sqrt(1-rho^2);
##      h2     = (h-rho.*k)./rho2;
##      k2     = (k-rho.*h)./rho2;
##      rSq    = h^2 + k2^2;
##      rDens  = exp(-rSq/2)/(2*pi);
##      NOTE: h=rSina1,k2=rCosa1, k=rSina2,h2=rCosa2.
##      NOTE: k2>=0 implies one can take -pi/2<=angle1<=pi/2
##            h2>=0 implies one can take -pi/2<=angle2<=pi/2 */




##   /* SPECIAL (SCALAR) CASES: */
##/* if     maxc(maxc(abs(rho))) > 1;
##          print "abs(correlation coefficient)>1"; Shift=-1;
##   elseif rho==0;  Shift= Uph.*Upk;
##   elseif rho==1;  Shift= Uph.*(h.>k) + Upk.*(h.<=k);
##   elseif rho==-1; Shift= Uph+Upk-1; Shift= Shift.*(Shift.>0);
##   else; */
##     /* GENERAL CASES: */

 polyord <- 31;##  /* order of approximating polynomial plus 1 = 11,21,31,41 */

if (all(polyord==11)){
 vweight <- matrix(0,nrow=21,ncol=1)
 vweight[1] <- 1.253298042;   vweight[2]<- -0.9997316607;
 vweight[3]<- 0.6250192459;   vweight[4]<- -0.3281915667;
 vweight[5]<- 0.1470331965;   vweight[6]<- -5.494856177E-2;
 vweight[7]<- 1.629827794E-2; vweight[8]<- -3.591257830E-3;
 vweight[9]<- 5.406619903E-4; vweight[10]<- -4.890254061E-5;
 vweight[11]<- 1.984741031E-6
}else if (all(polyord==21)){
 vweight<- matrix(0,nrow=21,ncol=1)
 vweight[1]<-  1.253314137222678;     vweight[2]<-  -0.999999996322367;
 vweight[3]<-  0.6266570144788616;    vweight[4]<-  -0.333332913298024;
 vweight[5]<-  0.1566622679316691;    vweight[6]<-  -6.666028411513941E-2;
 vweight[7]<-  2.609623581290622E-2;  vweight[8]<-  -9.499531778787081E-3;
 vweight[9]<-  3.232824755470388E-3;  vweight[10]<- -1.027330989192324E-3;
 vweight[11]<- 3.020025482658464E-4;  vweight[12]<- -8.068099400127647E-5;
 vweight[13]<- 1.912225364901334E-5;  vweight[14]<- -3.912070127260670E-6;
 vweight[15]<- 6.711727217209611E-7;  vweight[16]<- -9.365778519087949E-8;
 vweight[17]<- 1.026548544344371E-8;  vweight[18]<- -8.449983408562188E-10;
 vweight[19]<- 4.887156068638151E-11; vweight[20]<- -1.7641960413010695E-12;
 vweight[21]<- 2.982043919845545E-14
}else if(all(polyord==31)){
 vweight<- matrix(0,nrow=31,1)
 vweight[1]<-  1.253314137315500;     vweight[2]<-  -0.9999999999999754;
 vweight[3]<-  0.6266570686571231;    vweight[4]<-  -0.3333333333248598;
 vweight[5]<-  0.1566642670935284;    vweight[6]<-  -6.666666626444445E-2;
 vweight[7]<-  2.611070955316157E-2;  vweight[8]<-  -9.523804503777532E-3;
 vweight[9]<-  3.263827020322102E-3;  vweight[10]<- -1.058178793302071E-3;
 vweight[11]<- 3.263502098286853E-4;  vweight[12]<- -9.615835274820773E-5;
 vweight[13]<- 2.715573937406461E-5;  vweight[14]<- -7.363023032007412E-6;
 vweight[15]<- 1.915812162821294E-6;  vweight[16]<- -4.766013394595678E-7;
 vweight[17]<- 1.125080846091322E-7;  vweight[18]<- -2.491722870193576E-8;
 vweight[19]<- 5.103898425245527E-9;  vweight[20]<- -9.517178862491669E-10;
 vweight[21]<- 1.589460355702903E-10; vweight[22]<- -2.339412414181092E-11;
 vweight[23]<- 2.985763266418183E-12; vweight[24]<- -3.249122467023523E-13;
 vweight[25]<- 2.958566939223091E-14; vweight[26]<- -2.204022324089611E-15;
 vweight[27]<- 1.304589139468861E-16; vweight[28]<- -5.887558091126596E-18;
 vweight[29]<- 1.899221719494747E-19; vweight[30]<- -3.893896136973921E-21;
 vweight[31]<- 3.807675366775568E-23
}else{
 vweight<- matrix(0,nrow=41,ncol=1)
 vweight[1]<- 1.253314137315500;     vweight[2]<-  -1;
 vweight[3]<- 0.6266570686577501;    vweight[4]<-  -0.3333333333333332;
 vweight[5]<- 0.1566642671644364;    vweight[6]<-  -6.666666666665691E-2;
 vweight[7]<- 2.611071119401255E-2;  vweight[8]<-  -9.523809523527702E-3;
 vweight[9]<- 3.263838898234508E-3;  vweight[10]<- -1.058201055227967E-3;
 vweight[11]<- 3.263838829052882E-4;  vweight[12]<- -9.620008249423181E-5;
 vweight[13]<- 2.719863508376179E-5;  vweight[14]<- -7.399976375870648E-6;
 vweight[15]<- 1.942724550073706E-6;  vweight[16]<- -4.932964393369089E-7;
 vweight[17]<- 1.213895545910440E-7;  vweight[18]<- -2.899419757469565E-8;
 vweight[19]<- 6.728519264552478E-9;  vweight[20]<- -1.517131942551575E-9;
 vweight[21]<- 3.319108847287039E-10; vweight[22]<- -7.022141565916787E-11;
 vweight[23]<- 1.428684402852183E-11; vweight[24]<- -2.773496141377873E-12;
 vweight[25]<- 5.088537684255269E-13; vweight[26]<- -8.730508302470239E-14;
 vweight[27]<- 1.385595302709042E-14; vweight[28]<- -2.012439906068720E-15;
 vweight[29]<- 2.647206396468263E-16; vweight[30]<- -3.122000822598359E-17;
 vweight[31]<- 3.267760772243843E-18; vweight[32]<- -3.003526377626578E-19;
 vweight[33]<- 2.396165106310595E-20; vweight[34]<- -1.637042915419164E-21;
 vweight[35]<- 9.422353450786448E-23; vweight[36]<- -4.474458136625314E-24;
 vweight[37]<- 1.704342515274355E-25; vweight[38]<- -4.999746081902462E-27;
 vweight[39]<- 1.059282556071995E-28; vweight[40]<- -1.441337407752791E-30;
 vweight[41]<- 9.450854298284086E-33
}



 ###    /* ANGLE REDUCTION TO (-PI/2, PI/2]:
###        L(h,k;rho)
###      = cdfnc(h) - L(h,-k;-rho)                  if ((k2.<0) .and (h2.>=0))
###      = cdfnc(k) - L(-h,k;-rho)                  if ((h2.<0) .and (k2.>=0))
###      = -1 + cdfnc(k) + cdfnc(h) + L(-h,-k; rho) if ((h2.<0) .and (k2.<0)). */

 ###    /* APPLY ANGLE REDUCTION TO (-PI/2, PI/2] IF: */
     jkodd <- h2<0;
     jkeven <- k2<0;

###     /* ADD 'Shift' AFTER APPLYING ANGLE REDUCTION */
     Uph <- 1 - pnorm(h) + 1E-300;##    /* +1E-300 to prevent underflow message */
     Shift <- Uph*jkeven;
     Upk <- 1 - pnorm(k) + 1E-300;###    /* +1E-300 to prevent underflow message */
     Shift <- Shift + Upk*jkodd;
     Shift <- Shift - (jkodd & jkeven);

###     /* GET COSINUS AND MINUS SINUS OF REDUCED ANGLES */
     h <- h*(1 - 2*jkeven);
     k <- k*(1 - 2*jkodd);
     h2 <- abs(h2);
     k2 <- abs(k2);

  ###   /* ADJUST CORRELATION */
     corr0 <- abs(rho)<1E-5;

   ###  /* GET REDUCED ANGLE */
     angle <- -atan(rho2/(rho+corr0));
###/*   angle <- -atan2(rho2,rho); */
     angle <- angle + pi*( ((h>0) & (k>0) & (angle<0)) -
                           ((h<0) & (k<0) & (angle>0)) );
###/*   NECESSARY IF YOU REDUCED ANGLE1, ANGLE2 TO ONLY PI/2 IN MODULUS */


###     /* INITIALIZE THE POLYNOMIAL APPROXIMATION */
     jkeven <- angle;  ###                                    /* now jkeven=j0 */
     jkodd <- k+h; ###                                    /* now jkodd=j1 */

     rSqTrunc <- rSq<1400;
     h <- h*rSqTrunc;
     tan2 <- h/(k2+1E-300); ###                  /* = tan(angle2) */
     k <- k*rSqTrunc;
     mtan1 <- k/(h2+1E-300);  ###                /* = -tan(angle1) */
###/*     rho=0; rho2=0;  clear */


###     /* GET THE POLYNOMIAL APPROXIMATION */
     wx <- jkeven - vweight[1]*jkodd;
     rSqMax <- 70;
##/*   rSqTrunc = rSq.<=rSqMax;
##     h        = h.*rSqTrunc; k=k.*rSqTrunc; */
     rSqTrunc <- rSq*rSqTrunc;
     i <- 3;
     ind <- seq(from=3, to=polyord, by=2)
     for(n in ind){
        h <- h*k2; k=k*h2;
        jkeven <- (h+k + (i-2)*rSqTrunc*jkeven)/(i-1);
        h <- h*k2; k=k*h2;
        jkodd <- (h+k + (i-1)*rSqTrunc*jkodd)/i;
        wx <- wx-vweight[i]*jkodd-vweight[i-1]*jkeven;
      }
       

###/*     h=0; k=0; h2=0; k2=0; rSqTrunc=0;  clear */

###     /* APPLY AN ASYMPTOTIC EXPANSION IF NEEDED:
###          rDens.*(a2-a1 - int_a1^a2 1 - 1/(Rcos a)^2 + 3/(Rcos a)^4 -
###                                       15/(Rcos a)^6 da
###        = rDens.*((tan a)./R^2 - (3./R^4).*[(1/3)*(tan a)^3 + tan(a)]_a1^a2 +
###                  (15./R^6).*[(1/5)*(tan a)^5 + (2/3)*(tan a)^3 + tan(a)]_a1^a2
###        = rDens.*((tan a2 - tan a1).*(1./R^2 -3./R^4 + 15./R^6)  -
###                  (tan^3 a2 - tan^3 a1).*(1./R^4 - 10./R^6) +
###                  (tan^5 a2 - tan^5 a1).*(3./R^6) */

     jkodd <- tan2^2; jkeven <- mtan1^2;
     rSq <- rSq+(rSq==0);

 polyord <- 9; ###/* 1,3,5,7 or 9 */

if (polyord==1){
###/*   ONE-TERM ASYMPTOTIC EXPANSION */
     park <- (tan2+mtan1)/rSq;
}else if( polyord==3){
###/*   THREE-TERM ASYMPTOTIC EXPANSION */
     park <- rSq^2 - 3*rSq + 15;
     park <- (tan2+mtan1)*park;
     tan2 <- tan2*jkodd; mtan1=mtan1*jkeven;
     park <- park - (tan2+mtan1)*(rSq-10);
     tan2 <- tan2*jkodd; mtan1=mtan1*jkeven;
     park <- park + (tan2+mtan1)*3;
     park <- park/rSq^3
}else if (all(polyord==5)){
###/*   FIVE-TERM ASYMPTOTIC EXPANSION */
     park <- 945 - rSq*(105-rSq*(15-rSq*(3-rSq)))
     park <- (tan2+mtan1)*park
     tan2 <- tan2*jkodd; mtan1=mtan1*jkeven
     park <- park + (tan2+mtan1)*(1260 - rSq*(105-rSq*(10-rSq)))
     tan2 <- tan2*jkodd; mtan1 <- mtan1*jkeven
     park <- park + 3*(tan2+mtan1)*(378-rSq*(21-rSq))
     tan2 <- tan2*jkodd; mtan1=mtan1*jkeven
     park <- park + 15*(tan2+mtan1)*(36-rSq)
     tan2 <- tan2*jkodd; mtan1=mtan1*jkeven
     park <- park + (tan2+mtan1)*105
     park <- park/rSq^5
}else if( polyord==7){
###/*   SEVEN-TERM ASYMPTOTIC EXPANSION */
     park <- 135135 - rSq*(10395-rSq*(945-rSq*(105-rSq*(15-rSq*(3-rSq)))))
     park <- (tan2+mtan1)*park
     tan2 <- tan2*jkodd; mtan1 <- mtan1*jkeven
     park <- park + (tan2+mtan1)*(270270 - rSq*(17325-rSq*(1260-rSq*(105-rSq*(10-rSq)))))
     tan2 <- 3*tan2*jkodd; mtan1 <- 3*mtan1*jkeven
     park <- park + (tan2+mtan1)*(135135 - rSq*(6930-rSq*(378-rSq*(21-rSq))))
     tan2 <- 5*tan2*jkodd; mtan1 <- 5*mtan1*jkeven
     park <- park + (tan2+mtan1)*(25740 - rSq*(990-rSq*(36-rSq)));
     tan2 <- 7*tan2*jkodd; mtan1 <- 7*mtan1*jkeven;
     park <- park + (tan2+mtan1)*(2145-rSq*(55-rSq))
     tan2 <- 9*tan2*jkodd; mtan1 <- 9*mtan1*jkeven
     park <- park + (tan2+mtan1)*(rSq-78)
     tan2 <- 11*tan2*jkodd; mtan1 <- 11*mtan1*jkeven
     park <- (park+tan2+mtan1)/rSq^7
}else{
###/*   NINE-TERM ASYMPTOTIC EXPANSION */
     park <- 34459425- rSq*(20270-rSq*(135135-rSq*(10395-rSq*(945-rSq*(105-rSq*(15-rSq*(3-rSq)))))))
     park <- (tan2+mtan1)*park
     tan2 <- tan2*jkodd; mtan1 <- mtan1*jkeven
     park <- park + (tan2+mtan1)*(91891800 -rSq*(4729725-rSq*(270270-rSq*(17325-rSq*(1260-rSq*(105-rSq*(10-rSq)))))))
     tan2 <- 3*tan2*jkodd; mtan1 <- 3*mtan1*jkeven
     park <- park + (tan2+mtan1)*(64324260 -rSq*(2837835-rSq*(135135-rSq*(6930-rSq*(378-rSq*(21-rSq))))))
     tan2 <- 5*tan2*jkodd; mtan1 <- 5*mtan1*jkeven
     park <- park + (tan2+mtan1)*(18378360 -rSq*(675675-rSq*(25740-rSq*(990-rSq*(36-rSq)))))
     tan2 <- 7*tan2*jkodd; mtan1 <- 7*mtan1*jkeven
     park <- park + (tan2+mtan1)*(2552550 -rSq*(75075-rSq*(2145-rSq*(55-rSq))))
     tan2 <- 9*tan2*jkodd; mtan1 <- 9*mtan1*jkeven
     park <- park + (tan2+mtan1)*(185640 -rSq*(4095-rSq*(78-rSq)))
     tan2 <- 11*tan2*jkodd; mtan1 <- 11*mtan1*jkeven
     park <- park + (tan2+mtan1)*(7140 -rSq*(105-rSq))
     tan2 <- 13*tan2*jkodd; mtan1 <- 13*mtan1*jkeven
     park <- park + (tan2+mtan1)*(136-rSq)
     tan2 <- 15*tan2*jkodd; mtan1 <- 15*mtan1*jkeven
     park <- (park+tan2+mtan1)/rSq^9
   }

     wx <- wx - (wx-park)*(rSq>rSqMax);

   ###  /* GET THE RETURN VALUE FOR NONZERO CORRELATION */
     wx <- Shift + rDens*wx;

   ###  /* REPLACE THE RETURN VALUE FOR ZERO CORRELATION */
     wx <- wx + (Uph*Upk-wx)*corr0;

###/* endif; */
   return(wx)
}
##/* ----------------------------------------------------------------
##** cdfbvnorma(Bb,Bw,sigb,sigw,rho)
##** cdf of the bivariate truncated normal with bounds at 0,1, 0,1
##** with home grown numerical approximation. slower, but usually more precise,
##** for very small values, this function is not smooth.
##**
##** Bb = px1 vector of means of var 1
##** Bw = means of var 2
##** sigb = stan dev of var 1
##** sigw = stan dev of var 2
##** rho = correlation
##**
##*/
cdfbvnorma <- function(mu1,mu2,s1,s2,rho){

  l1 <- 0;##					@ bounds  @
  u1 <- 1;
  l2 <- 0;
  u2 <- 1;

  k <- 40;##					@ num of numerical intervals @
  r <- nrow(as.matrix(mu1));
  x <- seqase(l1,u1,k);
  
  s1s <- s1^2;
  omr <- (1-rho^2)*s2^2;
  
  b <- a <- matrix(0,nrow=k,ncol=r);
  
  for( i in 1:k){
    a[i+0,] <- as.vector(pdfnorm(x[i+0],mu1,s1s))
    cdf1 <- mu2+(rho*(s2/s1)*(x[i+0]-mu1)) 
    b[i+0,] <- cdfnormdif(cdf1,omr,l2,u2)
  }
  
  aa <- cdfnormdif(mu1,s1s,l1,u1);
  
  res <- aa*(colSums(a*b)/colSums(a))
  return(res)
}
##/************
##** PDF'S ****
##*************
##*/
##/*----------------------------------------------------------------
##** pdf of a normal distribution (not standardized)
##**
##** usage:  p = pdfnorm(y,mu,sigma2);
##**
##** INPUTS:
##** mu = mean
##** sigma2 = variance
##** y = normal variate
##**
##** OUTPUT:
##** p value of the probability density
##**
##*/
pdfNorm <- pdfnorm <- function(y,mu,s2){
  
  s <- sqrt(s2);
  res <- dnorm((y-mu)/s)/s;
  return(res);
}

###/* support proc for cdfbvnorma */  
cdfnormdif <- function(mu,s2,l,u){
  
  a <- cdfnorm(u,mu,s2)-cdfnorm(l,mu,s2);
  return(a)
}

    
##/* ----------------------------------------------------------------
##**    y = cdfbvnunit(bb,bw,sb,sw,rho);
##**
##** integral of the bivariate normal distribution above the unit square
##** slower and fairly precise.  Function is not smooth for very small values.
##**
##** INPUTS:  (of the bivariate normal)
##** bb,bw = means
##** sb,sw = standard deviations
##** rho   = correlation
##**
##** OUTPUT
##** y = probability of being in the unit square
##*/
cdfbvnunit <- function(bb,bw,sb,sw,rho,evbase=NULL){
  if(!length(evbase))
    evbase <- get("evbase", env=parent.frame())
  
   cdfbvnunit.paramsB <- cbind(bb,bw);
   cdfbvnunit.paramsS <- rbind(sb,sw,rho);
   intord <- 40;
   assign("cdfbvnunit.paramsB", cdfbvnunit.paramsB, env=evbase)
   assign("cdfbvnunit.paramsS", cdfbvnunit.paramsS, env=evbase)
   assign("intord",intord, env=evbase)
  lst <- integrate(cdfbvnunit.proc,lower=0, upper=1, subdivisions=100,evbase)
  res <- lst$value
  EcdfTol <- as.vector(get("EcdfTol", env=evbase))
  res <- recode(res,cbind((res<EcdfTol),(res>1)),rbind(EcdfTol,1));
  return(res);
}
  
##/*
##**  y = _cdfbvnunit_proc(betaw);
##**
##** support proc for cdfbvnunit(), written so that only one integral is
##** necessary to get double integral on unit square
##**
##** y = value of the density to be integrated for:
##**
##** betab,betaw = bivariate normal variates, with:
##** bb,bw       = means
##** sb,sw       = standard deviations
##** rho         = correlation
##*/
cdfbvnunit.proc <- function(betaw,evbase){
##  evbase <- get("evbase", env=parent.frame())
  cdfbvnunit.paramsB <- get("cdfbvnunit.paramsB", env=evbase)
  cdfbvnunit.paramsS <- get("cdfbvnunit.paramsS", env=evbase)
  bb <- cdfbvnunit.paramsB[,1];
  bw <- cdfbvnunit.paramsB[,2];
  sb <- cdfbvnunit.paramsS[1];
  sw <- cdfbvnunit.paramsS[2];
  rho <- cdfbvnunit.paramsS[3];
  sb2 <- sb^2;
  sw2 <- sw^2; 
  sbw <- rho*sb*sw;
  sbw2 <- sbw^2;
  mu <- bb+(sbw/sw2)*(betaw-bw);
  var <- sb2-sbw2/sw2;
  
  res <- pdfnorm(betaw,bw,sw2);
  res <- res*(cdfnorm(1,mu,var)-cdfnorm(0,mu,var));

  return(res);
}
##/* ----------------------------------------------------------------
##**  y = lcdfbvnorma(Bb,Bw,sigb,sigw,rho)
##**    
##** ln of cdf of the bivariate truncated normal with bounds at 0,1, 0,1
##** with home grown numerical approximation, EXPERIMENTAL
##**
##** Bb = px1 vector of means of var 1
##** Bw = means of var 2
##** sigb = stan dev of var 1
##** sigw = stan dev of var 2
##** rho = correlation
##**
##*/
lcdfbvnorma <- function(mu1,mu2,s1,s2,rho){
  mu1 <- as.matrix(mu1)
  l1 <- 0; ##					@ bounds  @
  u1 <- 1;
  l2 <- 0;
  u2 <- 1;

  k <- 20;	###			@ num of numerical intervals @
  r <- nrow(mu1);
  s1s <- s1^2;
  omr <- (1-rho^2)*s2^2;
  aa <- lcdfnormi(cbind(l1,u1),mu1,s1s);
  b <- a <- matrix(0,nrow=k,ncol=r);
 
  for (i in 1:k){
    t <- pnorm((l1-mu1)/s1)+(((2*i-1)/(2*k))*exp(aa))
    x <- mu1+s1*qnorm(t)
    cdf1 <- mu2+(rho*(s2/s1)*(x-mu1)) 
    b[i+0,] <- as.vector(lcdfnormi(cbind(l2,u2),cdf1,omr));
  }

  res <- aa+log(colMeans(exp(b)))
  return(res);
}

##/* ----------------------------------------------------------------
##** lncdfbvn2i(Bb,Bw,sigb,sigw,rho)
##** cdf of the bivariate truncated normal with bounds at 0,1, 0,1
##** using gauss's new lncdfbvn2 proc
##**
##** Bb = px1 vector of means of var 1
##** Bw = means of var 2
##** sigb = stan dev of var 1
##** sigw = stan dev of var 2
##** rho = correlation
##*/
lncdfbvn2i <- function(Bb,Bw,sigb,sigw,rho){
 
  evbase <- get("evbase", env=parent.frame())
  EcdfTol <- get("EcdfTol", env=evbase)
  res <- cdfbvn2(-Bb/sigb,1/sigb,-Bw/sigw,1/sigw,rho);
 
  res <- recode(res,cbind((res<EcdfTol),(res>1)),rbind(EcdfTol,1));
 
  return(log(res))
}


##/*
##** g = lpdfbvn(betab,betaw,bb,bw,sb,sw,rho);
##**                                        
##** ln of the pdf of the bivariate normal distribution
##*/
lpdfbvn <- function(betab,betaw,bb,bw,sb,sw,rho){
  evbase <- get("evbase", env=parent.frame())
  Evtol <- EvTol <- get("EvTol", env=evbase)
 
  z1 <- (betab-bb)/sb;
  z2 <- (betaw-bw)/sw;
  omr2 <- (1-rho^2);
  omr2 <- recode(omr2,omr2<Evtol,Evtol);
  sb <- recode(sb,sb<Evtol,Evtol);
  sw <- recode(sw,sw<Evtol,Evtol);
 
  w <- (z1^2+z2^2-2*rho*z1%dot*%z2)/omr2;
  res <- -1.83787706640934548-0.5*(w+log(omr2))-log(sb)-log(sw);
  return(res)
}

##/* ----------------------------------------------------------------
##**  a = lpdfnorm(y,mu,sig2);
##**
##** a more computationally precise method of computing ln(pdfnorm())
##**
##*/
lpdfnorm <- function(y,mu,sig2){
  
  x=y-mu;
  res <- -0.918938533204672741 - ( log(sig2) +  ( (x*x) / sig2 ) ) /2;
  return(res)
}
equalev <- function(vec,x,tol=1.e-4){
  res <- sapply(vec,function(ll){
    res <- TRUE
    err  <- ll-x
    err1 <- abs(err/ll)
    err2 <- abs(err/x)
    err <- ifelse(err1 > err2, err1, err2)
    if(err > tol) res <- FALSE
    return(res)
  })
  res <- unlist(res)
  if(any(res == FALSE)) return(FALSE)
  return(TRUE)
}
    
