##
##  This archive is part of the program EI
##  (C) Copyright 1995-2001 Gary King
##  All Rights Reserved.
##
## log-likelihood function for ecological inference model
##
##
#include ei.ext;
eiloglik <- function(b, dta,evbase=NULL,...){
  if(!length(evbase))
    evbase <- get("evbase", env=parent.frame())
  
evlocal <- getEnvVar(evbase, environment())
##  local sb2,sw2,sbw,x,y,llik,s2,bb,bw,mu,sb,sw,c0,c,c1,cT0,cT1,Zb,Zw,
##     rho,tt,bnds,R,omega,epsilon,Ebb,Vbb,res,prior,rs,z,o;

 lst <- pluckdta(dta,evbase);

 Zb <- lst$Zb
 Zw <- lst$Zw
 x <- lst$x
 y <- lst$y
 rs <- nrow(as.matrix(y));

###  /* reparameterize */

 lst <- eirepar(b,Zb,Zw,x,Ez,evbase=evbase);  ##tested 

### lst <- c(list(Bb=Bb), list(Bw=Bw), list(sb=sb), list(sw=sw), list(rho=rho))
 bb <- lst$Bb
 bw <- lst$Bw
 sb <- lst$sb
 sw <- lst$sw
 rho <- lst$rho
 sb2 <- sb^2;
 sw2 <- sw^2;

 ### /* divide up types of observations */
 llik <- matrix(0,nrow=rs,ncol=1);
 lst <- homoindx(x,get("EnumTol", env=evbase));
 c <- lst$c
 c0 <- lst$c0
 c1 <- lst$c1
 if(scalmiss(c))
   cT1 <- cT0 <- NA
 else{
   cT0 <- subset(c, subset=(y[c]< as.vector(EnumTol)));
   if(!length(cT0)) cT0 <- NA
   cT1 <- subset(c,subset=(y[c]>(1-as.vector(EnumTol))));
   if(!length(cT1)) cT1 <- NA
    c <- subset(c,subset=((y[c]>=as.vector(EnumTol))& (y[c]<=(1-as.vector(EnumTol)))));
    if(!length(c)) c <- NA
  }

 ### /* compute likelihood for different categories */
  if(!scalmiss(c0)){	###		@ X=0 @
        
    epsilon <- y[c0]-bw[c0];
    llik[c0] <- -0.5*(log(sw2)+(epsilon^2)/sw2) ### @ ln N(T|Bw,sigmaW) @
    bnds <- cbind(matrix(0,nrow =rows(c0),ncol=1), matrix(1,nrow=rows(c0),ncol=1));
    Ebb <- bb[c0]+rho *(sb/sw)*epsilon;
    Vbb <- sb2*(1-rho^2);
    res <- lcdfnormi(bnds,Ebb,Vbb);      ###     @ ln S'(Bu,Sigmau) @
    R <- lncdfbvnu(bb[c0],bw[c0],sb,sw,rho)###  @ ln R(Bu,Sigmau) @
    llik[c0] <- llik[c0]+res-R
  }

  if(!scalmiss(c1)){ ##			@ X=1 @
  
    epsilon <- y[c1]-bb[c1];
    llik[c1] <- -0.5*(log(sb2)+(epsilon^2)/sb2)### @ ln N(T|Bb,sigmaB) @
    bnds <- cbind(matrix(0,nrow=rows(as.matrix(c1)),ncol=1),matrix(1,nrow=rows(as.matrix(c1)),ncol=1))
    Ebb <- bw[c1]+rho%dot*%(sw%dot/%sb)%dot*%epsilon
    Vbb <- sw2%*%(1-rho^2)
    res <- lcdfnormi(bnds,Ebb,Vbb) ##           @ ln S'(Bu,Sigmau) @
    R <- lncdfbvnu(bb[c1],bw[c1],sb,sw,rho)##  @ ln R(Bu,Sigmau) @
    llik[c1] <- llik[c1]+res-R
  }

  if(!scalmiss(cT0)){###		@ T=0, 0<X<1 @
     
    z <- matrix(0,nrow=rows(as.matrix(cT0)),ncol=1)
    llik[cT0] <- lpdfbvn(z,z,bb[cT0],bw[cT0],sb,sw,rho) -lncdfbvnu(bb[cT0],bw[cT0],sb,sw,rho);
  }

  if(!scalmiss(cT1)){	###	@ T=1, 0<X<1 @
     
    o <- matrix(1,nrow=nrow(as.matrix(cT1)),ncol=1)
    llik[cT1] <- lpdfbvn(o,o,bb[cT1],bw[cT1],sb,sw,rho) -lncdfbvnu(bb[cT1],bw[cT1],sb,sw,rho)
  }
  
  if(!scalmiss(c)){ ###			@ 0<T<1, 0<X<1 @
   
    lst <- exvar(y[c],x[c],bb[c],bw[c],sb,sw,rho)  ##tested
  ###  {mu,s2,epsilon,omega,Ebb,Vbb} =
    
    mu <- lst$mu
    s2 <- lst$s2
    epsilon <- lst$epsilon
    omega <- lst$omega
    Ebb <- lst$Ebb
    Vbb <- lst$Vbb
###     llik[c] <- -0.5*(log(s2)+(epsilon^2)%dot/%s2) ###     @ ln N(T|mu,sigma) @
   
    llik[c] <- -0.5*(log(s2)+(epsilon^2)%dot/%s2) ###     @ ln N(T|mu,sigma) @
    lst <- bounds1(y[c],x[c],matrix(1, nrow=nrow(as.matrix(c)),ncol=1),EnumTol)
    bnds <- lst$bs
    tt <- lst$aggs
 
    res <- lcdfnormi(bnds[,1:2],Ebb,Vbb) ##  einormal.R           @ ln S(Bu,Sigmau) @
   
    R <- lncdfbvnu(bb[c],bw[c],sb,sw,rho)## einormal.R             @ ln R(Bu,Sigmau) @
   
    llik[c] <- llik[c]+res-R
  }

  ###/* priors */
  prior <- 0
  if (Esigma>0)
    prior <- prior-(1/(2*Esigma^2))*(sb2+sw2)	###      @ sb, sw @
  
  if(Erho[1]>0)
    prior <- prior+lpdfnorm(b[nrow(as.matrix(b))- 2],0,Erho[1]^2) ## einormal.R  @ rho @
  
  if(Ebeta>0){			       ### 	      @ bb, bw @
    prior <- prior+flatnorm(colMeans(as.matrix(bb)),Ebeta)
    prior <- prior+flatnorm(colMeans(as.matrix(bw)),Ebeta)
  }
  if(!scalmiss(EalphaB))###                         @ alphaB @
    prior <- prior+colSums(lpdfnorm(b[2:Ez[1]],EalphaB[,1],EalphaB[,2]^2))
  
  if(!scalmiss(EalphaW))###                         @ alphaW @
    prior <- prior +colSums(lpdfnorm(b[(Ez[1]+2):colSums(Ez)],EalphaW[,1],EalphaW[,2]^2))
 
  
  llik <- llik%plus%(prior/rs);
  
  res <- missrv(llik,0)
  
  return(res)

}

###
##  l = homoindx(x);
##  INPUT: is a matrix of one column or array.
##
##  OUTPUT: a list with 3 elements
##          c  = vector of index numbers for heterogeneous precincts
##          c0 = vector of index numbers for homogenous white (x=0) precincts
##          c1 = vector of index numbers for homogenous black (x=1) precincts
##          Translates code in Gauss written by Gary King
## Ferdinand Alhimadi & Elena Villalon (evillalon@iq.harvard.edu)
## Date August 17th, 2007

homoindx<-function(x, EnumTol=0.0001){
  
      ###  EnumTol <- try(get("EnumTol", env=get("evbase", env=parent.frame())))
        if(!length(EnumTol))
          EnumTol <- 0.0001
        x <- as.matrix(x)
        indx <- seq(1,rows(x), 1)
        
        res<-list()
        EnumTol <- as.vector(EnumTol)
        c0<- x< EnumTol
        c1<- x>(1- EnumTol)
        c<-(1-c0-c1)
        cl <- as.logical(c)
        res$c0<-indx[x <EnumTol]
        if(!length(res$c0)) res$c0 <- NA
        res$c1<-indx[x >(1-EnumTol) ]
        if(!length(res$c1)) res$c1 <- NA
        res$c<-indx[c==TRUE]
        if(!length(res$c)) res$c <- NA
        return (res)
}


##/*
##    {mu,s2,epsilon,omega,Ebb,Vbb} = exvar(T,x,bbetaB,bbetaW,sigb,sigw,rho);
##or call as:
##    {mu,s2,epsilon,omega,Ebb,Vbb} = exvar(T,x,eirepar(params,Zb,Zw));
##**
##**  reparameterize & compute expected value and variance
##**  support proc for psim1() and eiloglik()
##**
##** INPUTS:
##** T = dep var
##** x = explanatory variable
##** bbetaB,bbetaW,sigb,sigw,rho = output of eirepar
##**
##** OUTPUTS:
##** mu = E(T|X), all on untruncated scale
##** s2 = V(T|X)
##** epsilon = t-mu
##** omega = cov(t,betaB)
##** Ebb = E(betaB)
##** Vbb = V(betaB)
##**  
##*/
exvar <- function(t,x,bbetaB,bbetaW,sigb,sigw,rho,evbase=NULL){
  if(!length(evbase))
    evbase <- get("evbase", env=parent.frame())
  EvTol <- get("EvTol", env=evbase)
  sigb2 <- sigb^2; ###scalar
  sigw2 <- sigw^2; ####scalar as it is rho
  sigbw <- rho%*%sigb%*%sigw; ##scalar 
  
  omx <- 1-x; ### 1 column and p rows

  mu <- as.matrix(bbetaB)%dot*%as.matrix(x)+as.matrix(bbetaW)%dot*%as.matrix(omx);
  epsilon <- t-mu;
  
  ##here x*omx = px1 and is the same as x%dot*%omx
  s2 <- as.matrix((sigb2*(x^2)))+as.matrix(sigw2*(omx^2))+(2*as.vector(sigbw)*(as.matrix(x)%dot*%as.matrix(omx)));
  
  omega <- sigb2*x+sigbw%dot*%omx;
  
  Ebb <- bbetaB+((omega%dot/%s2)%dot*%epsilon);
  Vbb <- sigb2-((omega^2)%dot/%s2);
 
  EvTol <- as.vector(EvTol)
 
  
  if(any(Vbb < EvTol)){
  
    Vbb[Vbb < EvTol] <- EvTol
  }
 
 ### Vbb <- recode(Vbb,Vbb<EvTol,EvTol); ### @ fix numerical innacuracies @
  lst <- c(list(mu=mu), list(s2=s2), list(epsilon=epsilon), list(omega=omega), list(Ebb=Ebb), list(Vbb=Vbb))
  return(lst)
}



##/*
##**  prior = flatnorm(bb,sig);
##**  support proc for eiloglik()
##**
##** INPUT:
##** bb = parameter to put a flat normal prior on
##** sig = standard deviation of normal
##**
##** OUTPUT: log of flat normal prior
##**
##** the prior is flat in [0,1] and falls off with the normal distribution
##** outside that range
##*/
flatnorm <- function(bb,sig){
  
  sig2 <- sig^2;
  if (bb<0)
    prior <- -0.5*((bb^2)/sig2)
  else if( bb>1)
    prior <- -0.5*(((bb-1)^2)/sig2)
  else
    prior <- 0

  
  return(prior)
}

