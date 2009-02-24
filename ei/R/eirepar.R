##
##  This archive is part of the program EI
##  (C) Copyright 1995-2001 Gary King
##  All Rights Reserved.
##
##  various reparameterization procs
##
##
## {bb,bw,sb,sw,rho} = eirepar(params,Zb,Zw,x);
## 
## INPUTS: 
## params = params in scale of estimation (eta's at end)
## Zb,Zw = covariates (without constant term) or scalar 1
## x = aggregate-level classifying variable
## 
## OUTPUTS: on untruncated normal scale
## bb,bw = Px1 vectors of precinct means
## sb,sw = 1x1 standard deviations
## rho   = correlation 
##
## reparametrization from scale of estimation
## to the untruncated scale.  (see eirepart.g for truncated scale)
##
#include ei.ext;
eirepar <- function(b,Zb,Zw,x, Ez=NULL, evbase=parent.frame()){
 
  
 ## /* pluck off params */
  if(!length(Ez))
    Ez <- get("Ez", evbase)

  o <- matrix(1, nrow=rows(x),ncol=1);

  Bb0 <- b[1:Ez[1]];
  if (Ez[1]==1)
    Bb0v <- 0
  else
    Bb0v <- trimr(Bb0,1,0); ### @ vector of params for mean-adjusted Zb @
 
  Bb0 <- Bb0[1]*o;###		  @ constant term @
  
  Bw0 <- b[(Ez[1]+1):colSums(as.matrix(Ez))];
  if(Ez[2]==1)
    Bw0v <- 0
  else
    Bw0v <- trimr(Bw0,1,0); ### @ vector of params for mean-adjusted Zw @
  
  Bw0 <- Bw0[1]*o;	###	  @ constant term @
    
  r <- rows(b);
  sb0 <- b[r-4];
  sw0 <- b[r-3];
  rho0 <- b[r-2];
  etaB <- b[r-1];
  etaW <- b[r];
   
  ###/* reparameterize */
  sb <- exp(sb0);
  sw <- exp(sw0);

  m <- x-colMeans(as.matrix(x));
  
###  if( any(is.na(eiread.zb)))
###    Zb <- Zb-t(colMeans(as.matrix(Zb)));
 
###  if(any(is.na(eiread.zw)))
###    Zw <- Zw-t(colMeans(as.matrix(Zw)))
 
  
  Zb <- Zb-t(colMeans(as.matrix(Zb)));
  Zw <- Zw-t(colMeans(as.matrix(Zw)));
 
  
  
  Bb <- Bb0*(0.25+sb^2)+0.5+(Zb*Bb0v+etaB*m); ###matrix multp 
  Bw <- Bw0*(0.25+sw^2)+0.5+(Zw*Bw0v+etaW*m); ###matrix multp
  ### I need -rho0 to match Gauss 
  rho <- fisherzi(rho0);
  lst <- c(list(Bb=Bb), list(Bw=Bw), list(sb=sb), list(sw=sw), list(rho=rho))
  return(lst)
}
###/* ----------------------------------------------------------------
###  b = EireparT(params,Zb,Zw,x);
###**
##**  reparameterization to the ultimate Truncated scale of E() and V() of
##**  the precinct parameters  (see eirepar() for the intermediate 
##**  parameterization)
##**
##** params = bb|bw|sb|sw|rho|etaB|etaW (untruncated)
##**  b     = bb|bw|sb|sw|rho (truncated)
##**
##** uses 100*_Esims simulations for computations
##*/
 eirepart <- function(params,Zb,Zw,x, Ez, evbase=try(get("evbase", env=parent.frame())))
{

  Esims <- get("Esims", env=evbase)
  lst <- eirepar(params,Zb,Zw,x,Ez,evbase=evbase)
  bb <- lst[[1]]
  bw <- lst[[2]]
  sb <- lst[[3]]
  sw <- lst[[4]]
  rho <- lst[[5]]
  bb <- colMeans(as.matrix(bb))
  bw <- colMeans(as.matrix(bw))
  sims <- Esims*100
  bounds <- matrix(c(0,1,0,1),nrow=2, ncol=2,byrow=TRUE)
  
  b <- rndbtn(bb,bw,sb,sw,rho,bounds,sims)
  t <- cor(b)
  b <- colMeans(as.matrix(b))
  std <- sd(as.data.frame(b))
  b <- as.matrix(c(b,std,t[2, 1]))
  return(b)
}
