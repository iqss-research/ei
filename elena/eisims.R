##/*
##**  This archive is part of the program EI
##**  (C) Copyright 1995-2001 Gary King
##**  All Rights Reserved.
##*/
##/*
##   {betaBs,betaWs} = psim1(T,X,N,Zb,Zw,MLpsi,MLvc);
##**
##** precinct-level parameter simulations
##**
##** INPUTS:
##** T    = dependent variable (p x 1)
##** X    = explanatory variable (p x 1)
##** N    = denominator of T and X (e.g., total voting age population) (p x 1)
##** Zb   = 1 or explanatory vars to predict beta-b 
##** Zw   = 1 or explanatory vars to predict beta-w 
##** MLpsi,MLvc = output from quadcml()
##**
##** GLOBAL:
##** _EbetaWs = 1 compute betaws as output
##**            0 set to missing (default) (can use eiread to compute)
##**            
##** OUTPUT:
##** betaBs,betaWs  =  (p x sims) matrices of simulations
##**
##*/
#include ei.ext;
psim1 <- function(T,X,tvap,Zb,Zw,MLpsi,MLvc,evbase=get("evbase", env=parent.frame())){
##  local PSIsims,k,betaBs,betaWs,Bbeta,aggs,Ebb,Vbb,tt,dataset,
##  etaB,etaW,etaBs,etaWs,c,c0,c1,i,rs,r;
  Eprt <- get("Eprt", env=evbase)
  Eeta <- get("Eeta", env=evbase)
  Erho <- get("Erho", env=evbase)
  Esims <- get("Esims", env=evbase)
  EisChk <- get("EisChk", env=evbase)
  Ebetaws <- get("Ebetaws", env=evbase)
  ##  /* ESTIMATION VARIABILITY */
  if(Eprt>=2)
    message("Simulating estimation variation...");
 
 ## /* _Eselect is applied within packdta, so dataset is _Eselect'd but
  ##    x,Zb,Zw,t, which are used for fundamental variability, are not */
  dataset <- packdta(x,Zb,Zw,t);
  
 ## /* importance sampling for b's and normal for eta's:  sims x rows(MLpsi) */
  etaB <- MLpsi[rows(MLpsi)-1];
  etaW <- MLpsi[rows(MLpsi)];
  if (rows(Eeta)==1){
    etaBs <- 0
    etaWs <- 0
  }else if(rows(Eeta)==4){
    etaBs <- Eeta[3]
    etaWs <- Eeta[4]
  }else if(Eeta[1]==4){
    etaBs <- 0
    etaWs <- Eeta[3]
  }else if(Eeta[1]==5){
    etaBs <- Eeta[3]
    etaWs <- 0
  }
  if(Erho[1]==0){
    EisFix <- MLpsi[(rows(MLpsi)-2):rows(MLpsi)]
    yy1 <- rndisamp(eiloglik,trimr(MLpsi,0,3), MLvc[1:rows(MLpsi)-3, 1:rows(MLpsi)-3],
                    dataset,Esims) ###COMPUTE
    yy2 <-  rndmn(Erho[2],0,Esims)  ##COMPUTE
    yy3 <-  rndmn(etaB,etaBs,Esims)
    yy4 <- rndmn(etaW,etaBs,Esims)
    PSIsims <- cbind(yy1, yy2,yy3, yy4)
                    
                    
  }else{
    EisFix <- MLpsi[(rows(MLpsi)-1):rows(MLpsi)]
    yy1 <- rndisamp(eiloglik,trimr(MLpsi,0,2),
                    MLvc[1:rows(MLpsi)-2,1:rows(MLpsi)-2],dataset,Esims) ##COMPUTE
    yy2 <- rndmn(etaB,etaBs,Esims)  ##COMPUTE
    yy3 <-   rndmn(etaW,etaBs,Esims)
    PSIsims <- cbind(yy1, yy2, yy3)
                   
                   
  }
  psisims <-  as.matrix(PSIsims)
  if(EisChk)                          
    Eres <- vput(Eres,psisims,"PhiSims") ###  @ save name corresponds to book @ 
  else 
    Eres <- vput(Eres,cbind(colMeans(psisims)),sd(as.data.frame(psisims)),"PhiSims")
  
  
 ### /* FUNDAMENTAL VARIABILITY */
 ### /* Applies to all observations */
  if(Eprt>=2)
    message("Simulating fundamental variability...")
  
  rs <- rows(x);
  betaBs <- matrix(0, nrow=rs,ncol=Esims)
  betaWs <- betaBs
  lst <- bounds1(t,x,tvap) 
  Bbeta <- lst[[1]]
  aggs <- lst[[2]]
  for( k in (1:Esims)){
    y <- eirepar(PSIsims[k+0,],Zb,Zw,x,evbase)
    lst <- exvar(T,x,y)
    tt1 <- lst[[1]]
    tt2 <- lst[[2]]
    tt3 <- lst[[3]]
    tt4<- lst[[4]]
    Ebb <- lst[[5]]
    Vbb <- lst[[6]]
    betaBs[,k+0] <- rndtni(Ebb,Vbb,Bbeta[,1:2]) ###COMPUTE
  }

 ### /* compute betaWs from betaBs deterministically */
  if (Ebetaws)
    betaWs <- betab2w(t,x,betaBs,evbase)
  else
    betaWs <- NA
  
  lst <- c(list(betaBs=betaBs), list(betaWs=betaWs))
  return(lst)
}
##/*
## betaw = betab2w(t,x,betab);
##
## compute betaw deterministically given betab
##
## INPUT DIMENSIONS:
## rows of all inputs must be p
## ncol(t) must equal 1
## ncol(x) must be 1 (for ei) or Esims (for ei2)
## ncol(betab) must be 1 (for mean posterior) or Esims (for betaBs)
##/*

betab2w <- function(t,x,betab, evbase=NULL){
###  local c,c0,c1,betaw,i,col,x1;
  if(!length(evbase))
    evbase <- get("evbase", env=parent.frame())
  col <- cols(x);
  x1 <- 1-x;
  betaB <- betab
  betaW <- betaw <- betab;
  
  if(col==1){
     cs<-homoindx(x)
     c<-cs$c
     c0<-cs$c0
     c1<-cs$c1
    
    if(!scalmiss(c))
      betaW[c,] <- (t[c]/x1[c])-((betaB[c,]*x[c])/x1[c]);
    
    if(!scalmiss(c0))
      betaW[c0,] <- matrix(t[c0],nrow=rows(c0),ncol=cols(betaB));
    
    if(!scalmiss(c1))
      betaW[c1,] <- matrix(NA, nrow=rows(c1),ncol=cols(betaB));
    
    
  }else{
    if (col!=cols(betaB))
      stop("EI internal error, betab2w")
     
    for( i in 1:col){
      cs<-homoindx(x[,i])
      c<-cs$c
      c0<-cs$c0
      c1<-cs$c1
     
      if(!scalmiss(c))
        betaW[c,i] <- (t[c]/x1[c,i])-((betaB[c,i]*x[c,i])/x1[c,i]);
      
      if(!scalmiss(c0))
        betaW[c0,i] <- matrix(t[c0],nrow=rows(c0),ncol=1);
      
      if(!scalmiss(c1))
        betaW[c1,i] <- matrix(NA,nrow=rows(c1),ncol=1);
      
    }
    
  }
  return(betaW)
}
