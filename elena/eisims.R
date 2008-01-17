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
  t <- T
  x <- X
  
  Eprt <- get("Eprt", env=evbase)
  Eeta <- get("Eeta", env=evbase)
  Erho <- get("Erho", env=evbase)
  Esims <- get("Esims", env=evbase)
  EisChk <- get("EisChk", env=evbase)
  Ebetaws <- EbetaWs <- get("EbetaWs", env=evbase)
  Eres <- get("Eres", env=evbase)
  ##  /* ESTIMATION VARIABILITY */
  if(Eprt>=2)
    message("Simulating estimation variation...");
 
 ## /* _Eselect is applied within packdta, so dataset is _Eselect'd but
  ##    x,Zb,Zw,t, which are used for fundamental variability, are not */
  dataset <- packdta(x,Zb,Zw,t,evbase);
  
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
    assign("EisFix", EisFix, env=evbase)
    yy1 <- rndisamp(eiloglik,trimr(MLpsi,0,3), MLvc[1:(rows(MLpsi)-3), 1:(rows(MLpsi)-3)],
                    dataset,Esims,EisFix, evbase) 
    yy2 <-  rndmn(Erho[2],0,Esims)  
    yy3 <-  rndmn(etaB,etaBs,Esims)
    yy4 <- rndmn(etaW,etaBs,Esims)
    PSIsims <- cbind(yy1, yy2,yy3, yy4)
                    
                    
  }else{
   
    EisFix <- MLpsi[(rows(MLpsi)-1):rows(MLpsi)]
    assign("EisFix", EisFix, env=evbase)
    yy1 <- rndisamp(eiloglik,trimr(MLpsi,0,2),
                    MLvc[1:(rows(MLpsi)-2),1:(rows(MLpsi)-2)],
                    dataset,Esims,EisFix,evbase)
     
    yy2 <- rndmn(etaB,etaBs,Esims)  
    yy3 <-   rndmn(etaW,etaBs,Esims)
    PSIsims <- cbind(yy1, yy2, yy3)
                   
                   
  }
  psisims <-  as.matrix(PSIsims)
  if(EisChk)                          
    Eres <- vput(Eres,psisims,"PhiSims") ###  @ save name corresponds to book @ 
  else 
    Eres <- vput(Eres,cbind(colMeans(psisims),sd(as.data.frame(psisims))),"PhiSims")
  
  
 ### /* FUNDAMENTAL VARIABILITY */
 ### /* Applies to all observations */
  if(Eprt>=2)
    message("Simulating fundamental variability...")
  
  rs <- rows(x);
  betaBs <- matrix(0, nrow=rs,ncol=Esims)
  betaWs <- betaBs
  lst <- bounds1(t,x,tvap,get("EnumTol", env=evbase)) 
  Bbeta <- lst[[1]]
  aggs <- lst[[2]]
  for( k in (1:Esims)){
    y <- eirepar(PSIsims[k+0,],Zb,Zw,x,get("Ez", env=evbase),evbase)
    lst <- exvar(T,x,y$Bb,y$Bw,y$sb,y$sw, y$rho,evbase)
    tt <- lst[[1]]
    tt <- lst[[2]]
    tt <- lst[[3]]
    tt <- lst[[4]]
    Ebb <- lst[[5]]
    Vbb <- lst[[6]]
    betaBs[,k+0] <- rndtni(Ebb,Vbb,Bbeta[,1:2],evbase=evbase) 
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

  if(!length(evbase))
    evbase <- get("evbase", env=parent.frame())
  col <- cols(x);
  x1 <- 1-x;
  betaB <- betab <- as.matrix(betab)
  betaW <- betaw <- betab;
  
  if(col==1){
   
     cs<-homoindx(x)
     c<-cs$c
     c0<-cs$c0
     c1<-cs$c1
    
    if(!scalmiss(c))
      betaW[c,] <- (as.matrix(t[c])%dot/%as.matrix(x1[c]))%plus%-
        as.matrix(betaB[c,])%dot*% (as.matrix(x[c])%dot/%as.matrix(x1[c]));
    
    if(!scalmiss(c0))
      betaW[c0,] <- t[c0] %dot*% matrix(1,nrow=rows(c0),ncol=cols(betaB));
    
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
        betaW[c,i] <- (t[c]%dot/%x1[c,i])-((betaB[c,i]%dot*%x[c,i])%dot/%x1[c,i]);
      
      if(!scalmiss(c0))
        betaW[c0,i] <- t[c0] %dot*% matrix(1,nrow=rows(c0),ncol=1);
      
      if(!scalmiss(c1))
        betaW[c1,i] <- matrix(NA,nrow=rows(c1),ncol=1);
      
    }
    
  }
  return(betaW)
}
##/* -----------------------------------------------------------------------
##   y = rndisamp(&f,b,vc,dataset,sims);
##**
##** random numbers generated via importance sampling
##** using the multivariate normal or t as the first approximation.
##** 
##** inputs: f() = likelihood function with arguments b(kx1),dataset
##**         mu = kx1 means
##**         vc = kxk variance matrix
##**    dataset = input to f()
##**       sims = number of simulations
##**
##** output:  y = nxk matrix of dependent Multivariate Normal Random Variables
##**              each row of y is one 1xk simulation
##**
##** globals:  _Eprt if >=2 print progress reports
##**           _Eisn = sims*_Eisn is the number of normals to draw before 
##**                   resampling (default _Eisn=10)
##**           _EisFac = factor to multiply normal variance matrix by
##**                     or -1 to do normal approximation only
##**                     or -2 to use maximum posterior estimates
##**           _EisChk = 1 if save ln(imptce Ratio)~phi in lnir;
##**                     0 otherwise save means
##**           _Eist   = 0 (default) for normal random numbers; 
##**                    >2 for multivariate t, with _Eist degrees of freedom
##**           _EiLlikS = 1 save "lliksims" log-likelihood at each simulation 
##**                      in dbuf; 0 save mean(lliksims) in dbuf
##*/
rndisamp <- function(f,b,vc,dataset,sims,EFix,evbase=get("evbase",env=parent.frame())){

  if (rows(vc)!=rows(b) || cols(vc)!=rows(b))
    stop("rndisamp: input error")
   
  getEnvVar(evbase, environment())
  cml.bounds <- get("cml.bounds", env=evbase)
  EisFac <- get("EisFac", env=evbase)
  Ghactual <- GhActual <- get("GhActual", env=evbase)
  Eres <- get("Eres", env=evbase)
  if(as.vector(EisFac)==-1) { ###                     @ normal approximation only @
    Eres <- vput(Eres,NA,"lnir")
    Eres <- vput(Eres,NA,"meanIR")
    Eres <- vput(Eres,NA,"resamp")
    ei.vc <- get("ei.vc", env=evbase)
    if(ei.vc[Ghactual,1]!=-1)
      return(rndmn(b,vc,sims))
    else{
 
      return(rndtsn(b,vc,sims,cml.bounds[1:rows(b),],1e-3, Eprt))
    }
  }
  if (as.vector(EisFac)==-2) ###                @ use maximum posterior estimates @
    return(t(as.matrix(b))%dot*% matrix(1,nrow=Esims,ncol=1))
  else{
    if(ei.vc[Ghactual,1]!=-1)
      vc <- vc/as.vector(EisFac)  ###@ "vc" is -Hessian @
    else
      vc <- vc*as.vector(EisFac) ###  @ "vc" is a variance matrix @
  
  }
  k <- sims*Eisn;
  if (as.vector(Eist)==0){
    if (ei.vc[Ghactual,1]!=-1){
      
      psis <- rndmn(b,vc,k)
      ###print(psis)
    }else{ ###                     @ use when _ei_vc={-1 0} @
      psis <- rndtsn(b,vc,k,cml.bounds[1:rows(b),],1e-3, Eprt) 
      if (any(is.na(psis)) && as.vector(Eprt)>=2)
        message("ei: _EI_vc={-1 0} option failed.")
    }
   
  }else{
    vc <- vc*as.vector((Eist-2)/Eist)
    psis <- rndmt(b,vc,Eist,k) 
  }
  
  lik <- as.matrix(0,nrow=k)
  norm <- as.matrix(0,nrow=k)
  keepsi <- as.vector(0*b)
  keeplik <- 0
  if (ei.vc[Ghactual,1]!=-1){
   
    Eivc <- invpd(vc)
  }else
    Eivc <- vc
 
 
  for (i in 1:k){
    if ((i%%100)==0 && Eprt>=2)
     print( i/10);
       
    bb <- as.matrix(c(as.vector(psis[i,]),EFix))
    
    lik[i] <- colSums(as.matrix(na.omit(f(bb,dataset))))
      
    
    if (as.vector(Eist)==0){
      norm[i] <- lnpdfmn2(psis[i,],b,0,Eivc) ## var from _Eivc above @
    }else{
    
      norm[i] <- lnpdfmt(psis[i,.],b,0,Eist,Eivc) 
    }
  }
  lnir <- lik-norm;	###@ ln(imptce ratio)  @ 
  
  
  if (scalone(EisChk))
    Eres <- vput(Eres,cbind(lnir,psis),"lnir")
  else{
    max <- maxc(lnir)
    meanIR <- max+log(colMeans(as.matrix(exp(lnir-max)))) ###    @ = ln(meanc(exp(lnir)) @
    Eres <- vput(Eres,meanIR,"lnir")
  }

  if (any(is.na(lnir))){
    message("Importance sampling failed! Using maximum posterior estimates\n", 
    "(equivalent to _EisFac=-2). ")
    Eres <- vput(Eres,-1,"lliksims")
    Eres <- vput(Eres,-1,"resamp")
    return(res <- t(as.matrix(b))%dot*% matrix(1,nrow=Esims,ncol=1))
  }
  
  resamp <- 0;
  r <- 0;
  while( (r-1)<sims){
    lnir <- lnir-maxc(lnir)
    ir <- exp(lnir) ###	        @ imptce ratio  @
    rnd <- matrix(runif(rows(ir)), nrow=rows(ir), ncol=1, byrow=TRUE)
    tst <- rnd <=ir
    keepsi <- rbind(keepsi, selif(psis,tst))
    keeplik <- rbind(keeplik,selif(lik,tst))
    r <- rows(keepsi)
    lnir <- subset(lnir, subset=!tst)
    if(!length(lnir)) lnir <- NA
    psis <- subset(psis,subset=!tst)
    if(!length(psis)) psis <- NA
    lik <- subset(lik,subset=!tst)
    if(!length(lik)) lik <- NA
    if (Eprt>=2){
      resamp <- resamp+1
      message("Got ", r-1, ", resampling...")
    }
  }
  
  Eres <- vput(Eres,resamp,"resamp")
 
  keepsi <- keepsi[2:(sims+1),]
  keeplik <- keeplik[2:(sims+1),]
  if (as.vector(get("EiLliks", env=evbase))==1)
    Eres <- vput(Eres,keeplik,"lliksims")
  else
    Eres <- vput(Eres,colMeans(as.matrix(keeplik)),"lliksims")
 
  return(keepsi)
}
