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
    yy1 <- rndisamp(eiloglik,trimr(MLpsi,0,3), MLvc[1:rows(MLpsi)-3, 1:rows(MLpsi)-3],dataset,Esims,evbase) ###COMPUTE
    yy2 <-  rndmn(Erho[2],0,Esims)  
    yy3 <-  rndmn(etaB,etaBs,Esims)
    yy4 <- rndmn(etaW,etaBs,Esims)
    PSIsims <- cbind(yy1, yy2,yy3, yy4)
                    
                    
  }else{
    EisFix <- MLpsi[(rows(MLpsi)-1):rows(MLpsi)]
    yy1 <- rndisamp(eiloglik,trimr(MLpsi,0,2),
                    MLvc[1:rows(MLpsi)-2,1:rows(MLpsi)-2],
                    dataset,Esims,evbase) ##COMPUTE
    yy2 <- rndmn(etaB,etaBs,Esims)  
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
rndisamp <- function(f,b,vc,dataset,sims,evbase=get("evbase",env=parent.frame())){
 ## local f:proc;
 ## local lik,norm,i,max,psis,u,tst,r,lnir,ir,k,radd,keepsi,inf,ii,x,t,bb,meanir,resamp,vc0,keeplik;
  
  if (rows(vc)!=rows(b) || cols(vc)!=rows(b))
    stop("rndisamp: input error")
   
   evloc <- getEnvVar(evbase, environment())
  
  if(as.vector(EisFac)==-1) { ###                     @ normal approximation only @
    Eres <- vput(Eres,NA,"lnir")
    Eres <- vput(Eres,NA,"meanIR")
    Eres <- vput(Eres,NA,"resamp")
    Ghactual <- GhActual
    if(ei.vc[Ghactual,1]!=-1)
      return(rndmn(b,vc,sims))
    else
      return(rndtsn(b,vc,sims,cml.bounds[1:rows(b),],1e-3)) ###COMPUTE
  }
  if (as.vector(EisFac)==-2) ###                @ use maximum posterior estimates @
    return(as.vector(b)%dot*%as.matrix(1,nrow=Esims,1))
  
  if(ei.vc[Ghactual,1]!=-1)
      vc <- vc/EisFac  ###@ "vc" is -Hessian @
  else
      vc <- vc*EisFac ###  @ "vc" is a variance matrix @
   
  
  
  k <- sims*Eisn;
  if (as.vector(Eist)==0){
    if (ei.vc[Ghactual,1]!=-1)
      psis <- rndmn(b,vc,k)
    else{ ###                     @ use when _ei_vc={-1 0} @
      psis <- rndtsn(b,vc,k,cml.bounds[1:rows(b),],1e-3) ###COMPUTE
      if (any(is.na(psis)) && as.vector(Eprt)>=2)
        message("ei: _EI_vc={-1 0} option failed.")
    }
   
  }else{
    vc <- vc*as.vector((Eist-2)/Eist)
    psis <- rndmt(b,vc,Eist,k) ###COMPUTE
  }
  lik <- as.matrix(0,nrow=k)
  norm <- as.matrix(0,nrow=k)
  keepsi <- as.vector(0*b)
  keeplik <- 0
  if (ei.vc[Ghactual,1]!=-1)
    Eivc <- invpd(vc)
  else
    Eivc <- vc
 
  
  for (i in 1:k){
    if int(i/100)==i/100 and _Eprt>=2;
      i/10;
    endif;
    bb=(psis[i,.]'|_EisFix);
    lik[i]=sumc(packr(f(bb,dataset)));
    if _Eist==0;
      norm[i]=lnpdfmn2(psis[i,.]',b,0); @ var from _Eivc above @
    else;
      norm[i]=lnpdfmt(psis[i,.]',b,0,_Eist);
    endif;
  endfor;
  lnir=lik-norm;	@ ln(imptce ratio)  @

  if scalone(_EisChk);
    _Eres=vput(_Eres,lnir~psis,"lnir");
  else;
    max=maxc(lnir);
    meanIR=max+ln(meanc(exp(lnir-max)));    @ = ln(meanc(exp(lnir)) @
    _Eres=vput(_Eres,meanIR,"lnir");
  endif;

  if ismiss(lnir)==1;
    "Importance sampling failed! Using maximum posterior estimates";
    "(equivalent to _EisFac=-2). ";
    _Eres=vput(_Eres,-1,"lliksims");
    _Eres=vput(_Eres,-1,"resamp");
    retp(b'.*ones(_Esims,1));
  endif;
  
  resamp=0;
  r=0;
  do while r-1<sims;
    lnir=lnir-maxc(lnir);
    ir=exp(lnir);	        @ imptce ratio  @
    tst=(rndu(rows(ir),1).<=ir);
    keepsi=keepsi|selif(psis,tst);
    keeplik=keeplik|selif(lik,tst);
    r=rows(keepsi);
    lnir=delif(lnir,tst);
    psis=delif(psis,tst);
    lik=delif(lik,tst);
    if _Eprt>=2;
      resamp=resamp+1;
      "Got " r-1 ", resampling...";
    endif;
  endo;
  _Eres=vput(_Eres,resamp,"resamp");
  
  keepsi=keepsi[2:sims+1,.];
  keeplik=keeplik[2:sims+1,.];
  if _EiLlikS==1;
    _Eres=vput(_Eres,keeplik,"lliksims");
  else;
    _Eres=vput(_Eres,meanc(keeplik),"lliksims");
  endif;    
  retp(keepsi);
endp;  
