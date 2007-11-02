##/*
##  This archive is part of the program EI
##  (C) Copyright 1995-2001 Gary King
##  All Rights Reserved.
##*/
##/*
##    {b,vc}=quadcml(x,Zb,Zw,y);
##  Ecological Inference Likelihood Maximization via CML
##  With variance-covariance computed via global methods.
## 
##  MODEL:
##  mean function:  E(Y)=Bb*X + Bw*(1-x)
##  var  function:  V(Y)=vb*X^2 + c(b,w)*2*X*(1-X) + vw*(1-X)^2
##                  var params reparam'd: vars>0 and p.d.
##  distribution:   distribution implied on Y if Bb,Bw are bivariate 
##                  truncated normal
##  INPUT:
##  x = explanatory variable 
##  Zb = 1 or covariates for Bb (constant term included automatically)
##  Zw = 1 or covariates for Bw (constant term included automatically)
##  y = dependent variable
##  OUTPUT:
##  b = {Bb, Bw, sb, sw, rho}, where Bb, Bw are vectors if Zb,Zw have vars
##  vc = estimation global var cov matrix of b
##**
##  where params of interest are reparameterized as eirepar.g
##**
##  GLOBALS
##  _Erho[1]=standard deviation of normal prior on phi_5 (default=0.5)
##        0 fix it to _Erho[2] and don't estimate
##       <0 estimate without prior
##  _Ebounds=1 set bounds automatically unless z's are included
##           0 don't use bounds
##           kx2 or 1x2 matrix to indicate upper~lower bounds 
##  _Eprt=0 print nothing
##        1 print only final output from each stage
##        2 print everything useful plus friendly iteration numbers etc
##        3 print everything useful, iterations, and all sorts of checks
##  _Estval = 1 use default starting values or set to vector
##  _Eselect = scalar 1 to use all observations
##             vector of 1's to select and 0's to discard observations during
##             likelihood estimation.

#include ei.ext;
quadcml <- function(x,Zb,Zw,y,evbase=parent.frame()) {
   x <- matrix(x,ncol=1)
   y <- matrix(y,ncol=1)
  
   evloc <- getEnvVar(evbase, environment())###, vecvar=c("eimetar"))
   if (Eprt>=2)
    message("Likelihood estimation...");
 
  Edirtol <- EdirTol
  ###/* prepare data */
  dataset <- packdta(x,Zb,Zw,y);
  if(Eprt>=2){
    if (EselRnd !=1){
    ###  fmtt;
      message("EselRnd: Modifying _Eselect with random selection = ", fmtt(EselRnd))
    }
    tt <- colSums(as.matrix(1- Eselect))
    if (rows(Eselect)!= 1 && tt!=0){
  ###    fmtt;
      message("Eselect: deleting ", fmtt(tt) ," ", 
              "observations from estimation stage;", colSums(as.matrix(Eselect)),
              "observations remain")
    }
  }
  
###  /* starting values OR grid search */
  
   if (rows(Estval)==1 && !(scalone(Estval))){###  @ grid search @
     gridl <- ifelse (Estval==0, 5, Estval)
   }else{
     gridl <- 0 ###        @ no grid search @
   }
###    gridl <- 5 @ default number of gridlines for grid search @   
   
 
   stval <- Estval;
   if (scalone(Estval) || gridl !=0)
     stval <- as.matrix(c(rep(0,colSums(as.matrix(Ez))), -1.2,-1.2,0)) 
  
   if (rows(Eeta)==4)
    stval <- as.matrix(c(stval, Eeta[1:2]))
   else if(Eeta[1]==4)
    stval <- as.matrix(c(stval,c(0,Eeta[2])))
   else if(Eeta[1]==5)
    stval <- as.matrix(c(stval,c(Eeta[2],0)))
   else
    stval <- as.matrix(c(stval, c(0,0)))
 
  
 ### /* cml globals */
   title <- "EI Likelihood Maximization: "
   if(vin(Eres,"titl"))
     title <- paste(title,eiread(Eres,"titl"))
   cml.Active <- as.matrix(c(rep(1,rows(stval)-2),0,0))
   assign("cml.Active", cml.Active, env=evbase)
   cml.MaxIters <- Emaxiter <- EmaxIter;
   assign("cml.MaxIters", cml.MaxIters, env=evbase)
   cml.DirTol <- EdirTol;
   assign("cml.DirTol", cml.DirTol, env=evbase)
   cml.CovPar <- 0
   assign("cml.CovPar", cml.CovPar, env=evbase);
   cml.ParNames <- matrix("", nrow=7)
   spacer <- "   "
   
   cml.ParNames[1] <- paste(spacer,"Zb",ftosm(seq(from=0, by=1,length.out=Ez[1]),1,0),sep="")
   cml.ParNames[2] <- paste(spacer,"Zw", ftosm(seq(from=0,by=1,length.out=Ez[2]),1,1),sep="")
   cml.ParNames[3] <- paste(spacer,  "sigB", sep="")
   cml.ParNames[4] <- paste(spacer,"sigW",sep="")
   cml.ParNames[5] <- paste(spacer,"rho",sep="")
   cml.ParNames[6] <- paste(spacer,"etaB",sep="")
   cml.ParNames[7] <- paste(spacer, "etaW",sep="")
   assign("cml.ParNames",cml.ParNames, env=evbase)
###@ if change this code, change also eiread @
  if (scalzero(Ebounds))###     don't use bounds  
    cml.bounds <- cml.Bounds <- matrix(c(-1e+256, 1e+256), nrow=1, ncol=2)
  else if (cols(Ebounds)==2)
    cml.bounds <- cml.Bounds <- matrix(c(Ebounds,0,0.0001,0,0.0001), ncol=2,byrow=TRUE)
  else if (scalone(Ebounds)) ###    automatic bounds calculation 
    {
      bnds <- matrix(c(-10,10),nrow=1)
      nbnds <- matrix(c(-20,20),nrow=1)
      if (Ez[1]==1)
        cml.bounds <- cml.Bounds <- bnds
      else
        cml.bounds <- cml.Bounds <- nbnds %dot*% matrix(1,nrow=Ez[1],ncol=1)
    
      if (Ez[2]==1)
        cml.bounds <- cml.Bounds <- rbind(cml.Bounds,bnds)
      else
        cml.bounds <- cml.Bounds <- rbind(cml.Bounds,(nbnds %dot*% matrix(1,nrow=Ez[2],ncol=1)))
    
      cml.bounds <- cml.Bounds <- matrix(c(cml.Bounds,-6, 3,-6,3,-2,2), ncol=2, byrow=TRUE)
    }else{
    assign("cml.Bounds", cml.Bounds, env=evbase)
    stop(message="quadcml: problem with Ebounds")
    
  }
  assign("cml.Bounds", cml.Bounds, env=evbase)
  if (gridl==0 && cols(Ebounds)!=2){##  @ sneak past CML checks on parameters @
    tt <- rows(stval)
    tt <- matrix(c(tt-1,tt), nrow=2)
    cml.Bounds <- matrix(c(as.vector(cml.Bounds),stval[tt],stval[tt]), ncol=2,byrow=TRUE) 
    tt <- rows(cml.Bounds)
    tt <- matrix(c(tt-1,tt),nrow=2)
    cml.Bounds[tt,2] <- cml.Bounds[tt,2]+ Edirtol
    assign("cml.Bounds", cml.Bounds, env=evbase)
  }
  
  if (Eprt==0)
    output <- 0
  else if (Eprt<=2){
    cml.Diagnostic <- 2
    output <- 1
  }else if (Eprt==3){
    cml.Diagnostic <- 3
    output <- 2
  }else{
    assign("cml.Diagnostic", cml.Diagnostic, env=evbase)
    stop(message="quadcml: problem with Eprt")
  }
    assign("cml.Diagnostic", cml.Diagnostic, env=evbase)
  
  if (Erho[1]==0){
    r <- rows(stval)-2
    stval[r] <- Erho[2]
    cml.Bounds[r,] <- cbind((Erho[2]-macheps),(Erho[2]+macheps))
    assign("cml.Bounds", cml.Bounds, env=evbase)
    cml.Active[r] <- 0
    assign("cml.Active", cml.Active, env=evbase)
  }

  if (gridl==0){###         /* run CML */
   
 ###   lst <- cml(dataset,0,&eiloglik,stval)
 ###   {b,mlogl,grds,vc,ret}= cml(dataset,0,&eiloglik,stval)
 ###   b <- lst[[1]]
 ###   mlogl <- lst[[2]]
 ###   grds <- lst[[3]]
 ###   vc <- lst[[4]]
 ###   ret <- lst[[5]]
 ###  if (ret %inG% c(3,4,6,7,10,99)){
      message( "?")
      message("*********************************************")
      message("CML return code: ");
      message("CML (Constraint Maximum Likelihood) OPTION not avaliable")
###      print(ret)
      message("Restarting iterations with trust algorithm on")
      message("*********************************************")
###  cml.options <- matrix(trust)
###  {b,mlogl,grds,vc,ret}=cml(dataset,0,&eiloglik,stval);
###   lst <- cml(dataset,0,&eiloglik,stval)
###   b <- lst[[1]]
###   mlogl <- lst[[2]]
###   grds <- lst[[3]]
###   vc <- lst[[3]]
###   ret <- lst[[4]]
###   cml.options <- 0
###   assign("cml.options", cml.options, env=evbase)
###    }
###    Eres <- vput(Eres,ret,"retcode");
###    logl <- mlogl*cml.NumObs;
###    Eres <- vput(Eres,logl,"loglik");
###    if (Eprt>=2)
###      message("CML converged; Computing variance-covariance matrix...?")
    }
   
   if (gridl!=0){ ###/* run GRID search */
     if (Eprt>=2)
       message("Preliminary mean grid search (on 2 parameters)...")
     Edirtol <- EdirTol
     tt <- rows(cml.bounds)
     gr1 <- colMeans(t(cml.bounds[3:tt,]))
     
     gr1 <- rbind(cml.bounds[1:2,], matrix(c(gr1,gr1),ncol=2))
     lst <- eigrid(dataset,gr1,53,EdirTol) ###COMPUTE 
     b <- lst[[1]]
     logl <- lst[[2]]
     cml.bounds[1:2,] <- b[1:2]+matrix(c(-1,-1, 1,1),nrow=2, ncol=2)
    if (Eprt>=2)
      message("?; Main grid search (on all parameters)...")
     lst <- eigrid(dataset,cml.bounds,gridl,EdirTol)
     b <- lst[[1]]
     logl <- lst[[2]]
    
    ret <- 3333
    Eres <- vput(Eres,ret,"retcode")
    Eres <- vput(Eres,logl,"loglik")
    grds <- b*0+Edirtol
   }
  
  ###/* compute var cov matrix */
  if (EdoSim !=-1){
   
    GhFix <- as.matrix(c(rep(1,(rows(stval)-2)),0,0))
    if(Erho[1]==0)
      GhFix[r] <- 0
    assign("GhFix",GhFix, env=evbase)  
    vc <- gvc(eiloglik,b,dataset) ###COMPUTE
    if (Erho[1]==0){
      vc <- cbind(vc, as.matrix(rep(0,rows(vc))))
      vc <- rbind(vc, matrix(rep(0,cols(vc)), nrow=1))
      vc[rows(vc),cols(vc)]=0.00000000001
    }
   Eres <- vput(Eres,GhActual,"GhActual");
    if (scalmiss(vc)){
      message("quadcml: couldn't compute positive definite variance matrix")
      return(c(NA,NA))
    }
  }else
    vc <- matrix(0, nrow=(rows(b)-2),ncol=(rows(b)-2))
  
  
###  /* print likelihood results */
  if(Eprt>=1 || ret!=0){
    if (ret!=0 && ret!=3333)
      message("Overriding no print instructions.  See CML return code")
   
    if (ret==20)
      ret <- 1
  ###  formatC(x, digits=7, width=4)
  ###  format/ro 7,4;
    cat("?\n")
    message(title)
    message("CML Return code:          ", formatC(ret,digits=7, width=4))
    message("log-likelihood:           ", formatC(logl,digits=7, width=4))
    message("Number of observations:   ", formatC(cml.NumObs, digits=7, width=4))
    message("Number of iterations:     ",formatC(cml.IterData[1],digist=7,width=4))
    if(ei.vc[Ghactual,1]!=-1)
      message("Variance computed from:   ", formatC(ei.vc[Ghactual,],digits=7,width=4))
    else{
      message("-Hessian computed from:   ")
      pt <- sapply(ei.vc[Ghactual,], formatC, digits=7, width=4)
      print(pt)
    }
    cat("?\n")
    message("*## Parameter values in scale of estimation ***")
    mask <- c(0,1,1)
    fmt <- matrix("", nrow=3, ncol=3)
    fmt[1,] <- c("*.*s",  1,8)
    fmt[2,] <- c("*.*lf",10,4)
    fmt[3, ] <- c("*.*lf", 10, 4)
    if (all(ei.vc[Ghactual,]!=-1))
      se <- sqrt(matrix(c(extract.diag(vc), 0,0), ncol=1))
    else
      se <- matrix(NA, nrow=rows(b),ncol=1)
    
    printfm(c(cml.ParNames,b,se),mask,fmt) ###COMPUTE
    cat("?\n")
    lst <- eirepar(b,Zb,Zw,x) ###COMPUTE
    bb <- lst[[1]]
    bw <- lst[[2]]
    sb <- lst[[3]]
    sw <- lst[[4]]
    rho <- lst[[5]]
    vrs <- as.matrix(c(bb, bw, sb, sw, rho))
    message("Reparameterization back to the truncated scale, parameterized",
    "        according to the underlying untruncated distribution.")
    print(as.vector(vrs))
    pb <- c(colMeans(bb),colMeans(bw),sb,sw,rho)
    print(pb)
    pb <- as.matrix(pb)
    if (Eprt>=3){
      cat("?\n")
      message("Reparameterization back to ultimate truncated scale")
      pb <- eirepart(b,Zb,Zw,x)  ###COMPUTE
      print(as.vector(vrs))
      print(as.vector(pb))
    }
  }  

  Eres <- vput(Eres,cml.Parnames,"parnames")
  Eres <- vput(Eres,b,"phi")
  Eres <- vput(Eres,vc,"vcphi")
  lst <- c(list(b=b),list(vc=vc))
  return(lst)
}

