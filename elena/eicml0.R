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
quadcml <- function(x,Zb,Zw,y,evbase=parent.frame(),macheps=2.23e-16, algor="BFGS") {
   x <- matrix(x,ncol=1)
   y <- matrix(y,ncol=1)
   if(!length(macheps)) macheps <- .Machine$double.eps
  
  evloc <- getEnvVar(evbase, environment())###, vecvar=c("eimetar"))
   if (Eprt>=2)
     message("Likelihood estimation...");
 
  Edirtol <- EdirTol  ### = .Machine$double.eps * con$factr  
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

  
   if (scalone(Estval) || gridl !=0)
     stval <- as.matrix(c(rep(0,colSums(as.matrix(Ez))), -1.2,-1.2,0)) 
   else
      stval <- Estval;
   if (rows(Eeta)==4)
    stval <- rbind(stval, Eeta[1:2])
   else if(Eeta[1]==4)
    stval <- rbind(stval,0,Eeta[2])
   else if(Eeta[1]==5)
    stval <- rbind(stval,Eeta[2],0)
   else
    stval <- rbind(stval, 0,0)
 
  
 ### /* cml globals */
   title <- "EI Likelihood Maximization: "
   if(vin(Eres,"titl"))
     title <- paste(title,eiread(Eres,"titl"))
 ###  cml.Active <- as.matrix(c(rep(1,rows(stval)-2),0,0))
   cml.MaxIters <- Emaxiter <- EmaxIter; ###=con$maxit
   cml.DirTol <- EdirTol;
 ###  cml.CovPar <- 0
   cml.ParNames <- matrix("", nrow=7)
   spacer <- "   "
   
   cml.ParNames[1] <- paste(spacer,"Zb",ftosm(seq(from=0, by=1,length.out=Ez[1]),1,0),sep="")
   cml.ParNames[2] <- paste(spacer,"Zw", ftosm(seq(from=0,by=1,length.out=Ez[2]),1,1),sep="")
   cml.ParNames[3] <- paste(spacer,  "sigB", sep="")
   cml.ParNames[4] <- paste(spacer,"sigW",sep="")
   cml.ParNames[5] <- paste(spacer,"rho",sep="")
   cml.ParNames[6] <- paste(spacer,"etaB",sep="")
   cml.ParNames[7] <- paste(spacer, "etaW",sep="")
   cml.ParNames <- sapply(cml.ParNames,trim.blanks)


###@ if change this code, change also eiread @
  if (scalzero(Ebounds))###     don't use bounds  
    cml.bounds <- cml.Bounds <- matrix(c(-1e+256, 1e+256), nrow=1, ncol=2)
  else if (cols(Ebounds)==2)
    cml.bounds <- cml.Bounds <- rbind(Ebounds,cbind(0,0.0001),cbind(0,0.0001))
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
    
      cml.bounds <- cml.Bounds <- rbind(cml.Bounds,cbind(-6, 3),cbind(-6,3),cbind(-2,2))
    }else{
     
    stop(message="quadcml: problem with Ebounds")
    
  }
 
   
  if (gridl==0 && cols(Ebounds)!=2){##  @ sneak past CML checks on parameters @
    tt <- rows(stval)
    tt <- matrix(c(tt-1,tt), nrow=2)
    cml.bounds <- cml.Bounds <- rbind(cml.Bounds,cbind(stval[tt],stval[tt])) 
    tt <- rows(cml.Bounds)
    tt <- matrix(c(tt-1,tt),nrow=2)
    cml.bounds[tt,2] <- cml.Bounds[tt,2] <- cml.Bounds[tt,2]+ Edirtol
   
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
    
    stop(message="quadcml: problem with Eprt")
  }
    
  
  if (Erho[1]==0){
    r <- rows(stval)-2
    stval[r] <- Erho[2]
    cml.Bounds[r,] <- cbind((Erho[2]-macheps),(Erho[2]+macheps))
   cml.bounds <- cml.Bounds
    cml.Active[r] <- 0
   
  }
       
   assign("cml.Bounds", cml.Bounds, env=evbase)
   assign("cml.bounds", cml.Bounds, env=evbase)
  if (gridl==0){###         /* run CML */
   
 ###   lst <- cml(dataset,0,&eiloglik,stval)
 ###   {b,mlogl,grds,vc,ret}= cml(dataset,0,&eiloglik,stval)
 ###   theta --or-- par = stval
 ###   f=&eiloglik
 ###      constrOptim  
 
   ff <- function(x,y,ev=evbase){
     res <- eiloglik(x,y,ev)
     return(colSums(as.matrix(res)))
   }
     
   
   stval <- as.matrix(stval)
   rownames(stval) <- cml.ParNames
   sval <<- stval
   cmlb <<- cml.bounds
   ddta <<- dataset
   ###defaults
   par <- stval
   con <- list(trace = 0, fnscale = 1, parscale = rep.int(1, 
        length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100, 
        abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
        beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
        factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
   ###changes
   con$trace <- 1
   con$fnscale <- -1 ##maximizes
   con$REPORT <- 1
   con$factr <- 1e+11
   
   message("Covergence acuracy is ", .Machine$double.eps*con$factr);
  ###  delta.bounds <- cml.bounds[,2] - cml.bounds[,1]
  ###  ix <- which(delta.bounds == max(delta.bounds),arr.ind=TRUE)
  ###  con$parscale <-rep.int(1,length(par))
  ##  con$parscale[ix] <- 0.1
    ###faster convergence increase con$factr, i.e.  <- 1e+08
    ###tolerance is defined as .Machine$double.eps*con$factr
    ix <- match(c("reltol", "abstol"), names(con))
    con <- con[-ix]
    message("Optim in action....")
   optimlst <- optim(stval,ff,gr=NULL,dataset,evbase,method="L-BFGS-B",
               control=con, lower=cml.bounds[,1],upper=cml.bounds[,2],hessian=FALSE)
   assign("optimlst", optimlst, env=evbase)
   optimnm <- names(optimlst)
   b <-optimlst$par
   mlogl <-optimlst$value
   grds <- NA
   cntgrd <- optimlst$counts[[2]]
   cntfn <- optimlst$counts[[1]]
   ret <- optimlst$convergence
   hess <- if("hessian" %in% optimnm) optimlst$hessian 
   if(ret >=1 ){
     message("Optim did not converge or produce an error...")
     print(optimlst)
     message("con$factr= ", con$factr)
     message("con$parscale= ", con$parscale)
     message("con$fnscale= ", con$fnscale)
     stop("Change defaults")
   }
   print(optimlst)
   message("Calculating hessian ...")
   if(is.null(hess)){
     hess <- optimhess(par=b,fn=ff,gr=NULL,dataset,control=con,nm=cml.ParNames)
     assign("hessian", hess,env=evbase)
     vc <- inv(hess)
   } 
 }
 
###    Eres <- vput(Eres,ret,"retcode");
###    logl <- mlogl*cml.NumObs;
###    Eres <- vput(Eres,logl,"loglik");
###    if (Eprt>=2)
###      message("CML converged; Computing variance-covariance matrix...?")
    }
 
 ###  if (gridl!=0){ ###/* run GRID search */
     if (Eprt>=2)
       message("Preliminary mean grid search (on 2 parameters)...")
     Edirtol <- EdirTol
     tt <- rows(cml.bounds)
     gr1 <- colMeans(t(cml.bounds[3:tt,]))
     
     gr1 <- rbind(cml.bounds[1:2,], matrix(c(gr1,gr1),ncol=2))
   ###TESTING HERE
     lst <- eigrid(dataset,gr1,53,EdirTol)
    
     b <- lst[[1]]
     logl <- lst[[2]]
     cml.bounds[1:2,] <- b[1:2]%plus% matrix(c(-1,-1, 1,1),nrow=2, ncol=2)
    
    if (Eprt>=2)
      message("?; Main grid search (on all parameters)...")
     lst <- eigrid(dataset,cml.bounds,gridl,EdirTol) 
     b <- lst[[1]]
     logl <- lst[[2]]
     
    ret <- 3333
    Eres <- vput(Eres,ret,"retcode")
    Eres <- vput(Eres,logl,"loglik")
    grds <- b*0+Edirtol
 ###  }
  
  ###/* compute var cov matrix */
  if (EdoSim !=-1){
 
    GhFix <- as.matrix(c(rep(1,(rows(stval)-2)),0,0))
    if(Erho[1]==0)
      GhFix[r] <- 0
    assign("GhFix",GhFix, env=evbase)
   
    vc <- gvc(eiloglik,b,dataset)
 
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
    mask <- matrix(c(0,1,1),nrow=1)
    fmt <- matrix("", nrow=3, ncol=3)
    fmt[1,] <- c("*.*s",  1,8)
    fmt[2,] <- c("*.*lf",10,4)
    fmt[3, ] <- c("*.*lf", 10, 4)
    if (all(ei.vc[Ghactual,]!=-1))
      se <- sqrt(matrix(c(extract.diag(vc), 0,0), ncol=1))
    else
      se <- matrix(NA, nrow=rows(b),ncol=1)
    
  ###  printfm(c(cml.ParNames,b,se),mask,fmt) ###COMPUTE
    print(c(cml.ParNames,b,se))
    cat("?\n")
    lst <- eirepar(b,Zb,Zw,x,evbase)
   
    bb <- lst[[1]]
    bw <- lst[[2]]
    sb <- lst[[3]]
    sw <- lst[[4]]
    rho <- lst[[5]]
    vrs <- as.matrix(c(bb, bw, sb, sw, rho))
    message("Reparameterization back to the truncated scale, parameterized",
    "        according to the underlying untruncated distribution.")
    print(as.vector(vrs))
    pb <- c(colMeans(as.matrix(bb)),colMeans(as.matrix(bw)),sb,sw,rho)
    print(pb)
    pb <- as.matrix(pb)
    if (Eprt>=3){
      cat("?\n")
      message("Reparameterization back to ultimate truncated scale")
      pb <- eirepart(b,Zb,Zw,x)  
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

##/* 
##    vc = gvc( &f, b, dataset );)
##**    
##** Various methods of computing variance matrix of parameters.
##** Calculates a "global variance matrix".
##**
##** These procedures generate a positive definite hessian more often
##** than hessp.src for non-normal likelihoods, while giving identical results
##** for normal p.d. likelihoods.  In theory at least they should also be more
##** accurate summaries of the likelihood surface than hessp, which gives
##** the quadradic approximation to the likelihood but only within epsilon
##** of the maximum.  If a p.d. hessian isn't found on the first round, 
##** additional tries are made with different values of the globals.
##** Allows for wide step lengths, quadratic approximation, gradient
##** checks, generalized inverses, and generalized choleskey's.
##**
##** Calls ghessp() to do most of the heavy lifting.
##** 
##** 
##** Inputs:     &f -- pointer to a log likelihood function f(b,dataset),
##**                   a procedure, taking a Kx1 vector and Nxk matrix as
##**                   arguments (f:Kx1,Nxk -> 1x1).
##**
##**           b -- Kx1 vector specifying the point at which the Hessian 
##**                 of f() is to be computed. (usually the maximum likelihood
##**                 estimates)
##**
##**           dataset -- a Nxk set of data, the same input as in CML and MAXLIK
##** 
##** Output:   vc -- KxK variance covariance matrix
##** 
##** GLOBALS:  _EI_vc = Mx2, M>=1: Each row represents globals for 1 instruction 
##**                 for an attempt to compute a positive definite variance 
##**                 matrix.  The procedure exits after the first p.d. matrix
##**                 is found.  Options to include in various rows:
##**                {1 0}=usual hessian computation, 
##**                {1 _GhDelta}=use usual hessian and then adjust eigenvalues
##**                   simultaneously so they are > _GhDelta. 
##**   		  {2 _GhFall}=use wide step lengths at _GhFall falloff in
##**	  	     the likelihood function
##**		  {3 _GhFall}=use quadradic approximation to wide step lengths 
##**                   with falloff in likelihood function set at _GhFall
##**                {4 0}=use generalized inverse and cholesky approach
##**                {5 _GhFall}=use wide step lengths, but compute the length by
##**                   checking the gradient at each point.  Start at _GhFall
##**                   falloff in the likelihood function.
##**                {-1 0}=do not compute the variance covariance matrix and return 
##**                   -hessian instead which will be used for the multivariate 
##**                   normal sampling with singular value decomposition.
##**               (If the likelihood is normal, i.e., quadratic, all options
##**               will give nearly the same answers.)
##**              DEFAULT={1 0, 2 0.1, 2 0.05, 3 0.1, 1 0.1, 1 0.2};
##**
##**           _GhFix = vector of same size as b with 1 for the estimated coeffs
##**                    and 0 for those fixed and to be deleted prior to gvc
##**                    computation; or set to 1 for keep them all (default)
##**                (presently only works if fixed coeff's are at end of vector)
##**
##** OUTPUT GLOBALS:
##**           _GhActual = row of _ei_vc at which pos.def. hessian was found
##**
##*/
gvc <- function(fn,b,dataset,evbase=parent.frame())
{ ##f is the pointer to a function
 ## local hessian,vc,m,i,q1,q2,mq,eigvals,not_pd,j;
  ##  clearg _gvc_dataset,_gvc_procname,_GhActual,_gvc_fixKeep,
  ## _GhFall,_GhQuad,GhDelta;
 
  evloc <- getEnvVar(evbase, environment())##
  gvc.dataset <- dataset;
  gvc.ProcName <- f <- fn;
  if (rows(GhFix)!=rows(b) || cols(GhFix)!=1)
    stop("gvc: _GhFix error")
  if (scalone(Ghfix) || Ghfix==1){
    gvc.fixKeep <- NA
  }else{
    print(b)
    gvc.fixKeep <-  subset(b, subset=!GhFix)
    b <- subset(b,subset=Ghfix)
  }
  gvc.procedure <- function(b, gvc.FixKeep,gvc.dataset){
    return(gvc.func(b, gvc.FixKeep,gvc.dataset))}
      
  for (i in 1:rows(ei.vc))
    {
      if(ei.vc[i+0,1]==1)
        {    ###              @ use hessp @
          if (Eprt>=2)
            message("gvc: trying numerical hessian")
          
          hessian <- optim(b,gvc.procedure,method="BFGS",hessian=TRUE)$hessian; ###gvc.procedure pointer to func 
          GhDelta <- ei.vc[i+0,2]
          if (GhDelta>0)
            {
              if(Eprt>=2)
                message("     with eigenvalue floor of ",GhDelta, "\n")
              lst <- eighv(-hessian) ###COMPUTE
              q1 <- lst[[1]]
              q2 <- lst[[2]]
              mq <- minc(q1)
              if (mq < GhDelta){
                q1 <- q1+0.1-minc(q1)
                q2 <- as.matrix(q2)
                hessian <- -q2 %*% diagrv(diag(rows(q1)),q1)%*%t(q2);
              }
            }### if (GhDelta>0)
      
        }else if(ei.vc[i+0,1]==2)
          {	###@ use wide step lengths @
            GhFall <- ei.vc[i+0,2]
            if (Eprt>=2)
              message("gvc: trying wide step lengths with likelihood falloff of ",ftos(GhFall,"*.*lf",1,2))
            GhQuad <- 0
            GhDelta <- 0
            assign("GhQuad", 0, env=evbase)
            assign("GhDelta", 0, env=evbase)
            
            hessian <- optim(b,gvc.procedure,method="BFGS",hessian=TRUE)$hessian 
      
          }else if (ei.vc[i+0,1]==3)
            { ###		@ use quadratic approximation @
              if (Eprt>=2)
                message("gvc: trying quadratic approximation with likelihood falloff of ", ftos(GhFall,"*.*lf",1,2))
              
              GhFall <- ei.vc[i+0,2]
              GhQuad <- 1
              GhDelta <- 0
              assign("GhQuad", 0, env=evbase)
              assign("GhDelta", 0, env=evbase)
              assign("GhFall", 0, env=evbase)
            
            
              hessian <- optim(b,gvc.procedure,method="BFGS",hessian=TRUE)$hessian 
      
            }else if(ei.vc[i+0,1]==4)
              { ###           @ use generalized approach @
                if (Eprt>=1)
                  message("gvc: trying generalized inverse of -hessian matrix")
                hessian <- optim(b,gvc.procedure,method="BFGS",hessian=TRUE)$hessian
                ##try using solve and if it does not work then svd.inv
                vc <- inv(-hessian) 
                if(Eprt>=3){
                  print("hessian:",hessian)
                  print("first vc:",vc)
                }
                eigvals <- svd(vc,nu=0,nv=0)$d
                not.pd <- 0
                for (j in 1:length(eigvals)){
                  if(Eprt>=3)
                    message("gvc: vc eigenvalue ", eigvals[j+0])
                  if (eigvals[j+0] <= 0) not.pd <- 1
                }
                if(not.pd == 1 && Eprt>=1)
                  message("gvc:-hessian matrix not p.d.; trying generalized cholesky also")
                vc <- try(chol(vc), silent=TRUE)
                if(inherits(vc, "try-error"))
                  vc <- sechol(vc) 
                vc <- t(vc)%*%vc
                if (Eprt>=3){
                  print("new vc:")
                  print(vc)
                }
              
              }else if (ei.vc[i+0,1]==5)
              {###  @ use wide step lengths with gradient check @
                GhFall <- ei.vc[i+0,2]
                if (Eprt>=2)
                  message("gvc: trying wide step lengths with gradient checks",
                          " and likelihood falloff of ",ftos(GhFall,"*.*lf",1,2))
                GhQuad <- 0
                GhDelta <- 0
                hessian <- optim(b, gvc.procedure,method="BFGS",hessian=TRUE)$hessian 
                assign("GhDelta", 0, env=evbase)
                assign("GhFall", GhFall, env=evbase)
                assign("GhQuad", 0, env=evbase)
                
              }else if (ei.vc[i+0,1] == -1)
                { ###      @ do not compute the vc matrix @
                  hessian <- optim(b,gvc.procedure,method="BFGS",hessian=TRUE)$hessian
                  if (Eprt>=2)
                    message("gvc: avoided computing the variance covariance matrix.")
                  GhActual <- i+0
                  assign("GhActual", i+0, env=evbase)
                  return(-hessian)   ###             @ return -hessian instead @
                }                
    
      if( ei.vc[i+0,1] != 4) ### @ did not use generalized approach @
        vc <- try(invpd(-hessian), silent=T)
        
      if(!inherits(vc, "try-error")){
        GhActual <- i+0
        assign("GhActual", GhActual, env=evbase)
        if (Eprt>=2)
          message("gvc: success.")
        return(vc)
      }
  
}##end for loop:for (i in 1:rows(ei.vc))
  
  message("gvc: positive definite Hessian not found.  Change ei.vc and rerun")
  GhActual <- -1; 
  return(NA); 
  
}

###dummy procedure for input to ghessp */

gvc.func <- function(b, gvc.FixKeep,gvc.dataset){
 
  b <- na.omit(rbind(b,gvc.FixKeep))
  loglik <- gvc.ProcName
  res <- loglik(b,gvc.dataset)
  res <- colSums(as.matrix(res))
  return(res)
}

##/*
##**  y = sechol(A);
##**  by Jeff Gill, April 2002.
##**
##**  This procedure produces:
##**
##**  y = chol(A+E), where E is a diagonal matrix with each element as smaii
##**  as possible, and A and E are the same size.  E diagonal values are 
##**  constrained by iteravely updated Gerschgorin bounds.  
##**
##**  REFERENCES:
##**
##**  Jeff Gill and Gary King. 2002. "What to do When Your Hessian is Not
##**  Invertible: Alternatives to Model Respecification in Nonlinear
##**  Estimation," http://gking.harvard.edu.
##**
##**  Robert B. Schnabel and Elizabeth Eskow. 1990. "A New Modified Cholesky
##**  Factorization," SIAM Journal of Scientific Statistical Computating,
##**  11, 6: 1136-58.
##
###   cvswork/ei/p
##*/
sechol <- function(A, macheps=2.23e-16)
  {
 
    if(!length(macheps)) macheps <- .Machine$double.eps
    n <- rows(A);
    m <- cols(A);
    L <- matrix(0,nrow=n,ncol=n)
    deltaprev <- 0
    gamm <- maxc(abs(extract.diag(A))) 
    tau <-  macheps^(1/3);
   if (min(eigen(A)$values) > 0)
    ###  /* print("not pivoting"); */
      tau <- -1000000
  
   if (m != n){
      message("sechol: input matrix must be square")
      return(A)
    }

   norm.A <- maxc(colSums(as.matrix(abs(A))))
   gamm <- maxc(abs(extract.diag(A))) 
   delta <- maxc(maxc(cbind(macheps*norm.A,macheps)))
   Pprod <- diag(n);
  
   if (n > 2){ 
      for (k in (1:(n-2))){
          if (A[1,1] < 0 || (minc(extract.diag(A[(k+1):n,(k+1):n]) - A[k,(k+1):n]^2/A[k,k]) < tau*gamm 
		&& minc(eigen(A[(k+1):n,(k+1):n])$values)< 0)){
                   
            dmax <- maxindc(extract.diag(A[k:n,k:n]))
            if (A[(k+dmax-1),(k+dmax-1)] > A[k,k]){
	       ###/* print("pivot using dmax:"); print(dmax); */
               P <- diag(n)
               Ptemp <- P[k,]
               P[k,] <- P[(k+dmax-1),]
               P[(k+dmax-1),] <- Ptemp
               A <- P%*%(A%*%P)
               L <- P%*%(L%*%P)
               Pprod <- P%*%Pprod
             }
            g <- matrix(0,nrow=n-(k-1),ncol=1)
            for (i in (k:n)){
               
               if(i == 1)
                 sum1 <- 0
               else
                 sum1 <- colSums(as.matrix(abs(A[i,k:(i-1)])))
	       
               if (i == n)
                 sum2 <- 0
               else
                 sum2 <- colSums(as.matrix(abs(A[(i+1):n,i])))
	       
               g[i-(k-1)] <- A[i,i] - sum1 - sum2
             }### end for (i in (k:n))
            
            gmax <- maxindc(g)
            if (gmax != k){
	       ###/* print("gerschgorin pivot on cycle:"); print(k); */
               P <- diag(n);
               Ptemp <- P[k,]
               P[k,] <- P[(k+dmax-1),]
               P[(k+dmax-1),] <- Ptemp
               A <- P%*%(A%*%P)
               L <- P%*%(L%*%P)
               Pprod <- P%*%Pprod
             }
            normj <- colSums(as.matrix(abs(A[(k+1):n,k])))
            delta <- maxc(cbind(0,deltaprev,(-A[k,k]+maxc(cbind(normj,tau*gamm)))))
            if (all(delta > 0)){
               A[k,k] <- A[k,k] %plus% delta
               deltaprev <- delta
             }
         } ### end of if(...&& minc(eigen(A[(k+1):n,(k+1):n])$values)) < 0 || A[1,1] < 0)
         
         A[k,k] <- sqrt(A[k,k])
         L[k,k] <- A[k,k]
        
        
         for (i in ((k+1):n)){ 
            if (L[k,k] > macheps) A[i,k] <- A[i,k]/L[k,k]; 
   	    if (is.infinite(A[i,k])) A[i,k] <- 0
            L[i,k] <- A[i,k]
            A[i,(k+1):i] <- A[i,(k+1):i] - L[i,k]%*%t(L[(k+1):i,k])
            if (A[i,i] < 0) A[i,i] <- 0; 
          }### end for(i in ((k+1):n)
        
        
       }###end for (k in (1:(n-2))){
    }##end if (n > 2)
    
   A[(n-1),n] <- A[n,(n-1)]
   eigvals <- eigen(A[(n-1):n,(n-1):n])$values
 
   pt1 <- - minc(eigvals)+ tau*maxc( rbind(1/(1-tau)*(maxc(eigvals)%minus% minc(eigvals)),gamm))
   
   dlist <- rbind(0,deltaprev,pt1)
 
   if (dlist[1] > dlist[2]) 
      delta <- dlist[1]   
   else
      delta <- dlist[2]
 
   if (delta < dlist[3])
      delta <- dlist[3]
  
   if (delta > 0){
      A[(n-1),(n-1)] <- A[(n-1),(n-1)] + delta
      A[n,n] <- A[n,n] + delta
      deltaprev <- delta
    }
   A[(n-1),(n-1)] <- sqrt(A[(n-1),(n-1)])
   L[(n-1),(n-1)] <- A[(n-1),(n-1)]
   A[n,(n-1)] <- A[n,(n-1)]/L[(n-1),(n-1)]
   if (is.infinite(A[n,(n-1)])) A[n,(n-1)] <- 0
   L[n,(n-1)] <- A[n,(n-1)]
   A[n,n] <- sqrt(A[n,n] - L[n,(n-1)]^2)
   L[n,n] <- A[n,n]
    res <- t(Pprod) %*% t(L) %*% t(Pprod)
   return(res)
  }

##/* infinity checking function from Gauss manual, page 72. */
isinf <- function(x){return(is.infinite(x))}
  
test.sechol <- function(){
  AA <- matrix(1:16, nrow=4, ncol=4)
  testAA <- matrix(c(1,0,0,0,2,1.4142,0,0,3,0.7071,1.2248,0,4,0,0,0.009), nrow=4)
###  print(testA)
  sAA <- sechol(AA)
  test <- abs(testAA - sAA)
 
###  print(test)
 
  if(any(test > 5*1.e-3))
    warning("sechol does not get right result")

  m <- matrix(c(5,1,1,3),2,2)
  cm <- chol(m)
  scm <- sechol(m)
  test <- abs(cm -scm)
  if(any(test > 5*1.e-3))
    warning("sechol does not get right result")
  
}
