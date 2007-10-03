
ei <- function(t,x,tvap,Zb, Zw,...)
{
 ###  local res,et,MLpsi,MLvc,betaBs,betaWs,tst,Eselect,flat;
  evbase <- eiset(t,x,tvap,Zb,Zw,...)  ##environment 
  drvdot <- match.call(expand.dots=TRUE)
  drv  <-  match.call(expand.dots=FALSE)
  n <- tvap
 
  evbase <- expanddots(drvdot,drv,evbase)
 
  ###	@ timing start @
  et   <- proc.time()
  param  <- ls(env=evbase)

  if(length(grep("Eprt", param)) && length(grep("Eversion", param))){
    Eversion <- get("Eversion", env=evbase)
    message(Eversion)
  }
 
  ### copy variables from evbase to local environment evei
  evei <- getEnvVar(evbase, environment())  ##environment 
 
  if(!scalzero(Eres) && vin(Eres,"titl") && Eprt >0)
    message(vread(Eres,"titl"))
 
  Eres <- add.to.Eres(Eres, round=1, evbase)
 
###testing and debugging            
###  lapply(param, function(x){
###    print(x)
###    print(get(x, env=evbase))})
###    checkinputs(t,x,n,Zb,Zw);
###   return(evbase)
###  /* verify inputs */
  if(as.logical(Echeck)){
 
    tst <- ""
    tst <- checkinputs(t,x,tvap,Zb, Zw,evbase)
    if(tst != ""){
      print(tst)
  ###    "----- EI Aborted -----";
      return(Eres)
    }
    if(Eprt>0){
      if(EnonPar)
        print("Inputs ok, beginning nonparametric estimation...")
      else
      print("Inputs ok, beginning preliminary estimation...")
    }
  }
  
 ###  /* augment _Eselect if _EselRnd<1 */
  assign("Eselect", Eselect, env=evbase);  ###$@ save existing value @

  Eselect <- matrix(1,nrow=x,1)* as.vector(Eselect);
  
  if(EselRnd<1){
    vec <- runif(rows(x),min=0,max=1)
    vec <- matrix(vec)
    Eselect <- Eselect & (vec < EselRnd);
  }


###  /* nonparametric estimation */
  ###if(EnonPar){ ### check remove
    if(T){
       
    betaBs <- einonp(t,x, evbase);
    if(length(betaBs))
      Eres <- vput(Eres,betaBs,"betabs");
    Eres <- add.to.Eres(Eres, round=2, evbase)  
    return(timing(et,Eprt,Eres))
  }
  ### /* parametric estimation: */
  
  ###/* eta influence on Zb,Zw */
  if(rows(Eeta)!=4){
  if(Eeta[1]==1 ||  Eeta[1]==4){
    Zb <- x;
    Zw <- 1;
  }else if(Eeta[1]==2 || Eeta[1]==5){
    Zb <- 1;
    Zw <- x;
  }else if(Eeta[1]==3){
    Zb <- Zw <- x;
 
  }
}

  ### /* set internal global */

  assign("Ez", 0, env=evbase)
 ### clearg _Ez;	@ n of covariates, incl. implied constant term for Zb|Zw @
  Ez <- (cols(Zb)+1-scalone(Zb))|(cols(Zw)+1-scalone(Zw));


###  /* likelihood estimation */
  if(EdoML==1){
   
     lst  <- quadcml(x,Zb,Zw,T);
     MLpsi <- lst$Mlpsi
     MLvc <- lst$MLvc
     
     if(is.na(MLvc)){
       Eres <- vput(Eres,MLpsi,"phi");
       Eres <- vput(Eres,MLvc,"vcphi");
      return(Eres);
     }
  }else{
    message("Skipping likelihood estimation..");
     Eres <- add.to.Eres(Eres, round=3, evbase)  
  
    MLpsi <- mlpsi <- EdoML.phi;
    MLvc <-  mlvc <- EdoML.vcphi;
    if(rows(mlvc)!= rows(mlpsi))
      stop("ei: EdoML.phi or EdoML.vcphi input error");
    
  }

 ###/* simulation */
  if(EdoSim==1){
  ###  {betaBs,betaWs} = psim1(T,X,tvap,Zb,Zw,MLpsi,MLvc);
     lst <- psim1(T,X,tvap,Zb,Zw,MLpsi,MLvc);
     betaBs <- lst$betaBs
     betaWs <- lst$betaWs
     Eres   <- vput(Eres,betaBs,"betaBs"); ###@ no need to save betaWs; see eiread @
   }


 return(timing(et,Eprt,Eres));
 

}

          
### timing:
timing <- function(et, Eprt,Eres){
  if (Eprt>0){
    et <- proc.time() -et;###@ timing end  @
    et <- na.omit(et)
    vec <- c("User", "System", "Total1","Total2", "CumChild")
    ln <- length(et)
    vec <- vec[1:ln]
    
 ###   fmtt()
    et <- unlist(sapply(1:length(vec), function(n) paste(vec[n], et[n], "   ")))
    message("Done. Time =", et);
  ###  "----- EI Completed -----";
  }

tst <- paste("Run time: ", date(), "\nEversion", sep="")
  Eres <- vput(Eres,tst,"date");
  res <- Eres;
  Eres <- vput(list(),tst,"date");

##  assign("Eselect", Eselect, env=evbase)
 
  return(res)
}
             

checkinputs <- function(t,x,n,Zb, Zw,evbase=NULL){
 
###getting all global parameters from env=evbase
  if(!length(evbase))
    evbase <- get("evbase", env=parent.frame())
  ### assign the globals to local environment 
  evei <- getEnvVar(evbase, environment())
###  t <- as.matrix(t)
  t <- as.matrix(vread(Eres, "t"))
 
###  x <- as.matrix(x)
  x <- as.matrix(vread(Eres, "x"))

###  tvap <- as.matrix(n)
  tvap <- as.matrix(vread(Eres, "n"))
  
###  Zb <- as.matrix(Zb)
    Zb <- as.matrix(vread(Eres,"Zb"))
  
###  Zw <- as.matrix(Zw)
  Zw <- as.matrix(vread(Eres,"Zw"))

  zb <- Zb
  zw <- Zw
  ei.vc <- EI.vc 
  
  EIgraph.bvsmth <- eigraph.bvsmth
  
  ei2.m <- Ei2.m 
  
####
  if(any(t <0) || any(t > 1))
    stop("ei:'t' input must be between 0 and 1")
  if(any(x <0) || any(x > 1))
    stop("ei:'=x' input must be between 0 and 1")
  if(length(EnonPar) > 1)
    stop("ei:EnonPar must be scalar")
  if(!EnonPar %in% 0:1)
    stop("ei: EnonPar must be 0 or 1")
  r <- rows(t)
  if(!length(r)) r <- length(t)
 
  if( r != rows(x) || ( r != rows(Zb) && !scalone(Zb)) || (r != rows(Zw) && !scalone(Zw)))
    stop("ei: inputs do not have the rigth dimension")
  if(EnonPar == 0){
    if(any(tvap <= 0))
      stop("ei:'n' input must be greater than zero");
    
    if( r!=rows(tvap))
      stop("ei:'n' input does not have the right dimension");
  }
  if (any(is.na(x)) || any(is.na(t)) || any(is.na(tvap)) || any(is.na(zb)) || any(is.na(zw)))
    stop("ei:missing data detected. Delete and rerun");
  
  if(length(Esims)>1 || Esims < 1)
    stop("ei: Esims must be a scalar integer larger than 1");
  
  if(length(Eselect) > 1 ){
    
    if(!any(dim(Eselect)== c(r, 1)) || !all(as.vector(Eselect) %in% 0:1))
      stop("ei: Eselect must be scalar or px1 vector of 0's and 1's");
  }else if(!scalone(as.vector(Eselect)))
    stop("ei: Eselect must be px1 or scalar 1");
  
  
  if(length(EselRnd)>1)
    stop("ei: EselRnd must be a scalar");
  if (EselRnd> 1 || EselRnd <=0)
    stop("ei: EselRnd must in interval (0,1]");

  Eeta <- as.matrix(Eeta)
### Eeta ony one column and possibly multiple rows
  if(cols(Eeta) != 1)
    stop("ei: Eeta has the wrong dimensions");
### 4 rows
  if(rows(Eeta) == 4){
    if(!all(Eeta[3:4,]>=0))
      stop("ei: Eeta[3:4,] must be >=0");
  }
### either 1 (scalar) or 3 rows
  else if (!rows(Eeta) %in% c(1,3)){
    stop("ei: Eeta has the wrong dimensions");
###scalar, possible values are 0:3
  }else if(rows(Eeta) ==1){
    if(!Eeta %in% 0:3)
      stop("ei: Eeta is wrong")
    if(Eeta !=0 && !scalone(Zb) && !scalone(Zw))
      stop("ei: Eeta=1,2,or 3, Zb and Zw must be set to 1")
### 3 rows: restrictions
  }else if(Eeta[1] == 4 && !scalone(Zw))
    stop("ei: if Eeta[1]=4, Zw must be set to 1")
  else if(Eeta[1]==5 && !scalone(Zb))
    stop("ei: if Eeta[1]=5, Zb must be set to 1")

  if (rows(Eeta)==4 || all(Eeta==0)){
    a <- 5+cols(zb)+cols(zw)-scalone(zb)-scalone(zw);
  }else if(Eeta[1] %in% c(1,2, 4,5)){
    a <- 6;
  }else if(Eeta[1]==3){
    a <- 7;
  }
  
### done with Eeta

###EalphaB 
  if(!scalmiss(EalphaB)){
    if(scalone(Zb) && all(Eeta ==0))
      stop("ei:EalphaB should be specified only when Zb is not 1")
    if(cols(EalphaB) !=2)
      stop("ei: EalphaB must be missing or have 2 columns")
   
    if(any(EalphaB[, 2] <=0))
      stop("ei: Elements in the second column of _EalphaB must be > 0");
    
    if(rows(EalphaB) != cols(Zb) && all(Eetha ==0))
      stop("ei: nrow(EalphaB) must equal ncol(Zb)");
    
  }else{
    if(!scalone(Zb))
      warning("Including covariates (Zb) without priors works but is not generally recommended.  See _EalphaB")
  }
  
  
###EalphaW
  if(!scalmiss(EalphaW)){
    if(scalone(Zw) && all(Eeta==0))
      stop("ei: EalphaW should be specified only when Zw is");
    if(cols(EalphaW) != 2)
      stop("ei: EalphaW must be missing or have 2 columns");
     
    if(any(EalphaW[,2] <= 0))
      stop("ei: Elements in the second column of EalphaW must be > 0");
    
    if(rows(EalphaW)!=cols(Zw) && all(Eeta==0))
      stop("ei: nrow(EalphaW) must equal ncol(Zw)");
  }else{
      if(!scalone(Zw))
        warning("Including covariates (Zw) without priors works but is not generally recommended. See EalphaW");
    }

###Ebeta
    if(length(Ebeta) > 1)
      stop("ei: _Ebeta must be a scalar");
    if(Ebeta < 0)
      stop("ei: Ebeta cannot be negative");  
    
###Ecdfbvn
    
    if(any(Ecdfbvn > 6) || any(Ecdfbvn <1) || rows(Ecdfbvn) != 1)
      stop("ei: problem with Ecdfbvn");
###Esigma

    if(length(Esigma)>1)
      stop("ei: Esigma must be a scalar");
    if(Esigma > 0 && Esigma <0.000001)
      stop("ei: Esigma must be <= 0 (for no prior) or > 0.000001");
    if( (Erho[1]==0 && rows(Erho)!=2) || (Erho[1]!=0 && rows(Erho)!=1))
      stop("ei: problem with Erho");
    
###Estval  
    if(cols(Estval) !=1 )
      stop("ei: Estval may have only one column");
    if(length(Estval)==1){
### Estval needs to check
###      if(Estval != 0 && (Estval<0 || (Estval%/% 3==0) || ((Estval %%floor(Estval) )> 0)))
      if(Estval <0 || Estval == 2 || (Estval %%floor(Estval))> 0)
        stop("ei: Estval as a scalar must be 0 or an integer >=3");
###length(Estval) > 1
    }else if(rows(Estval) !=a)
      stop("ei: _Estval has wrong dimensions");
    
###Ebounds
    if(length(Ebounds) == 1){
      if(!Ebounds %in% c(0, 1))
        stop("ei: Ebounds must be 0, 1, 1x2, or kx2")
    }else{
      if(cols(Ebounds) != 2)
        stop("ei: Ebounds must have 1 or two columns")
      if(rows(Ebounds) != 1 && rows(Ebounds) != a)
        stop("ei: Ebounds must have 1 row or one row for each parameter") 
    }
###EdirTol
    if(length(EdirTol) > 1)
      stop("ei: EdirTol must be a scalar")
    if(EdirTol <= 0 || EdirTol > 1)
      stop("ei: EdirTol must be >0 and <1")

###EcdfTol
    if(length(EcdfTol) > 1)
      stop("ei: EcdfTol must be a scalar")
    if(EcdfTol <= 0 || EcdfTol > 1)
      stop("ei: EcdfTol must be >0 and <1")
    
### EvTol
    if(length(EvTol) > 1)
      stop("ei: EvTol must be a scalar")
    if(EvTol <= 0)
      stop("ei: EvTol must be >0")
###ei.vc
    if(cols(ei.vc) != 2)
      stop("ei: EI.vc must have 2 columns")
    if (min(ei.vc[,1]) < -1 || max(ei.vc[,1]) > 5 || abs(ei.vc %%floor(ei.vc)) >0)
      stop("ei: EI.vc may only have integers -1,1,...,5 in first column");
    if (min(ei.vc[,1])==-1 &&  eist !=0)
      stop("ei: EI.vc={-1 0} option is allowed only when EIsT = 0"); 
###EdoML
    if(!(as.vector(EdoML) %in% 0:1))
      stop("ei: EdoML must be scalar zero or one");
    
    if(scalzero(EdoML)){
      if( (rows(EdoML.phi) != rows(EdoML.vcphi)) ||
         (rows(EdoML.phi)!= cols(EdoML.vcphi)) )
        stop("ei: EdoML, EdoML.phi, or EdoML.vcphi are incorrect")
    }
###Edosim
    if(!(EdoSim %in% -1:1))
      stop("ei: EdoSim must be scalar -1, 0, or 1");
    
###eigraph.bvsmth
    if(length(eigraph.bvsmth) > 1)
      stop("ei: EIgraph.bvsmth must be a scalar");
    if (EIgraph.bvsmth<0.0000001)
      stop("ei: EIgraph.bvsmth must be greater than zero");
###Eisn
    if(length(Eisn) >1 || Eisn <1)
      stop("ei: Eisn must be a scalar or an integer >= 1");
###Eist
    if(length(Eist) >1)
      stop("ei: Eist must be a scalar")
    if(Eist <= 2.00001 && Eist != 0)
      stop("ei: EisT must be 0 or greater than 2");
    
###EisFac
    if(length(EisFac) > 1)
      stop("ei: EisFac must be scalar")
    if ( EisFac<0 && !EisFac %in% c(-1,-2))
      stop("ei: EisFac must be -1, -2 or greater than zero");
###EisChk
    if(length(EisChk) > 1 || !EisChk %in% c(0, 1))
      stop("ei: EisChk must be scalar zero or one");
    
### EmaxIter
    if (length(EmaxIter) >1 || EmaxIter <0)
      stop("ei: EmaxIter must be a scalar > 0")
###EnumTol
    if(length(EnumTol) > 1 || EnumTol <0)
      stop("ei: EnumTol must be a scalar >0");
###Eprt
    if (length(Eprt) >1 )
      stop("ei: Eprt must be a scalar");
    if(!Eprt %in% c(0,1,2,3))
      stop("ei: Eprt must be 0, 1, 2, or 3");
    
###EnonEval
    if(length(EnonEval) >1 || EnonEval<1)
      stop("ei: EnonEval must be a positive integer");
    
###EnonNumInt
    if(length(EnonNumInt) >1 || EnonNumInt < 1)
      stop("ei: EnonNumInt must be a positive integer");
    if(length(Ei2.m) >1 )
      stop("ei: Ei2.m must be positive");

    if((Ei2.m <1 || Ei2.m %% floor(Ei2.m) > 0) &&   Ei2.m != -1)
      stop("ei: Ei2.m must be -1 or a positive integer");
    
###EIMetaR
    if(length(EIMetaR) >1 ||  EIMetaR %% floor(EIMetaR) > 0)
      stop("ei: EiMetaR must be a scalar integer");
    
    if(!vin(Eres, "truth")){
    ###  message("No truth")
      str <- " "
      return(str)
    }
      ## eiread not written yet but returns vread     
  betaB <- eiread(Eres, "truthB")
  betaW <- eiread(Eres, "truthW")
  if(rows(betaB) != rows(x) || rows(betaW) != rows(x))
    stop("ei: stored 'truth' must have same dimensions as x & t")
  if (any(betaB >1)  || any(betaB <0))
    stop("ei: 'truthB' input must be between 0 and 1");
  if (any(betaW >1)  || any(betaW <0))
    stop("ei: 'truthW' input must be between 0 and 1");

  res <- bounds1(t, x, tvap)
  bnd <- res$bs
  a <- res$aggs
  a <- na.omit(cbind(betaB, bnd[,1]))
  tol <- Enumtol
  
  if(any((a[,1]+tol) < a[, 2])){
    evei <- settTruthBnds(betaB, betaW, bnd, truthB, truthW, bnds,  evei)
    str <- "ei: truthB < lower bound";
    return(str);
  }
    
        
  a <- na.omit(cbind(betaB, bnd[,2]));
  if(any((a[,1]-tol)> a[,2])){
    evei <- settTruthBnds(betaB, betaW, bnd, truthB, truthW, bnds,  evei)
    str <- "ei: truthB > upper bound"
    return(str)
  }
  a <- na.omit(cbind(betaW, bnd[,3]));
  if(any((a[,1]+tol) < a[,2])){
    evei <- settTruthBnds(betaB, betaW, bnd, truthB, truthW, bnds,  evei)
    str <- "ei: truthW < lower bound"
    return(str)
  }
      
  a <- na.omit(cbind(betaW, bnd[,4]));
  if(any((a[,1]-tol) >a[,2])){
    evei <- settTruthBnds(betaB, betaW, bnd, truthB, truthW, bnds,  evei)
    str <- "ei: truthW > upper bound";
    return(str)
  }

               
  str <- " "
  return(str)
}
### DESCRIPTION set the globals parameters and if something is passed with
###             the dots it will update the default values in eiset
###             As in the Gary King Gauss code
###
### INPUT the same as in ei() call
### OUTPUT the envoronment containing all globals parameters
###
###
eiset <- function(t,x,tvap,Zb,Zw,...){
  ## general
  
  
  Eversion="EI Version: 1.9, 2/8/2003";
##  Eres=vput("","Run time: "$+datestr(0)$+" "$+timestr(0)$+", "$+Eversion,
##              "date");
  Eres <- list()
  driver <- match.call()
  args <- names(driver)
    

  Eres <- c(Eres,list(eiversion=Eversion))
    Echeck <- as.matrix(0);###Echeck <- as.matrix(1); check remove
  Esims <- as.matrix(100);
  Eprt <- as.matrix(2);
  Eselect <- as.matrix(1);
  EselRnd <- as.matrix(1);
  EdoSim <- as.matrix(1);
  EdoML <- as.matrix(1);
  EdoML.phi <- as.matrix(0);
  EdoML.vcphi <- as.matrix(0);
  Ecdfbvn<- as.matrix(5);
  EnumTol<- as.matrix(0.0001);
  EnonEval<- as.matrix(11);
  EnonNumInt<- as.matrix(11);
  EnonPar<-  as.matrix(0);
  Ei2.m<- as.matrix(-1);
  eimetar<- EIMetaR <- as.matrix(1);
  ei2.mta<- as.matrix(0);
  
  ## priors 
  Erho<- as.matrix(0.5);
  Esigma<- as.matrix(0.5);
  Ebeta<- as.matrix(0);
  EalphaB<- as.matrix(NA);
  EalphaW<- as.matrix(NA);

  ## quadcml\
  ######check this procedure.......???????
 ### cmlset;
###  Estval<- as.matrix(1);
  Estval <- as.matrix(1)
  Ebounds<- as.matrix(1);
  Eeta<- as.matrix(0);
  EdirTol<- as.matrix(0.0001);
  EcdfTol<- as.matrix(3e-15);
  EvTol<- as.matrix(5e-307);
  Ez<- matrix(1,nrow=2, ncol=1);  ## internal, same as rbind(1,1) 
  EmaxIter<- as.matrix(500);

  ## gvc 
  EI.vc <-  ei.vc <- c(1, 0, 4, 0, 2, 0.1, 2, 0.05, 3, 0.1, 1, 0.1, 1, 0.2);
  EI.vc <- matrix(EI.vc, ncol=2,byrow=T)
  GhQuad<- as.matrix(0);
  GhFall<- as.matrix(0.1);
  GhDelta<- as.matrix(0);
  GhStart<- as.matrix(0);
  GhFix<- as.matrix(1);
  GhActual<- as.matrix(0);
  gvc.dataset<- as.matrix(0);
  gvc.procname<- as.matrix(0);
  gvc.fixKeep<- as.matrix(0);
  
  ## psim1 
  EbetaWs<- as.matrix(0);
  
  ## rndisamp 
  Eist<- as.matrix(0);
  Eisn<- as.matrix(10);
  EisFac<- as.matrix(4);
  EisChk<- as.matrix(0);
  EiLliks<- as.matrix(0);

  ## lnpdfmn 
  Eivc<- matrix(NA);

  ## internal global 
  EisFix<- matrix(NA);

  ## graphics Check procewdures?????????/
  ###graphset;
  ###graphgk;  

  ## eigraph 
  eigraphC<- 1;
  eigraph.thick<- 1;
  eigraph.Xlo<- 0;
  eigraph.Xhi<- 1;
  eigraph.Tlo<- 0;
  eigraph.Thi<- 1;
  eigraph.bblo<- 0;
  eigraph.bbhi<- 1;
  eigraph.bwlo<- 0;
  eigraph.bwhi<- 1;
  eigraph.x<- "X";
  eigraph.t<- "T";
  eigraph.bb<- "betaB";
  eigraph.bw<- "betaW";
  eigraph.loess<- 0;
  eigraph.eval<- 31;
  eigraph.bvsmth<- Eigraph.bvsmth <- 0.08;
  eigraph.smpl<- 1;
  eigraph.dbuf<- 0;
  tomogPct<- matrix(c(.5, .95));
  tomogClr<-matrix(c( 12, 9, 10, 11, 13, 5)); 
  eigraph.pro<- matrix(NA);
  
  ## eibias 
  Evc<- as.matrix(1);
  Etruth<- as.matrix(0);
  
  ## eicond 
  eicond.nums<- matrix(c( -1.5, 1.5, 50 ), ncol=1);
  
 ## internal globals 
  eigraph.circ <-  as.matrix(0);
  eigraph.psiu <-  as.matrix(0);
  
  ## dens 
  smth <-  as.matrix(0.03);
  strt <-  as.matrix(0);
  endd <-  as.matrix(0);
  pts  <-  as.matrix(100);
  Tleft <-  as.matrix(0);
  Tright <-  as.matrix(0);
  kern <-  as.matrix("E");
  whiskr <-  as.matrix(-1);
  jitter <-  as.matrix(0);
  output <-  as.matrix(1);
  
  ## loess 
  loess.Span <-  as.matrix(.6667);
  loess.NumEval <-  as.matrix(50);
  loess.Degree <-  as.matrix(1);
  loess.WgtType <-  as.matrix(2);
  
  ## reg 
  Rweight <-  as.matrix(1);
  Rxnames <-  as.matrix(0);
  Ryname  <-  as.matrix(0);
  Rfmt    <-  as.matrix(4);
  Rselect <-  as.matrix(1);
  Rrobust <-  as.matrix(-1);
  Rtheta  <-  as.matrix(0);
  Routput <-  as.matrix(1);
  Rconst  <-  as.matrix(1);
 
  ## token2 
  tokdel<- as.matrix(c(32,10,13,44,9))
  tokwds<- as.matrix(-1);

  ## eimodels 
  EImodels.save<- ""; 
  EI.bma.prior<- as.matrix(0);
  EI.bma.est<- as.matrix(1);

  ### gauss is case independent
  param <- ls(env=environment())
  paramlower <- sapply(param,tolower)
  ix <- 1:length(paramlower)
  ### assign the values in param to the lower case names
  evloc <- environment()
  res <- sapply(ix,function(n, ev =evloc){
    if(identical(paramlower[n], param[n]))
      return(NULL)
    assign(paramlower[n], get(param[n], env=evloc),env=ev)
    if(n >= length(ix))
      return(ev)
  })
  ### lower case only the first letter
  evloc <- unlist(res)[[1]]
  paraml <- sapply(param, function(m){
    mm <- strsplit(m, NULL)[[1]]
    mm[1] <- tolower(mm[1])
    return(paste(mm,collapse=""))})
  
  ix <- 1:length(paraml)
  res <- sapply(ix,function(n, ev =evloc){
    nml <- strsplit(paraml[[n]], NULL)[[1]]
    nm  <- strsplit(param[[n]],NULL)[[1]]
    if(identical(nml[1], nm[1]))
      return(NULL)
    assign(paraml[n], get(param[n], env=evloc),env=ev)
    if(n >= length(ix))
      return(ev)})
  res <- unlist(res)[[1]]
  return(as.environment(res))

}

