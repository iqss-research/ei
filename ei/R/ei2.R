##/* -----------------------------------------------------------------------
##   dbuf2 = ei2(V,dbuf,Gb,Gw);
##**
##** Ecological Inference Estimates for 2xC tables via multiple imputation.
##** Gives observation-level estimates of lambdaB and lambdaW given T & X in
##**     V = lambdaB.*x   + lambdaW.*(1-x), 
##** where x is estimated from a first stage EI() run from
##**     T = betaB.*X + betaW.*(1-X), and using x = HbetaB.*X./T 
##**     where HbetaB are multiply imputed estimates.
##**
##** INPUTS:
##** V    = outcome variable (p x 1) (Democratic fraction of vote)
##** dbuf = data buffer from first stage EI run (with inputs ei(t,x,n,Zb,Zw))
##** Gb   = 1 (for no covariates)
##**        or matrix of covariates to predict lambdaB (constant term is implied)
##** Gw   = 1 (for no covariates)
##**        or matrix of covariates to predict lambdaW (constant term is implied)
##**
##** OUTPUT:
##**  dbuf2 = packed data buffer. 
##**       Read contents with eiread(); graph with eigraph().
##**       Because these programs work for BOTH ei() and ei2(), reading
##**       estimates from dbuf2 requires using the names "betaB" and "betaW",
##**       even tho they might really be called "lambdaB" and "lambdaW".
##**       In addition, in EIREAD,
##**       _t        p x 1: variable that was originally t in dbuf
##**       _x        p x 1: variable that was originally x in dbuf
##**       _n        p x 1: variable that was originally n in dbuf 
##**       t         p x 1: redefined as V (dem proportion of two party vote)
##**       x2        p x _Esims: simulations of x from prior stage analysis
##**       x2rn      p x _Esims: x2 randomly horizontally permuted 
##**       x         p x 1: redefined as mean posterior
##**       Nb2       p x _Esims: denominator of x & t; x.*n (voting age blacks)
##**       Nw2       p x _Esims: (1-x).*n = n-Nb (voting age whites)
##**       see additional documentation under eiread and eigraph
##**        
##** GLOBALS:  all the globals from EI() work here too.  In addition,
##** _EI2_m   = -1, use the posterior mean of betabs to impute one data set.  
##**       or set to number of data sets to multiply impute (must be <=_Esims) 
##**       4 or so should be plenty
##**
##** OUTPUT GLOBAL:
##** _ei2_mta = a "meta-data buffer" with elements data buffers named dbuf#, 
##**            where #=1,2,...,_EI2_m.  Each is the output from a 
##**	      separate imputation run.  If ei2() fails, _ei2_mta will
##**            contain results from the imputation runs that completed.
##*/ 
ei2 <- function(V,dbuf,Gb,Gw){

  nm <- names(dbuf)
  if("evbase" %in% nm) evbase <- dbuf[["evbase"]]
  else evbase <- eiset()
  Eres <- underEres <- get("Eres", env=evbase)
  tst <- paste("Run time:", date())
  Eres <- vput(Eres,tst,"date")
  res <- Eres
 
  v <- V
  if (all(is.na(v))){
    message("ei2: 'V' should have no missing values")
    return(NA)
  }
 
  ei2.mta <- vput(list(),x="*MDB* Meta-Data Buffer from 2nd Stage *MDB*",xname="titl")
  ei2.m <- get("ei2.m", env=evbase)
  Esims <- get("Esims", env=evbase)
  underEsims <- Esims
  if (as.vector(ei2.m) != -1)
    Esims <- floor(as.vector(Esims)/as.vector(ei2.m))
 
 
 ### /* get data from first stage analysis data buffer */
  if (vin(dbuf,"x2"))
    x <- eiread(dbuf,"x2rn")
  else
    x <- eiread(dbuf,"x")
  
  underx <- x
  t <- eiread(dbuf,"t")
  undert <- t
  Enumtol <- EnumTol <- as.vector(get("EnumTol", env=evbase))
  t <- recode(t,as.vector(t)<=Enumtol,Enumtol/2)

  n <- undern <- eiread(dbuf,"n")
  betabs <- eiread(dbuf,"rnbetabs")
  mpx <- missrv(eiread(dbuf,"betab"),0)%dot*%colMeans(t(x),na.rm=FALSE)%dot/%t
  if(as.vector(ei2.m) == -1){### /* if _ei2_m = -1, use the posterior mean */
    x2 <- mpx
    ei2m <- 1
  }else{
    x2 <- (missrv(betabs[,1:ei2.m],0)%dot*%x)%dot/%t
    ei2m <- as.vector(ei2.m)
  }

  nt <- eiread(dbuf,"nt")
  nt <- recode(nt,as.vector(nt)<1,1)
  lambdaBs <- t*0
  x2k <- lambdaBs
  
 ### /* run second stage analyses */
  for (i in (1:ei2m)){
    if ((i+0)==ei2m){
      underEres <- Eres
      underEsims <- Esims-underEsims*(ei2m-1)
    }
    Eprt <- get("Eprt", env=evbase)
    if(Eprt>0)
      message("?;==================> Starting ei2 imputation ",(i+0),";?;")
    e <- cbind((x2[,i+0]>1),(x2[,i+0]<0))
    x2[,i+0] <- recode(x2[,i+0],e,as.matrix(c(1,0)))
    idbuf <- ei(V,x2[,i+0],nt,Gb,Gw)
    if (scalmiss(idbuf) || length(idbuf) <= 0 || !(vin(idbuf,"betabs"))){
      message("ei2: problem with EI run (noted above), terminating EI2.")
      return(NA)
    }
    
    ei2.mta <- vput(ei2.mta,x=idbuf,xname=paste("dbuf#",i+0,sep=""))
    lambdaBs <- cbind(lambdaBs,eiread(idbuf,"betabs"))
    x2k <- cbind(x2k,x2[,i+0]);
    Estval <- trimr(eiread(idbuf,"phi"),0,2)
  }
  
###/* prepare output */
###vput(list(),x,xname)
  idbuf <- vput(idbuf,lambdaBs[,(2:cols(lambdaBs))],"betaBs")
  idbuf <- vput(idbuf,x2k[,(2:cols(x2k))],"x2")
  mpx <- recode(mpx,cbind((as.vector(mpx)>1),(as.vector(mpx)<0)),rbind(1,0))
  idbuf <- vput(idbuf,mpx,"x")
  idbuf <- vput(idbuf,Esims,"underEsims")
  idbuf <- vput(idbuf,undert,"undert")
  idbuf <- vput(idbuf,underx,"underx")
  idbuf <- vput(idbuf,undern,"undernn")
  
###  "===== EI2 Estimation Complete =====";
  rm(Eres,env=evbase)
  return(idbuf)
}

