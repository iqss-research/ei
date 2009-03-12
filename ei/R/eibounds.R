##/*
##**  This archive is part of the program EI
##**  (C) Copyright 1995-2001 Gary King
##**  All Rights Reserved.
##*/
##/*
##    {Bs,aggs}=bounds1(t,x,n);
##**
##** computes bounds on parameters given aggregate data
##**
##** INPUTS: nx1 vectors, unit of analysis is precinct
##** see output of sumvar()
##**
##** OUTPUTS: bounds on precinct-level parameters
##** Bs   = cols: lower-black ~ upper-black ~ lower-white ~ upper-white 
##** aggs =  bounds on district aggregates
##**         cols: lower ~ upper
##**         rows: beta-b, beta-w
##*/
### Translation of the Gauss code by Gary King
### AUTHOR: Ferdinand Alhimadi & Elena Villalon
###         evillalon@iq.harvard.edu
###
#include ei.ext;
bounds1<-function(t,x,n,tol,eps=.Machine$double.eps){
        ## local LbetaB,UbetaB,LbetaW,UbetaW,aggs,omx,Nb,Nw,c,c0,c1,p,tx,tomx,z,o,m;
        x[x<=0] <- eps
        x[x>=1] <- 1-eps
        omx<-1-x;c
        Nb<-x*n;                       #nr of black people
        Nw<-omx*n;                     #nr of white people
        ##  {c,c0,c1} = homoindx(x);       # proc in eiloglik.src
        cs<-homoindx(x,tol)
        c<-cs$c
        c0<-cs$c0
        c1<-cs$c1
        p<-rows(x);                    # nr of precinct could by nrow(data) if the args are dataset
        
        LbetaB<-matrix(0, nrow=p,ncol=1)
        UbetaB<-matrix(0,nrow=p,ncol=1)
        LbetaW<- matrix(0, nrow=p,ncol=1)
        UbetaW<-matrix(0, nrow=p,ncol=1)
        z<-matrix(0, nrow=p,ncol=1)
        o<-matrix(1, nrow=p,ncol=1)
        m<-o*NA;
        
        c <- na.omit(c)
        if (length(c)){
                tx<-t[c]/(x[c])
                tomx<-t[c]/(omx[c])
                LbetaB[c]<-maxr(z[c],tx-(omx[c]/(x[c])))
                UbetaB[c]<-minr(tx,o[c])
                LbetaW[c]<-maxr(z[c],tomx-(x[c]/(1+eps-x[c])))
                UbetaW[c]<-minr(tomx,o[c])
        }
        c0 <- na.omit(c0)
        if (length(c0)){## homogeneously white 
                LbetaB[c0]<-m[c0]
                UbetaB[c0]<-m[c0]
                LbetaW[c0]<-t[c0]
                UbetaW[c0]<-t[c0]
        }
        c1 <- na.omit(c1)
        if (length(c1)){## homogeneously black @
                LbetaB[c1]<-t[c1]
                UbetaB[c1]<-t[c1]
                LbetaW[c1]<-m[c1]
                UbetaW[c1]<-m[c1]
        }
        
### fix rounding errors due to machine precision */
### basically change any negative value to 0 and any value >1 to 1
        vu <- as.matrix(c(0,1))
        
        LbetaBa <- recode(LbetaB,cbind((LbetaB<=0),(LbetaB>=1)),vu)
    ###      LbetaBa <- recode(LbetaB,cbind((LbetaB<0),(LbetaB>1)),vu)
        UbetaB  <- recode(UbetaB,cbind((UbetaB<=0),(UbetaB>=1)),vu)
        LbetaW  <- recode(LbetaW,cbind((LbetaW<=0),(LbetaW>=1)),vu)
        UbetaW  <- recode(UbetaW,cbind((UbetaW<=0),(UbetaW>=1)),vu)
        
        
        bs<-list(bs=cbind(LbetaB,UbetaB,LbetaW,UbetaW))
       
        
        
        aggs<-list(aggs=rbind(cbind(meanwc(LbetaB,Nb), meanwc(UbetaB,Nb)),
                        cbind(meanwc(LbetaW,Nw),meanwc(UbetaW,Nw))))
        res <- c(bs, aggs)
        return(res)
}

##
##  {bnds,aggs} = bounds2(v,t,x,n);
##
##  INPUTS:
##  v = democratic fraction of the two party vote
##  t = fraction of people turning out to vote
##  x = fraction of people who are black
##  n = number of voting age people in each precinct
##
##  OUTPUTS: 
##  Bounds on fraction of blacks (lambdaB) and whites (lambdaW) voting for the dems
##  bnds = lower_lambdaB ~ upper_lambdaB ~ lower_lambdaW ~ upper_lambdaW
##  aggs =  bounds on district aggregates
##         cols: lower ~ upper
##         rows: lambdaB, lambdaW
##
bounds2 <- function(v,t,x,n,tol){

###  local LlambdaB,UlambdaB,LlambdaW,UlambdaW,aggs,omx,Nb,Nw,c,c0,c1,p,tx,
###  tomx,z,o,m,d;
  omx <- 1-x;
  Nb <- x*n;
  Nw <- omx*n;
  lst <- homoindx(x,tol);
  c <- lst$c
  c0 <- lst$c0
  c1 <- lst$c1
  p <- rows(x);
  
  LlambdaB <- matrix(0, nrow=p,ncol=1);
  UlambdaB <- matrix(0,nrow=p,ncol=1);
  LlambdaW <- matrix(0,nrow=p,ncol=1);
  UlambdaW <- matrix(0,nrow=p,ncol=1);
  z <- matrix(0,nrow=p,ncol=1);
  o <- matrix(1, nrow=p,ncol=1);
  m <- o*NA;
  
  if(!scalmiss(c)){###			@ heterogeneous precincts @
    d <- v[c]*t[c];
    LlambdaB[c] <- maxr(z[c],d-(1-x[c]))/(maxr(z[c],d-(1-x[c]))+minr(t[c]-d,x[c]));
    UlambdaB[c] <- minr(d,x[c])/(minr(d,x[c])+maxr(z[c],(t[c]-d)-(1-x[c])));
    LlambdaW[c] <- maxr(z[c],d-x[c])/(maxr(z[c],d-x[c])+minr(t[c]-d,1-x[c]));
    UlambdaW[c] <- minr(d,1-x[c])/(minr(d,1-x[c])+maxr(z[c],(t[c]-d)-x[c]));
    
   ### /* fix unanimous districts */
    LlambdaB[c] <- missrv(LlambdaB[c],v[c]);
    UlambdaB[c] <- missrv(UlambdaB[c],v[c]);
    LlambdaW[c] <- missrv(LlambdaW[c],v[c]);
    UlambdaW[c] <- missrv(UlambdaW[c],v[c]);
  }
   ### 
  if(!scalmiss(c0)){ ###			@ homogeneously white @
    LlambdaB[c0] <- m[c0];
    UlambdaB[c0] <- m[c0];
    LlambdaW[c0] <- v[c0];
    UlambdaW[c0] <- v[c0];
  }
    
  if(!scalmiss(c1)){###			@ homogeneously black @
    LlambdaB[c1] <- v[c1];
    UlambdaB[c1] <- v[c1];
    LlambdaW[c1] <- m[c1];
    UlambdaW[c1] <- m[c1];
  }

  ###/* fix rounding errors due to machine precision */
  LlambdaB <- recode(LlambdaB,cbind(LlambdaB<0,LlambdaB>1),as.matrix(0:1));
  UlambdaB <- recode(UlambdaB,cbind(UlambdaB<0, UlambdaB>1),as.matrix(0:1));
  LlambdaW <- recode(LlambdaW,cbind(LlambdaW<0, LlambdaW>1),as.matrix(0:1));
  UlambdaW <- recode(UlambdaW,cbind(UlambdaW<0, UlambdaW>1),as.matrix(0:1));
  
  aggs <- rbind(cbind(meanwc(LlambdaB,Nb),meanwc(UlambdaB,Nb)), 
                cbind(meanwc(LlambdaW,Nw),meanwc(UlambdaW,Nw)));
  lst <- c(list(cbind(LlambdaB,UlambdaB,LlambdaW,UlambdaW)), list(aggs=aggs))
 
}
### support proc
## row maximum
###
maxr<-function(a,b){
  a <- as.matrix(a)
  b <- as.matrix(b)
  res<-apply(cbind(a,b),1,max)
  
  return (as.matrix(res))
}

### support proc
## row minimum
###
minr<-function(a,b){
  a <- as.matrix(a)
  b <- as.matrix(b)
  res<-apply(cbind(a,b),1,min)
  return (as.matrix(res))
}

