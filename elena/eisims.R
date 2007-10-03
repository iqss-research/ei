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
