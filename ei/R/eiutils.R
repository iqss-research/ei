## Some utilities that Gary King wrote in Gauss for library ei translated
## into R code and some other utilities that Elena wrote in R that has
## has some corresponding functions in Gauss
## Defined in original Gauss code are:
## scalzero,scalone, meanwc (meanWc), ismissm, seqas, seqase, vin,
## sortbyRow, mkmissm, fisherzi,fisherz, ftosm, fmtt, infrv  
## 
###
##  t = scalzero(y);
##
##  y = matrix
##  t = T if y is 0 and F otherwise
##

scalzero<-function(y){
    y <- as.vector(y)
    length(y) ==1  && !any(is.na(y)) && is.numeric(y) && all(y == 0)  
      
}

###
##    t = scalone(y);
##
##    y = matrix
##    t = 1 if y is a scalar one and 0 otw
##    as.numeric(TRUE) =1, as.numeric(FALSE) =0

scalone<-function(y){
  y <- as.vector(y)
  
 length(y) ==1  && !any(is.na(y))&& is.numeric(y)  && all(y == 1)  
       
}


### Wrapper around the R code weighted.mean 

meanwc<-meanWc <- function(x,wt,na.remove=TRUE){
  x <- as.matrix(x)
 
  wwt <- wt
  if(scalmiss(wt) || wt==1)
    wwt <- rep(1, rows(x))
  wwt[is.na(wwt)]<- 0
  ln <- -floor(-length(x)/length(wwt))
  
  if(ln >1) 
    wwt <- rep(wwt, ln)
  if(length(wwt) > length(x))
    wwt <- wwt[1:length(x)]
  wt <- as.matrix(wwt)
  dm <- dim(x)
 
  if(dm[2]==1)
    return(weighted.mean(x,wt,na.rm=na.remove))
  return(apply(x,2,weighted.mean,,na.rm=na.remove))
}
              
    
###
##  y = meanWc(x,wt);
##
##  x = NxM matrix to be avg'd
##  wt = scalar 1 or Nx1 or NxM weight used in averaging
##
##  works with missing values; packs rowwise.
## weighted.mean

meanwcFerdi<-meanWcFerdi <- function(x,wt){
     
     ###  lapply(as.data.frame(x), weighted.mean, wt)
        x <- as.matrix(x)
        if(all(is.na(wt)) || wt==1)
          wt<-rep(1,nrow(x))
        wwt<-wt
        wwt[is.na(wwt)]<-0
        xtmp<-x
        xtmp[is.na(xtmp)]<-0
        if(is.matrix(x))
          res<-(apply((xtmp*wwt),2,sum))/(apply((1-ismissm(x+wt))*wwt,2,sum))
        else
          res<-sum(xtmp*wwt)/(sum(1-ismissm(x+wt))*wwt)
        return(res)
}

##/*
##** sims = rndchi(r,c,v);
##**
##** inputs: r = row
##**         c = column
##**         v = df
##**
##** output: sim = kxk matrix of independent chi-square simulations with v
##**
##** 4/13/99 KS
##*/
rndchi <- function(r,c,v){
    return(2*rndgam(r,c,v/2))
  }



###
##  y = ismissm(x)
##
##      x = an n x m matrix
##
##      y = an n x m matrix of 0's (indicating not missing) and 1's (missing)
##

ismissm<-function(x){
       is.na(x)
}

### same as seqase but not including the end-points
seqas <- function(strt,endd,pts){
  res <- (endd-strt)/pts
  fst <- strt+0.5*res
  inc <- res
  n <- pts
  res <- seq(from=fst, by=res, length.out=pts)  
  
  return(res)
}
seqase<-function(strt,endd,pts){
  return(seq(from=strt,to=endd,length.out=pts))}

### dbuf is a named list 
### str is string.
### Check if str is part of the
### names of dbuf independently of case.

vin <- function(dbuf,str){
  str <- tolower(str)
  str <- paste("^",str,"$",sep="")
  cv <- names(dbuf)
  cv <- unlist(sapply(cv,tolower))
  res <- TRUE
  ix <- grep(str, cv, ignore.case=T)
  if(length(ix) <= 0)
    res <- FALSE

  return(res);
}




###DESCRIPTION Takes a matrix and sort all rows independently
###            from each other or only the rows in index ix
###
sortbyRow <- function(mat, ix=NULL){
  mad <- as.data.frame(t(mat))
  mm <- mat
  ind <- 1:length(mad)
  if(length(ix))
    ind <- ix
  for( n in ind){
    v <- sort(mad[[n]])
    mm[n, ] <- v
  }
  
  return(mm)
}



##    y = mkmissm(x,m);
##
##  x = n x k data matrix
##  m = n x k matrix of 0's and 1's
##
##  y[i,j] = missing if m[i,j]=1; else y[i,j]=x[i,j]
##
##  example:
##
##  m=ismissm(d);       /* remember which are missing */
##  d2=missrv(d,-9)     /* change missing to -9's */
##  /* do recodes to d2 as if no missing values were present */
##  d3=mkmissm(data2,m)  /* after recodes, return to missing*/
##
## From Gary's code

mkmissm <- function(x, m){
  mn <- as.numeric(m)
  if( !all(as.vector(mn) %in% 0:1))
    warning("mkmissm: m must have only 0 or 1 or T, F entries")
  m <- as.logical(m)
  y[m] <- NA
  return(y)

}

##  This archive is part of the program EI
##  (C) Copyright 1995-2001 Gary King
##  All Rights Reserved.
##
##  Utility procs
##
##   z = fisherzi(x);
##   inverse of fisher's z transformation
##
fisherzi <- function(x){
 
  t <- exp(2*x);
  t <- (t-1)%dot/%(t+1);
 
  return(t);
}
##
##   z = fisherz(x);

##   fisher's z transformation
##
fisherz <- function(x){
  t=0.5*log((1+x)%dot/%(1-x));
  return(t);
}

                      

###/* ftosm - matrix field to string
##**
##** y = ftosm(sym,num,num2);
###**
##**  input:  sym = a character vector
##**          num = number of characters to use
##**
##**  output: y = a character vector of numbers
###*/
ftosm <- function(sym,num,num2){

 
  res <- matrix(0, nrow=rows(sym),ncol=1)
  if(num2==0)num2 <- 1
  for (i in 1:rows(sym))
    res[i+0] <- ftos(sym[i+0],"d",digits=abs(num),width=num2)
 
  return(res);
  }

###digits total digist displayed including the 
fmtt <- function(x){
return(formatC(x, digits=7, width=3))
}
                 

##/* reverse infinities
##
## y = infrv(x,minus,plus);
## x = input vector
## minus, plus = scalars
## y = an ExE conformable matrix with -INF changed to minus and +INF changed
##       to plus
##
infrv <- function(x,m,p){
  x <- as.matrix(x)
  plus <- Inf
  minus <- -Inf;
  s <- seq(from=1,by=1, length.out=nrow(x))
 

  pinf <- subset(s,subset= !(x != plus))  ###=(.not(x ./= plus)));
  minf <- subset(s,subset=!(x != minus));
  pinf <- as.matrix(pinf)
  minf <- as.matrix(minf)
  res <- x
  
  if(!scalmiss(pinf))
    res[pinf] <- matrix(p,nrow=rows(pinf),ncol=1);
  
  
  if(!scalmiss(minf))
    res[minf] <- matrix(m,nrow=rows(as.matrix(minf)),ncol=1)
  
  return(res)
}

 


###/*
##   y = makefacn(vars,nums);
##**
##**  vars  = number of columns to make in y
##** nums   = (levels x 1) vector of numbers to use in place of 1,2,3... in
##**          makefac().
##**       OR (levels x vars) matrix
##**
##** example:
##**
##**  nums=(1~ 2 ~3)|
##**       (9~10~11);
##**  call makefacn(3,nums);
##**   
##**   1    2    3 
##**   9    2    3 
##**   1   10    3 
##**   9   10    3 
##**   1    2   11 
##**   9    2   11 
##**   1   10   11 
##**   9   10   11 
##**
##*/
makefacn <- function(vars,nums){
###  local x,k,i,c;
  c <- cols(nums);
  k <- rows(nums);
  
  if (c==1){
    x <- makefac(vars,k)
    for (i in (1:vars))
      x[,i+0] <- nums[x[,i+0]];
  
  }
  else if( c==vars){
    x <- makefac(vars,k)
    for (i in (1:vars))
      x[,i+0]=nums[x[,i+0],i+0];
    
  
  }else
    stop("makefacn: input error, cols(nums) must = 1 or vars")
    
  return(x)
}


##/*
##       y = makefac(vars,levels);
##**
##**  vars = number of columns to make in y
##** levels= number of levels of y
##**
##** y = (levels^vars x vars).  first column is 1,2,3,...,vars,
##**     second column is 1,1,1,1(vars times),2,2,2,, etc
##**    third column...
##**   
##** example:
##** y = makefac(2,4);
##**   y;
##**  1   1
##**  2   1
##**  3   1
##**  4   1
##**  1   2
##**  2   2
##**  3   2
##**  4   2
##**  1   3
##**  2   3
##**  3   3
##**  4   3
##**  1   4
##**  2   4
##**  3   4
##**  4   4
##**
##** NOTE: round() corrects for numerical inaccuracies.
##**       i+0 works around a Gauss bug
##*/

makefac <- function(vars,levels)
  {
    ## local i,tmp;
    if (length(vars) > 1 || length(levels) > 1)
      stop("makefac: arguments must be scalars");
   
 
    tmp <- matrix(1:levels,nrow=round(levels^vars),ncol=1)
  
    for (i in 1:vars){   
      ff <- as.matrix(facvec(i+0,vars,levels))
      tmp <- cbind(tmp,ff)
      
    }
    return(tmp[,2:cols(tmp)]);
  }
    
##/*  facvec: support proc for makefac
##**  
##   y = facvec(i,v,l);
##**
##** i = var number
##** v = value
##** l = level
##**
##** y = output column vector
##**
##** round() corrects for numerical inaccuracies.
##*/
##  v1 <- matrix(seq(from=1, by=1,length.out=l),nrow=1)
##  v1 <- t(v1 %dot*% matrix(1,nrow=round(l^(i-1)),ncol=l))
##  v1 <- v1 %dot*% matrix(1, nrow=1, ncol=round(l^(v-i)))

facvec <- function(i,v,l){
  if(i > v)
    stop("Bad arguments in facvec")
  mat <- matrix(0,nrow=l^v)
  if(i <= 1)
    return(vec <- as.matrix(rep(1:l,l^v/l)))
  vec <- NULL
  for(n in 1:l)
    vec <- c(vec, rep(n,l^(i-1)))
  
  if(length(vec) < l^v)
    vec <- rep(vec,l^v/length(vec))
  
  return(as.matrix(vec))
  
}

