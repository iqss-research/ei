###
##  Gauss functions and their corresponding R counterparts
##



###
##    m <- zeros(r,c)
##    r - nr of rows
##    c - nr of columns
##    m - matrix of zeros
##
zeros<-function(r,c){
      matrix(0, nrow=r, ncol=c)
}


###
##    m <- ones(r,c)
##    r - nr of rows
##    c - nr of columns
##    m - matrix of ones
## 

ones<-function(r,c){
      matrix(1, nrow=r, ncol=c)   
}

###
##    res <- menac(x)
##    x   - matrix
##    res - vector with mean of each col. in matrix x 
## 
meanc<-function(x, na.rm=F){
###       mean(as.data.frame(x))

   return(colMeans(x,na.rm=na.rm))
 

}




###DESCRIPTION calculate standard deviations along the columns of
###            a matrix x, same as the Gauss function
###
stdc <- function(x){
  x <- as.data.frame(x)
  return(sd(x))
}

###
## counts the number of elemetns of a vector that fall into a specific rage
## c = count(x,v)
## x = Nx1 vector containing the numbers to be counting
## v = Px1 vector (sorted in asc order) containing the ranges within which counts are to be made

counts<-function(x,v){
       subx <- x[x >= v[1] & x <= v[length(v)]]
       return(length(subx))
       
}


## remove t first and b last rows from x

trimr<-function(x,t,b){
        
        return(x[(t+1):(nrow(x)-b),])
}

###
## removes the last row of the original matrix AND add a
## row with "missing values" (as the first row of the original matrix)

lag<-function(x){
       x <- as.matrix(x)
       x <- rbind(NA, x) 
       x[-nrow(x),]
        
}

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
### test if y is scalar missing value
scalmiss <- function(y){
  length(y) ==1 && any(is.na(y))
}

###DESCRIPTION sort a matrix rows according to a column
###            as the corresponding Gauss
###
sortc <- function(mat, c=1, decreasing=FALSE){
  mat <- as.matrix(mat)
  ord <- order(mat[,c],na.last=TRUE)
  mat <- mat[ord,]
  return(mat)
}
###
##  y = meanWc(x,wt);
##
##  x = NxM matrix to be avg'd
##  wt = scalar 1 or Nx1 or NxM weight used in averaging
##
##  works with missing values; packs rowwise.
## weighted.mean

meanwc<-meanWc <- function(x,wt){
     
     ###  lapply(as.data.frame(x), weighted.mean, wt)

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

ismiss <- function(x){
  any(is.na(x))
}
###
##  create vector of PTS evenly spaced points between STRT and ENDD,
##  including the end points.

seqase<-function(strt,endd,pts){

        t<-(endd-strt)/(pts-1)
        res<-seq(strt,endd,t)
        return (res)
}

### dbuf is a named list 
### str is string.
### Check if str is part of the
### names of dbuf independently of case.

vin <- function(dbuf,str){
  str <- paste("^",str,"$",sep="")
  cv <- names(dbuf)
  res <- TRUE
  ix <- grep(str, cv, ignore.case=T)
  if(length(ix) <= 0)
    res <- FALSE

  return(res);
}

vread <- function(dbuf, str){
 
  cv <- names(dbuf)
  str <- paste("^",str,"$",sep="")
  ix <- grep(str, cv, ignore.case=T)
  if(length(ix) <=0){
     warning(paste("Variable", str, "is not in the data buffer"))
     return(list())
   }
  return(dbuf[[ix]]) 
  
}
### similar to vread but returns also the list dbufnew without
### the element str 
vget <- function(dbuf, str){
   str <- paste("^",str,"$",sep="")
  cv <- names(dbuf)
  ix <- grep(str, cv, ignore.case=T)
  if(!is.list(dbuf) || length(ix) <=0){
     warning(paste("ei: vget. Variable", str, "is not in data buffer"))
     return(list())
   }
  var <- dbuf[[ix]]
  dbufnew <- dbuf[-ix]
  lst <- c(list(var), dbufnew)
  return(lst)
}  
### Similar to the Gauss vput: It inserts
### a new element in the list dbuf and
### return the new list

vput <- function(dbuf=list(), x, xname=NULL){
  if(!is.list(dbuf))
    stop("Ei: vput the input buffer needs to be a list")
  nm <- NULL
  nc <- length(dbuf)
  if(length(dbuf))
    nm <- names(dbuf)
  dbuf[[nc+1]] <- x
  ### trying to figure out names
  xname1 <- NULL
  if(length(names(x)))
    xname1 <- names(x)
  else  if(length(colnames(x)))
    xname1 <- colnames(x)
  if(!length(xname))
    xname <- xname1 
  if(length(xname))
    names(dbuf) <- c(nm, xname)
  
  
return(dbuf)
}
vnamecv <- function(dbuf){
  return(names(dbuf))
}
### DESCRIPTION the gauss code is case insensitive and that it was
###             inG is doing, which is a wrapper around %in% but ignoring case
### flag =0, 1, 2 for character case sensitive, numeric , character case insensitive
### as in(y,cv,0)
###
"%inG%" <- function(y,vars){
  y <- unique.default(y)
  vars <- unique.default(vars)
  retp <- FALSE
  ###make it case insensitive
  if(!is.numeric(y)) 
    y <- toupper(y)
  if(!is.numeric(vars))
    vars <- toupper(vars)
    
  ##  mat <- setdiff(y, vars)
  ##  if(length(mat))
  ##    retp <- TRUE
  retp <- all(y %in% vars)
  return(retp)                 
    
}
### DESCRIPTION similar to the recode func in Gauss
###             x vector of N values
###             e matrix NxK of 1's and 0's
###             v vector of K elements
### OUTPUT y is Nx1 matrix containing recoded values of x
###        See recode.test() for an example
### AUTHOR Elena Villalon
### Date August 17th, 2007
###
recode <- function(x,e,v){
  x <- as.matrix(x)
  e <- as.matrix(e)
  dm <- dim(e)
 
  if(length(v) <= 1 && dm[[2]]> 1)
    v <- rep(v,dm[[2]])
  v <- as.matrix(v)
  if(dm[1] != nrow(x) || nrow(v) != dm[2] || !all(e %in% c(0,1)))
    stop("recode: check your inputs")
  y <- e %*% v
  ix <-  y != 0
  x[ix, ] <- y[ix, ]
  return(x)
  }
###Another version less efficient
recode1 <- function(x,e,v){
  x <- matrix(x)
  v <- matrix(v)
  dm <- dim(e)
  if(dm[1] != nrow(x) || nrow(v) != dm[2] || !all(e %in% c(0,1)))
    stop("recode: check your inputs")
  
  for(r in 1:nrow(e)){
  
   ev <- e[r, ]
   if(sum(ev) <=0) next;
   ev <-  as.numeric(ev)
   ix <- match(1, ev)
   x[r,] <- v[ix, ]
 }
  return(x)
  } 
  

  

recode.test <- function(){
   x <- matrix(c(20, 45,32,63,29))
   e1 <- 20 < x & x < 30
   e2 <- 30 < x & x < 40
   e3 <- 40 < x & x < 50
   e4 <- 50 < x & x < 60
   e <- cbind(e1, e2, e3, e4)
  v <- matrix(1:4)
   recode(x,e,v)
   
 }
selif <- function(x, e){
  subset(x, subset=e)
}

delif <- function(x, e){
  subset(x, subset=!e)
}
test.selif <- function(){
  x <- matrix(c(0, 30, 60, 10, 40, 70, 20, 50, 80), nrow=3)
  e <- (x[, 1] > 0 & x[,3] < 100)
}


sumc <- function(x){ return(colSums(x))}

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
### DESCRIPTION finds the index (row number) of the smallest element
###             in each column of a matrix as the Gauss function
minindc <- function(mat){
  mad <- as.data.frame(mat)
  ind <- 1:length(mad)
  res <- ind
  for(n in ind){
   v <- mad[[n]]
   mn <- min(v)
   ix <- grep(mn, v)
   if(length(ix) > 1)
     warning("minindc finds multiple values of min in col ", n)
   
   res[n] <- ix[1]
 }
  return(res)
}
### DESCRIPTION Correspond to Gauss function
###             For each column of x finds the values equal to the
###             corresponding v entry and turn it into NA
###
miss <- function(x, v){
  x <- as.matrix(x)
  
  if(length(v) <= 1) {
    v <- paste("^",v, "$", sep="")
    ind <- grep(v, x)
    if(length(ind))
      x[ind] <- NA
    return(x)
  }
  v <- sapply(v, function(m) paste("^", m, "$", sep=""))
  xx <- as.data.frame(x)
  for(n in 1:length(v)){
    if(n > length(xx)) break; 
    ind <- grep(v[[n]], xx[[n]])
    x[ind, n] <- NA
  }
 
  return(x)
}
### DESCRIPTION Correspond to Gauss function
###             For each column of x finds the values that are NA
###             and substitute them with the corresponding v entry.
###
missrv <- function(x, v){
  x <- as.matrix(x)
  if(length(v) <= 1){
   ind <- is.na(x)
   if(length(ind))
     x[ind] <- v
   return(x)
 }
  xx <- as.data.frame(x)
  
  for(n in 1:length(v)){
    if(n > length(xx)) break; 
    ind <- is.na(xx[[n]])
    if(any(ind))
      x[ind, n] <- v[[n]]
  }
 
  return(x)
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
  if( !all(as.vector(m) %in% 0:1))
    warning("mkmissm: m must have only 0 or 1 or T, F entries")
  y <- miss(m, 1) + x
  return(y)

}
###DESCRIPTION the function diag in Gauss returns a 1 column vector
###            with the diagonal values of mat, but diag in R builds
###            a diagonal matrix.  To avoid confusion we call extract.diag
###            to the R function that extracts the diagonal elemnts of a matrix
###            and returns a matrix with one column. 
extract.diag <- function(mat){
  v <- mat[col(mat)==row(mat)]
  return(as.matrix(v))
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
  t <- (t-1)/(t+1);
  return(t);
}
##
##   z = fisherz(x);
##   fisher's z transformation
##
fisherz <- function(x){
  t=0.5*log((1+x)/(1-x));
  return(t);
}
###DESCRIPTION x is a matrix compute variance covariance along cols
vcx <- function(x){return(var(x))}
cdfni <- function(x){qnorm(x)} ###inverse cumulative density  
cdfn <- function(x){pnorm(x)} ### cumulative density 
corrx <- function(x){cor(x)} ###correlation matrix
rndn <- function(r, c){
  return(matrix(rnorm(r*c, mean=0, sd=1), nrow=r, ncol=c, byrow=T))}
rndu <- function(r, c){
  matrix(runif(r*c), nrow=r, ncol=c, byrow=TRUE)}

                       
###DESCRIPTION If v is scalar finds the indices of elements in x == v
###            If v is length two find the indices of elements in x
###            such that x > v[1] & x <= v[2]
###            From Gaus utility 
indexcat <- function(x, v){
  if(length(v) <= 1)
    return(grep(v, x))
  v <- as.vector(sort(v))
  x <- as.vector(x)
  xx <- x[x> v[1] & x <= v[2]]
  xx <- unique.default(xx)
  ind <- unlist(sapply(xx, grep, x))
  return(ind)
}
dotfeq <- function(x,y, tol=NULL){
  x <- as.matrix(x)
  y <- as.matrix(y)
  if(!all(dim(x) == dim(y))){
    message("x and y have different dimensions")
    return(FALSE)
  }
  if(!length(tol))
    return(x==y)
  x <- floor(x/tol)
  y <- floor(y/tol)
  return(x==y)
}
  
###DESCRIPTION Computes the cdf of the standardized bivariate normal
###            with lower limits in -Inf, i.e. lower tail. 
###            x and t are the upper limits for the two variables
###            and rho is the correlation coefficients
###            Wraps pmvnorm of of package mvtnorm
###            
cdfbvn <- function(x,t,rho, maxpts=25000, abseps=0.001, releps=0){
  if(!require(mvtnorm))
    stop("ei:To compute bivariate normal you need to install package mvtnorm")
  v  <- c(as.vector(x), as.vector(t))
  ln <- length(v)
  low <- rep(-Inf, ln)
  p00 <- pmvnorm(lower=low, upper=v,mean=rep(0, ln), sigma=rho, maxpts=maxpts,abseps=abseps, releps=releps);
  return(p00)
 }
ftos <-  function(x){
  nc  <- nchar(as.character(floor(x))) 
      fmt <- formatC(x,width=nc,digits=0, format="f")
}
loess <- function(depvar, indvars,data, loess.span, loess.wgtType){
  y.loess <- loess(depvar~indvars, data, weights=loess.wgType, span=loess.span)
  yhat <- y.loess$fitted
  ys <- y.loess$y
  xs <- y.loess$x
  lst <- c(list(yhat=yhat), list(ys=ys), list(xs=xs))
  return(lst)
}
