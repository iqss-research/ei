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
    xx <- as.data.frame(x)
  for(n in 1:length(v)){
    if(n > length(xx)) break; 
    ind <- is.na(xx[[n]])
    x[ind, n] <- v[n]
  }
 
  return(x)
  }
