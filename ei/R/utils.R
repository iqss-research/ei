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
  
       length(y) ==1  && !any(is.na(y)) && is.numeric(y) && all(y == 0)  
      
}

###
##    t = scalone(y);
##
##    y = matrix
##    t = 1 if y is a scalar one and 0 otw
##    as.numeric(TRUE) =1, as.numeric(FALSE) =0

scalone<-function(y){
 length(y) ==1  && !any(is.na(y))&& is.numeric(y)  && all(y == 1)  
       
}
### test if y is scalar missing value
scalmiss <- function(y){
  length(y) ==1 && any(is.na(y))
}

###
##  y = meanWc(x,wt);
##
##  x = NxM matrix to be avg'd
##  wt = scalar 1 or Nx1 or NxM weight used in averaging
##
##  works with missing values; packs rowwise.
## weighted.mean

meanwc<-function(x,wt){
     
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
  
  cv <- names(dbuf)
  res <- TRUE
  ix <- grep(str, cv, ignore.case=T)
  if(length(ix) <= 0)
    res <- FALSE

  return(res);
}

vread <- function(dbuf, str){
 
  cv <- names(dbuf)
  ix <- grep(str, cv, ignore.case=T)
  if(length(ix) <=0){
     warning(paste("Variable", str, "is not in the data buffer"))
     return(list())
   }
  return(dbuf[[ix]]) 
  
}

vget <- function(dbuf, str){
 
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

vput <- function(dbuf=list(), x, xname){
  if(!is.list(dbuf))
    stop("Ei: vput the input buffer needs to be a list")
  nm <- NULL
  nc <- length(dbuf)
  if(length(dbuf))
    nm <- names(dbuf)
  dbuf[[nc+1]] <- x
  names(dbuf) <- c(nm, xname)
  
return(dbuf)
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
  x <- matrix(x)
  v <- matrix(v)
  dm <- dim(e)
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
