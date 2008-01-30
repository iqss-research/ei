###
##  Gauss functions and their corresponding R counterparts
##  Contains: 
##  zeros, ones, meanc, stdc, counts, trimr,lag, lag1, lagn
##  scalmiss,sortind,sortc,sortr,ismiss,seqa,vread, vget
##  vput, vnamecv, recode,selif, delif, sumc,minindc,miss,
##  missrv, vcx, corrx, indexcat,dotfeq,ftos,loess,strput,
##  rows, cols, vec, intquad1,substute,cumsumc,     
##  diagrv, eye, eig, inv, invpd, pinv, eig0, eighv, reshape

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
### DESCRIPTION: Returns vector with max along each column
###
maxc <- function(mat){
  mat <- as.matrix(mat)
  if(ncol(mat) <= 1)
    return(max(mat))
  
  return(apply(mat, 2, max))}
  
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
  v <- sort(v)
  vst <- v[1]
  res <- v
  if(length(v) <= 1)
    return(length(x[x <= v[1]]))
  res[1] <- length(x[x <= v[1]])
  
  for(n in 2:length(v)){     
    subx <- x[x > vst & x <= v[n]]
    res[n] <- length(subx)
    vst <- v[n]
  }
  return(res)
       
}

counts.test <- function(x=1:9, v=c(4,5,8)){
  res <- counts(x,v)
}

## remove t first and b last rows from x

trimr<-function(x,t,b){
  x <- as.matrix(x)
        
  return(x[(t+1):(rows(x)-b),])
}

###
## removes the last row of the original matrix AND add a
## row with "missing values" (as the first row of the original matrix)

lag1 <- lag<-function(x){
       x <- as.matrix(x)
       x <- rbind(NA, x) 
       x[-rows(x),]
        
}
###Same as lag specifying the number of rows to lag
lagn <- function(x, p){
  x <- as.matrix(x)
  r <- nrow(x)
  if(p >= r)
    return(x*NA)
  vec <- c(rep(1, r-p), rep(NA, p))
  res <- na.omit(x*vec)
  res <- rbind(matrix(rep(NA,p)*ncol(x),ncol=ncol(x), nrow=p),res)
  return(res)
}

### test if y is scalar missing value
scalmiss <- function(y){
  length(y) ==1 && any(is.na(y))
}
sortind <- function(x){
  x <- as.vector(x)
  return(as.matrix(order(x)))
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
sortr <- function(mat, r=1, decreasing=FALSE){
  return(t(sortc(t(mat),c=r,decreasing=decreasing)))

}

ismiss <- function(x){
  any(is.na(x))
}


##Additive sequance first point=st, increment =inc, number of points=n

seqa <- function(st,inc,n){
  return(seq(from=st, by=inc,length.out=n))
}


### dbuf is a named list 
### str is string.
### Check if str is in dbuf 
### If that is the case returns the list element
### associated with names str.
###
vread <- function(dbuf, str){
 
  cv <- names(dbuf)
  cv <- unlist(sapply(cv,tolower))
  str <- tolower(str)
  str <- paste("^",str,"$",sep="")
  ix <- grep(str, cv, ignore.case=T)
  if(length(ix) <=0){
     warning(paste("Variable", str, "is not in the data buffer"))
     return(NA)
   }
  res <- as.matrix(dbuf[[ix]])
  return(res) 
  
}
### similar to vread but returns also the list dbufnew without
### the element str 
vget <- function(dbuf, str){
   str <- paste("^",str,"$",sep="")
  cv <- names(dbuf)
  ix <- grep(str, cv, ignore.case=T)
  if(!is.list(dbuf) || length(ix) <=0){
     warning(paste("ei: vget. Variable", str, "is not in data buffer"))
     return(NA)
   }
  var <- as.matrix(dbuff[[ix]])
   
  dbufnew <- dbuf[-ix]
  lst <- c(list(var), dbufnew)
  return(lst)
}  
### Similar to the Gauss vput: It inserts
### a new element, x,in the list dbuf with list name = xname
### and return the new list with element inserted at end. 

vput <- function(dbuf=list(), x, xname=""){
  if(!is.list(dbuf))
    stop("Ei: vput the input buffer needs to be a list")
  nm <- NULL
  nc <- length(dbuf)
  if(length(dbuf))
    nm <- names(dbuf)
  
  dbuf[[nc+1]] <- x
  ### trying to figure out names
  xname1 <- xname
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
###Similar to gauss returns the names of list dbuf
vnamecv <- function(dbuf, trim=T){
  res <- names(dbuf)
  if(trim)
    res <- unlist(sapply(res, trim.blanks))
  return(res)
}

###DESCRIPTION same as the Gauss code returns the names of
###            objects stored in list dbuf and a description of them
vlist <- function(dbuf){
  obj <- vnamecv(dbuf)
  objdesc <- lapply(dbuf, mode)
  objdim <- lapply(dbuf, function(m){
    if(rows(m) > 1 || cols(m) >1)
      return(v <- paste("rows=", rows(m)," cols=", cols(m)))
    else
      return("Scalar")
  })
                
                               
  lst <- lapply(1:length(obj), function(n){
    nm <- obj[[n]]
    desc <- objdesc[[n]]
    dm <- objdim[[n]]
    return(c(nm,desc,dm))})
  names(lst) <- unlist(obj)

  return(lst)
}
  
### DESCRIPTION: Same as Gauss recode (see manual)
### Changes the value of existing matrix x elements to the new values
### that are stored in vector v
### x Nx1 matrix to recode
### e NxK matrix of 0 (FALSE) and 1 (TRUE)
### v Kx 1 matrix with the new values
### OUTPUT Nx1 vector with the recoded values of x
###
### AUTHOR: Elena Villalon (evillalon@iq.harvard.edu)
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
  
  for(r in 1:rows(e)){
  
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
  res <- subset(x, subset=e)
  if(!length(res)) res <- NA
  return(as.matrix(res))
}

delif <- function(x, e){
  res <- subset(x, subset=!e)
  if(!length(res)) res <- NA
  return(as.matrix(res))
}
test.selif <- function(){
  x <- matrix(c(0, 30, 60, 10, 40, 70, 20, 50, 80), nrow=3)
  e <- (x[, 1] > 0 & x[,3] < 100)
}


sumc <- function(x){ return(colSums(as.matrix(x)))}

### DESCRIPTION minimum along the columns of mat
minc <- function(mat) {
  if(ncol(as.matrix(mat)) <=1)
    return(min(mat))
  md <- as.data.frame(mat)
  return(sapply(md, min))
}

### DESCRIPTION finds the index (row number) of the smallest element
###             in each column of a matrix as the Gauss function
minindc.old <- function(mat){
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
minindc <- function(mat){
  return(maxindc(mat,decrease=FALSE))}

maxindc <- function(mat, decrease=TRUE){
  mm <- as.data.frame(mat)
  ix <- apply(mat,2, function(x){
    ord <- order(x, decreasing=decrease)
    return(ord[1])})
  return(ix)
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
    if(length(ind))
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


###DESCRIPTION x is a matrix compute variance covariance along cols
vcx <- function(x){return(var(x))}
corrx <- function(x){cor(x)} ###correlation matrix


###DESCRIPTION: log(mat!) 
   lnfact <- function(mat){ return(lfactorial(mat))}
    
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
  ln <- length(x)==1 || length(y)==1
  if(!all(dim(x) == dim(y)) && !ln)
    stop("x and y have different dimensions")
  if(length(x) ==1) x <- as.vector(x)

  if(length(y) ==1) y <- as.vector(y)
  
  if(!length(tol))
    return(x==y)
  x <- floor(x/tol)
  y <- floor(y/tol)
  return(x==y)
}
  

### digits how many digits total including the left and rigth of decimal point
### width= the total field width. If it is smaller than digits + decimal point
### then it will default to width=digits+1 or width=digits (if not decimal point)
### if width > 0 then right-justified otherwise left-justified.  

ftos <-  function(x,fmat="f", digits=NULL, width=1){
  nc  <- nchar(as.character(floor(x)))
  nm <- nchar(x%%floor(x))
  dc <- ifelse(nm > 5, 3, nm-2)
   
  if(!length(digits))
      return(fmt <- formatC(x,digits=(nc+dc), format="f"))
  return(fmt <- formatC(x,width=digits,digits=digits,format=fmat))
}
     

                
loess <- function(depvar, indvars,data, loess.span, loess.wgtType, deg){
  y.loess <- loess(depvar~indvars, data, weights=loess.wgType, span=loess.span,degree=deg )
  yhat <- y.loess$fitted
  ys <- y.loess$y
  xs <- y.loess$x
  lst <- c(list(fitted=yhat), list(ys=ys), list(xs=xs))
  return(lst)
}
strput <- function(substr, str,off){
  if(off > 1) stp <- off-1
  strp <- substr(str, 1, stp)
  return(paste(strp,substr, sep=""))
}


rows <- function(mat){
  mat <- as.matrix(mat)
  return(nrow(mat))
}
 

cols <- function(mat){
  mat <- as.matrix(mat)
  return(ncol(mat))
}

vec <- function(mat){
  v <- as.matrix(as.vector(mat))
  return(v)
}

int <- function(x){ return(floor(x))}
    
  
stof <- function(x){ return(as.matrix(as.numeric(strsplit(x, " ")[[1]])))}

###DESCRIPTION: Integral of function f with limits in vector v
  intquad1 <- function(f,v){
    lst <- integrate(f, lower=v[2], upper=v[1])
    return(lst$y)
  }

##DESCRIPTION as in the Gauss function. Substitutes old values
##            for new values in matrix x, according to  the
##            outcome of the logical expression e and values in v.
###
## Elena Villalon (evillalon@iq.harvard.edu)
###
substute <- function(x, e, v){
  x <- as.matrix(x)
  ee <- e <- as.matrix(e)
  e <- as.logical(e)
  if(!is.matrix(e))
    e <- matrix(e,nrow=rows(ee),ncol=cols(ee))
  v <- as.matrix(v)
  ix <- grep(TRUE, e)
  if(length(v)< length(x)){
    mat <- matrix(FALSE, nrow=rows(x), ncol=cols(x))
    mat[e] <- v
  }else
  mat <- v
 
    x[ix] <- mat[ix]
  return(x)
}
substute.test <- function(){
  x <- matrix(c("Y", "N", "Y", "N", "Y", "N", 1:12), nrow=6)
  e <- as.matrix(c(1,0,0,0,1,0))
  v <- c("R", "S")
###  v <- rep("R", length(x))
  return(substute(x,e,v))
}
cumsumc <- function(mat){
  cumsum(as.data.frame(mat))}


###DESCRIPTION inserts the vector v in the diagonal of mat      
diagrv <- function(mat, v){
  mat[row(mat) == col(mat)] <- 0
  v <- as.vector(v)
  mat <- mat + diag(v)
  return(mat)
}
###DESCRIPTION creates the identity matrix of n rows and n columns
eye <- function(n){ return(diag(n))}
    
###DESCRIPTION obtaines the eigenvalues of a matrix or array of matrices
###
eig <- function(arr){
  dm <- dim(arr)
  if(length(dm) <= 2)
    return(svd(vc,nu=0,nv=0)$d)
  varr <- array(,c(dm[1],dm[2],1))
  for(n in 1:dm[1]){
    mat <- arr[n,,]
    veg <- svd(mat,nu=0,nv=0)$d
    varr[n,,] <- as.matrix(veg)
  }
 return(varr)   
}
###DESCRIPTION inverting matrix mat
###            If the matrix is not invertable then
###            it uses singular value descomposition
###            and calculates the pseudo-inverse.
###            It has the same name as the corresponding
###            Gauss procedure but the algoritms may be different
         
inv <- function(mat,svdtol=1e-10)
{
  S <- try(solve(mat),silent=T)
  if(inherits(S, "try-error"))
    S <- try(solve(mat, LINPACK=TRUE),silent=T)
        
  if(!inherits(S, "try-error"))
    return(S)
  message(paste("Cannot invert matrix ",
                "\n", "Using singular value descomposition", sep=""))
    
  S <- svd.inv(mat,svdtol=svdtol)
  
  return(S)
}

###DESCRIPTION inverting a symmetric, positive definite matrix from Choleski decomposition
###
invpd <- function(mat){
  x <- chol(mat,pivot=TRUE)
  return(chol2inv(x))}

###DESCRIPTION pseudo-inverse in Gauss corresponds to
###            svd.inv in the YourCast software
###
pinv <- function(x,svdtol=1e-10){
  return(svd.inv(x,svdtol=svdtol))}

###DESCRIPTION another version the eig but using function eigen
###            instead of svd (see eig in this file)
###
eig0 <- function(mat){
  return(eigen(mat)$values)}
### For hessian use optim, where fn=&f is pointer to a function
hessp <- function(fn,b){
  return(res <- optim(b,fn, hessian=TRUE)$hessian)}

###DESCRIPTION computes eigenvalues and eigenvectors of a symmetric(hermitian) matrix
eighv <- function(x){ eigen(x,symmetric=TRUE,only.values=FALSE)}
  
reshape <- function(x,r,c){
  x <- as.matrix(x)
  x <- vector(x)
  return(matrix(x,nrow=r,ncol=c,byrow=TRUE))
}
