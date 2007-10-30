## Some utilities that Gary King wrote in Gauss for library ei translated
## into R code and some other utilities that Elena wrote in R that has
## has some corresponding functions in Gauss
## Defined in original Gauss code are:
## scalzero,scalone, meanwc (meanWc), ismissm, seqas, seqase, vin,
## sortbyRow, mkmissm, fisherzi,fisherz, ftosm, fmtt, infrv,  

## Added to provide same functionality as in Gauss
## %inG% (same as in Gauss), extract.diag (same as diag in Gauss)
## trim.blanks, %dot/%, %dot*%, %-%
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

  return(weighted.mean(x,wt,na.rm=na.remove))
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

### DESCRIPTION the gauss code is case insensitive and that it was
###             inG is doing, which is a wrapper around %in% but ignoring case
###             Note that "in" may return a vector of T and F but inG applies
###             all to the results of in and it is only T if all vector elemnts
###             are T.
###
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
###
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

 
###DESCRIPTION : trims leading and trailing blanks for any string but not
###              inter-words blanks.  For example, trim.blanks("   abc  ") = "abc";
###              but trim.blanks("   abc    def ")= "abc    def".
### AUTHOR: Elena Villalon
##          evillalon@iq.harvard.edu
###
trim.blanks <- function(x) {
### at the beginning of string"^" gets anny number (+) of white spaces
  f <- x
  if(length(x))
    f <- na.omit(x)
  
  if(length(f) <= 0)
    return(x)
  if(length(f)>1)
    print(f)
  if(f=="" )
    return(x)
  x <- sub("^[[:space:]]*(.*)", "\\1",x) ###get \n\t
  x <- sub('^ +', '', x) ###get white spaces
  
### at the ending of string"$" gets anny number (+) of white spaces
  
  x <- sub("(.*)[[:space:]]*$", "\\1", x)
  x <- sub(' +$', '', x)
  return(x)
}



### DESCRIPTION this is the element division equivalent
###             to the Gauss and Matlab operation "./"
### INPUT v=vector and mat=matrix
###       v may have only one element or length(v) = ncol(mat) or
###       length(v) = nrow(mat)
###
### OUTPUT divide elemnt by element mat/v or t(t(mat)/v)
###        depending if nrow(mat) or ncol(mat) = length(v)
###
###        Gauss mat./v -or- v./mat 
###
###AUTHOR Elena Villalon
###       evillalon@iq.harvard.edu
###
"%dot/%" <- function(mat, v){

  if(length(v) <= 1 || length(mat) <= 1)
    return(mat/v)
  
  mat <- as.matrix(mat)
  
  if(length(v)==length(mat)){
    res <- as.vector(mat) /as.vector(v)
    res <- matrix(res,nrow=rows(mat), ncol=cols(mat))
    return(res)
  }
  mat0 <- mat
  v0 <- v
  invrt <- FALSE
  trnps <- FALSE
  if(length(v) > length(mat)){
    v <- as.vector(mat0)
    mat <- as.matrix(v0)
    invrt <- T
  }else
  v <- as.vector(v)
    
  if(length(v) != ncol(mat) && length(v) != nrow(mat))
    stop("Arrays are non-conformable")
  if(length(v) == rows(mat)) {
    mat <- t(mat)
    trnps <- T
  }
  if(!is.data.frame(mat)){
    mdf <- as.data.frame(mat)
    res <- as.matrix(t(t(mdf)/v))
    colnames(res) <- NULL
    rownames(res) <- NULL
    if(invrt)
      res <- 1/res
    if(trnps)
      res <- t(res)
    return(res)
  }else{
    if(!invrt)
      res <- t(mat/v)
    else
      res <- 1/t(mat/v)
  if(trnps)
    res <- t(res)
  return(res)
  }
}
### DESCRIPTION this is the element multiplication equivalent
###             to the Gauss and Matlab operation ".*"
### INPUT v=vector and mat=matrix
###       v may have only one element or length(v) = ncol(mat) or
###       length(v) = nrow(mat)
###
### OUTPUT multiply elemnt by element mat*v or t(t(mat)*v)
###        depending if nrow(mat) or ncol(mat) = length(v)
###
###        Gauss v.*mat -or- mat*.v 
###
###AUTHOR Elena Villalon
###       evillalon@iq.harvard.edu
###
"%dot*%" <- function(v, mat){

  if(length(v) <= 1 || length(mat) <=1)
    return(mat*v) 
  mat <- as.matrix(mat)
  if(length(v)==length(mat)){
    res <- as.vector(v) *as.vector(mat)
    res <- matrix(res,nrow=rows(mat), ncol=cols(mat))
    return(res)
  }
mat0 <- mat
v0 <- v

trnps <- FALSE
if(length(v) > length(mat)){
  v <- as.vector(mat0)
  mat <- as.matrix(v0)
}else
v <- as.vector(v)
     

   if(length(v) != ncol(mat) && length(v) != nrow(mat))
    stop("Arrays are non-conformable")

  if(length(v) == rows(mat)) {
    mat <- t(mat)
    trnps <- T
  }

  if(!is.data.frame(mat)){
    mdf <- as.data.frame(mat)
    res <- as.matrix(t(t(mdf)*v))
    colnames(res) <- NULL
    rownames(res) <- NULL
    if(trnps) res <- t(res)
    return(res)
  }else{
    res <- (t(t(mat)*v))
    if(trnps) res <- t(res)
    return(res)
  }
}
###DESCRIPTION difference two matrices
### mat = Mx 1
### v= Vx 1
### res = M x V
### Gauss (mat-m'), where mat and m are 1 column matrices and m'=t(m)
## 

 "%-%" <- function(mat, m){
     res <- matrix(as.vector(mat),nrow=length(mat), ncol=length(m))
     red <- t(as.data.frame(t(res))-m)
     res <-  as.matrix(red)
     colnames(res) <- rownames(res) <- NULL
     return(res)
   }
     
     
         
    
    

  
