### extract.diag, %inG%,%dot/%, %dot*%,%-%, %plus%, %minus%
### svd.inv, trim.blanks, user.prompt, messout
###

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
### DESCRIPTION this is the element division equivalent to the 
###             Gauss and Matlab operation "./" (the invisible dot)
###             
### INPUT v=vector and mat=matrix
###       v may have only one element or length(v) = ncol(mat) or
###       length(v) = nrow(mat)
###       If mat = Mx1 and v =Kx1 (outer product)
###
### OUTPUT divide elemnt by element mat/v or t(t(mat)/v)
###        depending if nrow(mat) or ncol(mat) = length(v)
###
###        Gauss mat./v -or- v./mat 
###
###AUTHOR Elena Villalon
###       evillalon@iq.harvard.edu
###
"%./%" <- "%dot/%" <- function(mat, v){
  if(length(v) <= 1 )
    return(mat*as.vector(v))
  if(length(mat) <=1)
    return(as.vector(mat)*as.vector(v))
   
  mat <- as.matrix(mat)
  
  if(all(dim(v)==dim(mat))){
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
  
  if((rows(mat0) <= 1 || cols(mat0) <= 1) && (rows(v0) <= 1 || cols(v0) <= 1))
    return(res <- outer(as.vector(mat0),as.vector(v0), FUN="/"))   
 
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
     stop("Arrays are non-conformable")
}
### DESCRIPTION this is the element multiplication equivalent to the 
###             Gauss and Matlab operation ".*" (the invisible dot)
###
### INPUT v=vector and mat=matrix
###       v may have only one element or length(v) = ncol(mat) or
###       length(v) = nrow(mat)
###       If mat = Mx1 and v =Kx1 (outer division)
###
### OUTPUT multiply elemnt by element mat*v or t(t(mat)*v)
###        depending if nrow(mat) or ncol(mat) = length(v)
###
###        Gauss v.*mat -or- mat*.v 
###
###AUTHOR Elena Villalon
###       evillalon@iq.harvard.edu
###
"%.*%" <- "%dot*%" <- function(v, mat){

  if(length(v) <= 1 )
    return(mat*as.vector(v))
  if(length(mat) <=1)
    return(as.vector(mat)*as.vector(v))
  mat <- as.matrix(mat)
  if(all(dim(v)==dim(mat))){
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
     
  if((rows(mat0) <= 1 || cols(mat0) <= 1) && (rows(v0) <= 1 || cols(v0) <= 1)){
   
    return(res <- outer(as.vector(v0),as.vector(mat0), FUN="*")) 
  }

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
   stop("NON conformable arrays")
}
###DESCRIPTION difference two matrices
### mat = Mx 1
### v= Vx 1
### res = M x V
### Gauss (mat-m'), where mat and m are 1 column matrices and m'=t(m)
###
##  This is a praticular case mat %minus% m  

 "%-%" <- function(mat, m){
     res <- matrix(as.vector(mat),nrow=length(mat), ncol=length(m))
     red <- t(as.data.frame(t(res))-m)
     res <-  as.matrix(red)
     colnames(res) <- rownames(res) <- NULL
     return(res)
   }
     
###DESCRIPTION adding two matrices
### mat = M x K, and m = M x K or m = 1 x 1
### output = M x K 
### mat = M x K, and m= M x 1 or m= K x 1
### output = M x K
### mat = M x 1 and m = 1 x K
### output= M x K
###  
### AUTHOR: Elena Villalon (evillalon@iq.harvard.edu)
###    
"%plus%" <- function(mat, tam)
{
  mat0 <- mat
  tam0 <- tam
  mat <- as.matrix(mat)
  tam <- as.matrix(tam)
  nrm <- nrow(mat)
  ncm <- ncol(mat)
  nrt <- nrow(tam)
  nct <- ncol(tam)
  if(nrm==nrt && ncm==nct)
    return(res <- mat + tam)
  
  if(length(mat) <= 1 || length(tam) <=1)
    return(res <- mat + tam)
  if((ncm <=1 || nrm <= 1) && (length(mat)==nrt || length(mat) == nct)){
    vec <- as.vector(mat)
    return(res <- as.data.frame(tam) + vec)
  }
  if((nct <=1 || nrt <= 1) && (length(tam)==nrm || length(tam) == ncm)){
      vec <- as.vector(tam) 
      return(res <- as.data.frame(mat) + vec)
    }
 
  if((nrm <= 1 && (nrt <= 1 || nct <= 1)) ||
     (ncm <= 1 && (nrt <= 1 || nct <= 1))){
   
    return(res <- outer(as.vector(mat0),as.vector(tam0),FUN="+"))
  }
  
   stop("NON conformable arrays")
}
    
### DESCRIPTION Substracting matrices (same as %plus%)

"%minus%" <- function(mat, tam){
 
  return(res <- mat %plus% (-tam))
 }

### expression like (delta > 1) in Gauss means (all(delta > 1)
### If relational operators are not preceded with a dot="."
### the result is either TRUE (1) or FALSE(0) based upon
## a comparison of all elemnts in x and y, e.g. x > y
## If relational op preceded with dot="." then the result is
## either a matrix or a vector of elemnts 
## TRUE(1) and FALSE(0), for each element comparison,e.g. x .> y

"%>%"  <- function(x, y) {return(all(x > y))}  
"%.>%" <- "%dot>%" <- function(x, y) {return (x > y)}

"%>=%"  <- function(x, y) {return(all(x >= y))}
"%.>=%" <- "%dot>=%" <- function(x, y) {return (x >= y)}

"%<%" <- function(x, y) {return(all(x < y))}
"%.<%" <- "%dot<%" <- function(x, y) {return (x < y)}

"%<=%" <- function(x, y) {return(all(x <= y))}
"%.<=%" <- "%dot<=%" <- function(x, y){return(x <= y)}

"%==%" <- function(x, y) {return(all(x == y))}
"%.==%"  <- function(x, y) {return(x == y)}


   
  
###DESCRIPTION pseudo-inverse of matrix mat
###            using the singular value descomposition.
###
### Elena Villalon (evillalon@iq.harvard.edu)

svd.inv <-  function(x,svdtol=1e-10){
  s <- svd(x);
  w <- s$d;
  U <- s$u;
  V <- s$v;
  w.inv <- 0*w;
  w.inv[w > svdtol] <- 1/w[w > svdtol];
  if(length(w.inv) <= 1)
    w.inv <- as.matrix(w.inv)
  W.inv <- diag(w.inv);
  return(V %*% W.inv %*% t(U)); 
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
###Request input from user 
user.prompt <- function (verbose=TRUE){
  verb <- T
 if(length(verbose) <= 0 || verbose){
   silent <- readline("\nPress <return> to continue or Ctrl-c Ctrl-c to quit: ")
   
 }else{
   
   answer <- readline("\nPress 'Y' for verbose output or enter 'N' otherwise. \nPress Ctrl-c Ctrl-c to quit: " )
   
   if(substr(answer,1,1)=='N')
     verb <- F
   else if(substr(answer,1,1)=='Y')
     verb <- T
   
 }
   return(verb)
     
}
### Used instead of message to control verbose output
messout <- function(str, verbose=T, obj=NULL){
  if(verbose)
  message(str);
  if(length(obj) >0 && verbose)
    print(obj)
}
 

