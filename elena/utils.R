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
  
  if(length(mat) <= 1 && length(v) <= 1)
    return(as.vector(mat)/as.vector(v))
  if(length(v) <= 1 )
    return(mat/as.vector(v))
  if(length(mat) <=1)
    return(as.vector(mat)/as.vector(v))
  if(all(dim(as.matrix(v)) == dim(as.matrix(mat))))
    return(as.matrix(mat)/as.matrix(v))
  mat <- as.matrix(mat)
  
  
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
    return(res <- as.matrix(outer(as.vector(mat0),as.vector(v0), FUN="/")))   
  
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

 if(length(mat) <= 1 && length(v) <= 1)
   return(as.vector(mat)*as.vector(v))
 
  if(length(v) <= 1 )
    return(mat*as.vector(v))

  if(length(mat) <=1)
    return(as.vector(mat)*as.vector(v))
 if(all(dim(as.matrix(v)) == dim(as.matrix(mat))))
    return(as.matrix(as.matrix(mat)*as.matrix(v)))
  mat <- as.matrix(mat)
  
  
  mat0 <- mat
  v0 <- v
  trnps <- FALSE
  if(length(v) > length(mat)){
    v <- as.vector(mat0)
    mat <- as.matrix(v0)
  }else
  v <- as.vector(v)
 
  if((rows(mat0) <= 1 || cols(mat0) <= 1) && (rows(v0) <= 1 || cols(v0) <= 1)){
   
    res <- outer(as.vector(v0),as.vector(mat0), FUN="*")
  ###  res <- t(res)
    return(as.matrix(res))
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
    return(res <- as.matrix(mat) + as.matrix(tam))
  if(length(mat) <= 1 || length(tam) <=1){
   if(length(mat) <= 1) mat <- as.vector(mat)
   if(length(tam) <= 1) tam <- as.vector(tam)
    return(res <- mat + tam)
  }
  if((length(mat) != length(tam)) && (ncm==1 || nrm==1) && (nrt==1 || nct==1))
     stop("NON conformable arrays")
    
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
   
    return(res <- as.matrix(outer(as.vector(mat0),as.vector(tam0),FUN="+")))
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
###DESCRIPTION: uses optim with cml.bounds to estimater parameters with maximum liklehodd
 ###INPUT        stval inital value for the paramneters (px1)
##               lower and upper bounds (px2) for each in stval
###              dataset input to function fn
###              envarinment to stored and obtainbed results
###OUTPUT        the list with the resul;ts of running optim
###
cml.optim <- function(par,cml.bounds,dataset, fn,fctr=1.e+11,hess=TRUE,evbase=get("evbase",env=evbase))
      { 
        ff <- fn
        stval <- par
        stval <- as.matrix(stval)
###defaults
        par <- stval
        con <- list(trace = 0, fnscale = 1, parscale = rep.int(1,length(par)),
                    ndeps = rep.int(0.001, length(par)), maxit = 100, 
                    abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1, 
                    beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5, 
                    factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
###changes
        con$trace <- 1
        con$fnscale <- 1 ##maximizes because function fn returns -res 
        con$REPORT <- 1
        con$factr <- fctr
   
        message("Optim: Covergence acuracy is ", .Machine$double.eps*con$factr);
        message("con$parscale= ", con$parscale)
        message("con$fnscale= ", con$fnscale)

###faster convergence increase con$factr, i.e.  <- 1e+08
###tolerance is defined as .Machine$double.eps*con$factr
        ix <- match(c("reltol", "abstol"), names(con))
        con <- con[-ix]
        message("Optim in action....")
        optimlst <- optim(stval,ff,gr=NULL,dataset,evbase,method="L-BFGS-B",
                          control=con, lower=cml.bounds[,1],upper=cml.bounds[,2],
                          hessian=hess)
        hbool <- hess
        optimnm <- names(optimlst)
      
        ret <- optimlst$convergence
        hess <- if("hessian" %in% optimnm) optimlst$hessian 
        
        lst <- list()
        lst$par <- optimlst$par
        lst$objective <- optimlst$value
        lst$evaluations[[1]] <-  optimlst$counts[[2]]
        lst$evaluations[[2]] <-  optimlst$counts[[1]]
        lst$convergence <- optimlst$convergence
        lst$message <- optimlst$message
        lst$iterations <- NA
        lst$gradient <- NA
        lst$hessian <- hess
        ret <-  lst$convergence
        mess <- paste(lst$message," ;with control$factr= ", fctr) 
    
     ### optim did not converge use nlm 
     if(ret >=1 ){
     
      message("Error with optim...Change control$factr")
      message(mess)
     
      message("Running nml...")
    
      lstnlm <- nlm(fn,par,gradtol=1.e-4, hessian=hess, dataset,evbase)
      lst <- list()
      lst$par <- lstnlm$estimate
      lst$objective <- lstnlm$minimum
      lst$gradient <- lstnlm$gradient
      lst$iterations <- lstnlm$iterations
      lst$hessian <- if(hbool) lstnlm$hessian
      ret <- lstnlm$code
      lst$convergence <- ret
      if(ret >= 3){
        message("code from nlm...", ret)
        stop("Problem occur with nlm and optim change default tolerance")
      }
    }  
        return(lst)
      }

### DESCRIPTION: wrapper around  optim to calculate the hessian only
###              It assumes that the parameters par are already optimized
###              with either optim or nmlinb. Invokes .Internal optimhess
###              Also prints out the time used for the computation
### 
optimhess <- function(par,fn,gr,...,control,Etol=1.e-4,nm=NULL){
 
  initime <- proc.time()
 
  con <- control
  fn1 <- function(par) fn(par,...)
  gr1 <- NULL
  if(!is.null(gr))
     gr1 <- function(par) gr(par,...)
  cml.bounds <- find.bounds(par,Etol)
   hess <- optim(par,fn,gr,...,method="L-BFGS-B",
                 control=con, lower=cml.bounds[,1],upper=cml.bounds[,2],hessian=TRUE)$hessian
###  hess <- .Internal(optimhess(par, fn1, gr1, control))
###  hess <- 0.5 * (hess + t(hess))
  if (!is.null(nm)) 
    dimnames(hess) <- list(nm, nm)
  dtime <- proc.time() - initime
  wd <- sapply(dtime,round)
  wd <- sapply(as.character(wd),nchar)
  dtime <- sapply(1:length(dtime), function(n){
    return(formatC(dtime[[n]], digits=wd[[n]]+3))})
  message("Time consume ...", dtime)
  return(hess)
}
    
### DESCRIPTION : wrapper around nlm to calculate the Hessian
###               It assumes that the parameters p have already being
###               optimized with either optim or nmlinb
###               Returns an object of class nml with the gradient and the
###               the hessian; it also prints out the time used for the computtation
###               It hess=FALSE, it will return the gradient and will also 
###               optimize the parameters p.

nlmhess <- function(fn,p,hess=TRUE,...){
  initime <- proc.time()
 
  res <- nlm(fn,p,hessian=hess,gradtol=1.e-3,steptol=1.e-3,...)
  dtime <- proc.time() - initime
  message("Time consume ...", dtime)
  return(res)
}
relError <- function(x, y){
  x <- as.matrix(x)
  y <- as.matrix(y)
  x[x==0] <- NA
  y[y==0] <- NA
  if(any(dim(x) != dim(y)))
     stop("Inputs need to have same dimensions")
  delx <- abs(x - y)/abs(x+.Machine$double.eps)
  dely <- abs(x-y)/abs(y+.Machine$double.eps)
  errx <- log10(delx+.Machine$double.eps)
  erry <- log10(dely+.Machine$double.eps)
  lst <- c(list(errx=errx), list(erry=erry), list(delx=delx), list(dely=dely))
  return(lst)
}
     
###DESCRIPTION : It calculates the lower and upper bounds to be used with
###              optim and method="L-BGFS-B" after obtaining 
###              estimates for param with nlminb
###              we need the hessian using those values of param
###              and the bounds are narrow limits around param
###
### INPUT: param vector of estimates from nmlinb: px1
###        Edirtol scalar global smaller than 1, i.e 1.e-4
###        perctg percentage of param (cml.bounds) to estimate the bounds
###
### OUPUT : estimated cml.bounds (p x2) around param, which comes from nlminb
###
find.bounds <- function(param,Edirtol=1.e-4,perctg=0.2){
###  delta.bounds <- abs(cml.bounds[,2] - cml.bounds[,1])
###  delta.bounds <- delta.bounds*perctg[1]
###  delta.bounds[delta.bounds< Edirtol] <- Edirtol
###  delta.bounds <- cbind(-abs(param)-delta.bounds,abs(param)+delta.bounds)
  ###bounds

  Edirtol <-  as.vector(Edirtol)
  bounds <- matrix(0, nrow=length(param), ncol=2)
  param <- as.vector(param)
  for(n in 1:length(param)){
  bounds[n,1] <- param[n] - abs(param[n])*perctg[1]
  
 
  bounds[n,2] <- param[n] + abs(param[n])*perctg[1]
  
}

  bounds[bounds[,1] > -Edirtol & bounds[,1] <= 0, 1] <- -Edirtol
  bounds[bounds[,2] < Edirtol & bounds[,2] >= 0, 2] <- Edirtol
  
  return(bounds)
}
