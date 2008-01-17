## random numbers from truncated normal distribution
## via sample rejection method with hard cases done via CDF method
##
## USAGE:  r = rndtni(m,v,bounds);
##
## INPUT: all inputs have N rows
## m = vector of means
## v = vector of variances
## bounds = upper-bound ~ lower-bound
##
## OUTPUT:
## r = nx1 vector of independent random numbers with means m and variances v
##

 rndtni <- function(m,v,bnds,evbase=get("evbase", env=parent.frame())){
###  local r,t,sigma,i,lb,ub,inds;
   if(!length(evbase))
     evbase <- get("envbase", env=parent.frame())

   lb <- as.matrix(bnds[,1]);
   ub <- as.matrix(bnds[,2]);
   if (nrow(lb)==1)
     lb <- as.vector(lb) * matrix(1, nrow=rows(m),ncol=1);
 
   if(nrow(ub)==1)
     ub <- as.vector(ub) * matrix(1, nrow=rows(m),ncol=1);
  
   if (any(v<0)){
     rndtni.v <- v;
     assign("rndtni.v", rndtni.v, env=evbase)
  ###  @v=recode(v,v.<0,1e-10);@
     stop("rndtni: negative variance; see rndtni.v")
   }
   t <- (lb>ub);
   if(any(t))
     stop("rndtni: upper bound less than lower bound!")
 
   sigma <- sqrt(v);

   fcmptol <- 1e-12;
   t <- 1-dotfeq(lb,ub,tol=fcmptol); ###1 -(lb==ub)
  
   sigma <- as.matrix(t) %dot*% as.matrix(sigma);
  
   m <- as.matrix(t) %dot*% as.matrix(m) +as.matrix(1-t) %dot*% lb;
   
   r <- m+matrix(rnorm(rows(m), mean=0, sd=1), nrow=rows(m), ncol=1) %dot*% as.matrix(sigma);
   t <- (r<lb)| (r>ub);
  
   for(i in 1:5){
###   /* sample rejection method */
    if(colSums(t) ==0) break; 
     inds <- grep(1, t)
     if(length(inds)){
       r[inds] <- m[inds]+ matrix(rnorm(rows(inds), mean=0, sd=1), nrow=rows(inds), ncol=1, byrow=T) %dot*% as.matrix(sigma[inds]);
       t <- (r<lb)|(r>ub);
     }
   
    
   }
   if(colSums(as.matrix(t))!=0){
###  /* sample rejection fails for some elements; try CDF method */
     inds <- grep(1, t)
     if(length(inds))
     
      r[inds]=invcdftn( as.matrix(runif(length(inds))),m[inds],v[inds],lb[inds],ub[inds]);
   }
  return(r);
 }

##    a =  rndbtn(bb,bw,sb,sw,rho,bounds,sims);
##
## bivariate truncated normal random numbers
##
## inputs: bb = 1st mean
##         bw = 2nd mean
##         sb = 1st standard deviation
##         sw = 2nd standard deviation
##        rho = correlation
##     bounds = 2x2 lower~upper for 1st mean in 1st row and 2nd in 2nd
##       sims = number of simulations
##
## output:  a = sims x 2 matrix of  BIvariate Truncated Normal Random Variables
##              each row of a is one 1x2 simulation
##

rndbtn <- function(bb,bw,sb,sw,rho,bounds,sims, evbase=parent.frame()){
  if(!length(evbase))
    evbase <- get("envbase", env=parent.frame()) 
  o <- matrix(1, nrow=sims,ncol=1);
  sb2 <- sb^2;
  sw2 <- sw^2;
  sbw <- rho*sb*sw;
  mat <- as.vector(o)%dot*%bounds[2,]  ##outer product in Gauss language

  bwsims <- rndtni(bw*o,sw2*o,mat)

  m <- bb+(sbw%dot/%sw2)%dot*%(bwsims-bw);
  v <- sb2-((sbw^2)%dot/%sw2);
   mat <- as.vector(o)%dot*%bounds[1,] 
  bbsims <- rndtni(m,v%dot*%o,mat);
  mat <- cbind(bbsims,bwsims)
  return(mat);
}


## inverse of the truncated normal CDF
##
## USAGE:    y = invcdftn(p,mu,sigma2,left,right);
##
## INPUTS: all Nx1
## mu = mean
## sigma2 = variance
## p = Prob(Y<y|mu,sigma2), where y is a realization
##                          of the random variable Y
## left = left truncation bound
## right = right truncation bound, left <= Y <= right
##
## OUTPUT: Nx1
## y = normal variate

invcdftn <- function(p,mu,sigma2,lft,rgt){
###  local t,res,ok,tL,tR,clft,crgt;
  if( lft>=rgt)
    stop("invcdftn: left must be < right")


  if(any(p>=1) ||  any(p<=0))
    stop("invcdftn: input must be (0,1)")
  if(sigma2 <=0){
    message("variance must be positive")
    return(NULL)
  }
  sigma <- sqrt(sigma2)
  clft <- pnorm(lft,mean=mu,sd=sigma);
  crgt <- pnorm(rgt,mean=mu,sd=sigma);
  
  res <- qnorm(p%dot*%(crgt-clft)+clft,mean=mu,sd=sigma);
  tL <- (res<lft);
  tR <- (res>rgt);
  ok <- (res>=lft) & (res<=rgt);
  res <- res%dot*%ok+lft%dot*%tL+rgt%dot*%tR;
  tL <- (tL+tR);
  t <- colSums(as.matrix(tL));
  if(t!=0){
 ###   /*t=seqa(1,1,rows(p));
 ###   t=selif(t,tL);*/
    nc  <- nchar(as.character(floor(t))) 
    fmt <- formatC(t,width=nc,digits=0, format="f")
    message("invcdftn: Warning: Some bounds are very far from distribution mean.", 
    "          Forcing ", fmt, " simulations to their closest bound")
  }
  return(res);
}
##/*
##** Random draws from a truncated singular multivariate normal density
##**
##**  res = rndtsn(mean,negHess,sims,bounds,tol);
##**
##** INPUTS:
##**   mean    = k x 1 vector of means
##**   negHess = k x k matrix, the negative of the hessian
##**   sims    = scalar, the number of simulations
##**   bounds  = k x 2 matrix, upper bounds ~ lower bounds
##**   tol     = scalar, the tolerance level for eigenvalues (e.g., 1e-12)
##**
##** OUTPUTS:
##**   res = sims x k matrix of truncated multivariate singular normal random draws
##**
##*/
rndtsn <- function(mean,invvc,sims,bounds,tol,Eprt=2){
   
   k <- rows(mean);
  ### /* input checks */
   
   if (any(bounds[,2] <bounds[,1])){
     message("error(rndmnsvd): Upper bounds must be greater than lower bounds.")
     return(NA)
   }
   if( rows(invvc)!=cols(invvc)){
     message("error(rndmnsvd): the -Hessian must be square.")
     return(NA)
   }
   if( k!=rows(invvc)){
     message("error(rndmnsvd): the dimensions of the -Hessian matrix and ", 
             "mean vector must be the same.")
     return(NA)
   }
    lst <- eigen(invvc,symmetric=TRUE,only.values=FALSE)
   s <- lst[[1]]
   u <- lst[[2]]
   ### @ spectral decomposition: u*diag(s)*u'=invvc @
   s[s<tol] <- tol
 
   snum <- sims*10 ###          @ # of draws at a time @
 
   ###/* shift to the middle of the bounds, create (k x snum) matrices */
   midbounds <- colSums(t(as.matrix(bounds)))/2
   tmpmean <- t(u)*(mean-midbounds) ###        @ transformed mean @
   tmpmean <- as.matrix(as.vector((tmpmean %dot*% matrix(1,nrow=k,ncol=snum))))
   v <- as.matrix(as.vector((1./s)%dot*% matrix(1,nrow=k,ncol=snum))) ###        @ transformed variance @
   temp <- bounds-midbounds
   temp <- maxc(t(abs(temp)));
   dist <- sqrt(t(temp)%*%temp)
   bnds <- cbind(-dist, dist) ###                 @ transformed bounds @

   limsim <- 1              
   res <- matrix(0,nrow=1,ncol=k) 
   bnds1 <- bnds%dot*% matrix(1,nrow=rows(v),ncol=2)
   while (rows(res)<=sims){
     temp <- rndtni(tmpmean,v,bnds1)
     temp <- t(matrix(as.vector(temp),nrow=rows(temp)/k,ncol=k,byrow=TRUE))
     temp <- u%*%temp+midbounds;
     t <- colSums(as.matrix((bounds[,1] >as.vector(temp))| (bounds[,2] <as.vector(temp))));
     temp <- subset(t(temp), subset=(t==0))
     if(!length(temp)) temp <- NA

     if( !scalmiss(temp))
       res <- rbind(res,temp)
    
     if( as.vector(Eprt) >=2)
       message("rndtsn: trying ", limsim, "th time...")
     
     limsim <- limsim+1
     if (limsim==100){
       message("error(rndtsn): the sampling method failed. adjust the bounds.")
       return(NA)
     }
   }
   return(res[2:sims+1,])
 }
