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

 rndtni <- function(m,v,bnds,evbase=parent.frame()){
###  local r,t,sigma,i,lb,ub,inds;
   if(!length(evbase))
     evbase <- get("envbase", env=parent.frame()) 
   lb <- as.matrix(bnds[,1]);
   ub <- as.matrix(bnds[,2]);
   if (nrow(lb)==1)
     lb <- matrix(lb, nrow=rows(m),ncol=1);
 
   if(nrow(ub)==1)
     ub <- matrix(ub, nrow=rows(m),ncol=1);
  
   if (any(v<0)){
     rndtni.v <- v;
     assign("rndtni.v", rndtni.v, env=evbase)
  ###  @v=recode(v,v.<0,1e-10);@
     stop("rndtni: negative variance; see rndtni.v")
   }
   t <- lb>ub;
   if(any(t))
     stop("rndtni: upper bound less than lower bound!")
 
   sigma <- sqrt(v);
   fcmptol <- 1e-12;
   t <- 1-dotfeq(lb,ub,tol=fcmptol); ###1 -(lb==ub)
   sigma <- t*sigma;
   m <- t *m +(1-t)*lb;

   r <- m+matrix(rnorm(rows(m), mean=0, sd=1), nrow=rows(m), ncol=1)* sigma;
   t <- (r<lb)| (r>ub);
   i <- 1;
   while( i<5 | colSums(t)!=0){
###   /* sample rejection method */
     inds <- grep(1, t)
     if(length(inds)){
       r[inds] <- m[inds]+ matrix(rnorm(rows(inds), mean=0, sd=1), nrow=rows(inds), ncol=1, byrow=T)*sigma[inds];
       t <- (r<lb)|(r>ub);
     }
     i <- i+1;
   }
   if(colSums(t)!=0){
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
  
  bwsims <- rndtni(bw*o,sw2*o,bounds[2,]*o);

  m <- bb+(sbw/sw2)*(bwsims-bw);
  v <- sb2-((sbw^2)/sw2);
  bbsims <- rndtni(m,v*o,bounds[1,]*o);
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
  
  res <- qnorm(p*(crgt-clft)+clft,mean=mu,sd=sigma);
  tL <- (res<lft);
  tR <- (res>rgt);
  ok <- (res>=lft) & (res<=rgt);
  res <- res*ok+lft*tL+rgt*tR;
  tL <- (tL+tR);
  t <- colSums(tL);
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
