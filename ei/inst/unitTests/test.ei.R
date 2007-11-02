.setUp <-function(){
require(ei)
}

###
## runit function for testing bounds1()
##
test.bounds1 <- function()
{
  tol <- 1.e-4
  data(sample)
  bnd <- bounds1(t,x,1,1.e-4)      
  load("bounds.RData")
  bs <- y.bounds[[1]]
  agg <- y.bounds[[2]]
  ## testing ...
       
  checkEquals(bs,bnd[[1]],tolerance=tol)
  checkEquals(agg,bnd[[2]],tolerance=tol)
      
}

##
test.nonbiv <- function()
{
  tol <- 1.e-4
  data(sample)
           
  load("betab.RData")
  load("betaw.RData")
  pz <- nonbiv(t,x,betab,betaw, evbase=NULL, eigraph.bvsmth=0.08,Enumtol=1.e-4)
  load("nonbiv.RData")
  checkEquals(y.nonbiv,pz,tolerance=tol)
        
           
}
##
test.betas <- function()
{
      
  tol <- 1.e-4
  data(sample)
           
  load("betab.RData")
  load("betaw.RData")
  load("bounds.RData")
  nobs <- rows(t);
  EnonEval <- 11
  bnds <- bounds[[1]]
  tt <- bounds[[2]]

  betab.test <- matrix(0, nrow=EnonEval,ncol=nobs);
  for (i in 1:nobs)
    betab.test[,i+0] <- seq(from=bnds[i+0,1],to=bnds[i+0,2],length.out=EnonEval)
      
  x1 <- 1-x;
 
  betaw1 <- betab.test
  for(n in 1:ncol(betab.test))
    betaw1[,n] <- (x/x1)[n]* betab.test[,n]
  for(n in 1:ncol(betaw1))
    betaw1[,n] <- (t/x1)[n] - betaw1[,n]
  
  checkEquals(betab.test,betab,tolerance=tol)
  checkEquals(betaw1,betaw,tolerance=tol)
           
}

test.unitarea <- function()
{
  tol <- 1.e-4
  data(sample)
  uar <- unitarea(t,x, evbase=NULL,EnonNumInt=11, Enumtol=1.e-4,eigraph.bvsmth=0.08)
  load("unitarea.RData")
  checkEquals(uar,y.unitarea,tolerance=tol)
 
}

  
