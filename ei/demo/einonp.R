
einonp <- function()
{
 
  message("Running non-parametric estimation")
 
###  verb <- user.prompt()
  message("Loading the data sample")
  res <- data(sample)
  res <- get(res, env=environment())
  t <- res[[1]]
  x <- res[[2]]
  n <- res[[3]]
###    verb <- user.prompt()
  message("Running non-parametric estimation")
  dbuf <- ei(t,x,n,1,1,EnonPar=1)
  print(names(dbuf))
  
  
}

einonp()
