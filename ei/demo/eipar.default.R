einonp <- function()
{
 
  message("Running parametric estimation:Ecdfbvn=5")
 
###  verb <- user.prompt()
  message("Loading the data sample")
  res <- data(sample)
  t <- res[[1]]
  x <- res[[2]]
  n <- res[[3]]
###    verb <- user.prompt()
  message("Running non-parametric estimation")
  dbuf <- ei(t,x,n,1,1,EnonPar=1)
  print(names(dbuf))
  
  
}

einonp()
