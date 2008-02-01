 
  message("Running non-parametric estimation")
 
###  verb <- user.prompt()
  message("Loading the data sample")
  res <- data(sample)
 
  t <- sample[[1]]
  x <- sample[[2]]
  n <- sample[[3]]
###    verb <- user.prompt()
  message("Running non-parametric estimation")
  dbuf <- ei(t,x,n,1,1,EnonPar=1)
  print(names(dbuf))
  message("Obtaining beta blacks")  
  betab <- dbuf$betabs
  message("Calculating beta whites")
  betaw <- betab2w(t,x,betab)



