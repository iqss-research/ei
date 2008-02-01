
 
  message("Running parametric estimation:Ecdfbvn=1")
 
###  verb <- user.prompt()
  message("Loading the data sample")
  res <- data(sample)
 
  t <- sample[[1]]
  x <- sample[[2]]
  n <- sample[[3]]
###    verb <- user.prompt()
  message("Running non-parametric estimation")
  dbuf <- ei(t,x,n,1,1,Ecdfbvn=1)
  print(names(dbuf))
  
  message("Obtaining beta blacks")  
  betab <- dbuf$betaBs
  message("Calculating beta whites")
  betaw <- betab2w(t,x,betab)
