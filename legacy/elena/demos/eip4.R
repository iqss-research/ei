
 
  message("Running parametric estimation:Ecdfbvn=4")
 
###  verb <- user.prompt()
  message("Loading the data sample")
  res <- data(sample)
  res <- get(res, env=environment())
  t <- res[[1]]
  x <- res[[2]]
  n <- res[[3]]
###    verb <- user.prompt()
  message("Running non-parametric estimation")
  dbuf <- ei(t,x,n,1,1,Ecdfbvn=4)
  print(names(dbuf))
  
  message("Obtaining beta blacks")  
  betab <- dbuf$betaBs
  message("Calculating beta whites")
  betaw <- betab2w(t,x,betab)
