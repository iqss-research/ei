##files are: sample.asc, sample2.asc
readvars <- function(datapath="~/cvswork/ei/ei/data", file="sample.txt"){
  
  datstring  <- paste(datapath,"/",file,sep="")
  message(paste("Reading file ", datstring, sep=""))
          dat0 <- scan(file=datstring,
                       na.strings="NA",
                       multi.line=T,quiet=T)
          dat0 <- matrix(dat0,ncol=3,byrow=T)
return(dat0)
}
getSample <- function(){

mat <- readvars()
t0 <<- mat[, 1]
x0 <<- mat[,2]
tvap <<- mat[,3]
message("Names of variables are t0 for t, x0 for x, tvap for tvap and n")
}
