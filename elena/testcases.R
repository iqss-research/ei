##files are: sample.asc, sample2.asc
readvars <- function(datapath="~/cvswork/ei/ei/data", file="sample.asc"){
  
  datstring  <- paste(datapath,"/",file,sep="")
  message(paste("Reading file ", datstring, sep=""))
          dat0 <- scan(file=datstring,
                       na.strings="NA",
                       multi.line=T,quiet=T)
          dat0 <- matrix(dat0,ncol=3,byrow=T)
return(dat0)
}
mat <- readvars()
t <- mat[, 1]
x <- mat[,2]
tvap <- mat[,3]
