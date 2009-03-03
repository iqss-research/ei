##files are: sample.asc, sample2.asc
readat <- function(datapath="~/GAUSSCODE/ei/data", file="pa90.asc",skp=1){
  
  datstring  <- paste(datapath,"/",file,sep="")
  message(paste("Reading file ", datstring, sep=""))
          dat0 <- scan(file=datstring,na.strings="NA",skip=skp, 
                       multi.line=T,quiet=T)
          dat0 <- matrix(dat0,ncol=5,byrow=T)
return(dat0)
}
getASCII <- function(eps =.Machine$double.eps ){

mat <- readat()
t0 <- mat[, 3]
assign("t0",t0, env=.GlobalEnv)
x0 <- mat[,1]
assign("x0",x0, env=.GlobalEnv)

invtvap <- mat[,5]
tvap <- 1/(invtvap +eps )
assign("tvap",tvap, env=.GlobalEnv)
message("Names of variables are t0 for t, x0 for x, tvap for tvap and n")
}
