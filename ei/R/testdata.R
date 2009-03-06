##files are: sample.asc, sample2.asc
readat <- function(datapath="~/GAUSSCODE/ei/data", file="pa90.txt",skp=1,ncl=5){
  
  datstring  <- paste(datapath,"/",file,sep="")
  message(paste("Reading file ", datstring, sep=""))
          dat0 <- scan(file=datstring,na.strings="NA",skip=skp, 
                       multi.line=T,quiet=T)
          dat0 <- matrix(dat0,ncol=ncl,byrow=T)
return(dat0)
}
getASCII <- function(eps =.Machine$double.eps,file="pa90.txt",nocol=5 ){

mat <- readat(file=file,ncl=nocol)
if(nocol==5){
  t0 <- mat[, 3]
  x0 <- mat[,1]
  invtvap <- mat[,5]
  tvap <- 1/(invtvap +eps )
}else if(nocol==2){
  t0 <- mat[, 2]
  x0 <- mat[,1]
  tvap <- rep(1, length(t0))
}
xind <- which(x0 <= 0 | x0 >= 1)
tind <- which(t0 <= 0 | t0 >= 1)
nind <- which(tvap<=0)
ind <- unique.default(c(xind,tind,nind))
if(length(ind)) {
  x0 <- x0[-ind]
  t0 <- t0[-ind]
  tvap <- tvap[-ind]
}
assign("t0",t0, env=.GlobalEnv)
assign("x0",x0, env=.GlobalEnv)
assign("tvap",round(tvap), env=.GlobalEnv)
message("Names of variables are t0 for t, x0 for x, tvap for tvap and n")
}
