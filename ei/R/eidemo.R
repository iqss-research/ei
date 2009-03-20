###DESCRIPTION to run demos for non-parametric estimations
###
###INPUT: dat the data as data.frame or mat
###       tind,xind,nind numeric indeces to identify x,t,n in dat
###       invn, boolean to invert dat[[nidn]]
### Run the demo
###
eidemononpar <- function(dat,tind=1,xind=2,nind=3,invn=FALSE,...){
res <- dat
t <- res[[tind]]
x <- res[[xind]]
invtvap <- res[[nind]]
if(invn)
  tvap <- 1/(invtvap +.Machine$double.eps)
else
  tvap <- invtvap
ind <- which(x <= 0 | x >= 1 |t <= 0 | t >= 1 |tvap<=0)
if(length(ind)){
  ind <- unique.default(ind)
  x <- x[-ind]
  t <- t[-ind]
  tvap <- tvap[-ind]
}
n <- round(tvap)
message("Running non-parametric estimation")
###user.prompt()
dbuf <- ei(t,x,n,1,1,EnonPar=1,dbug=FALSE,...)
print(names(dbuf))
message("Running graphics:") 
eigraph(dbuf,"tomog")
user.prompt()
eigraph(dbuf,"tomogp")
user.prompt()
eigraph(dbuf,"tomoge")
user.prompt()
eigraph(dbuf,"tomogci")
user.prompt()
eigraph(dbuf,"tomogci95")
user.prompt()
eigraph(dbuf,"estsims")
user.prompt()
eigraph(dbuf,"tomogs")
user.prompt()
eigraph(dbuf,"nonpar")
user.prompt()
eigraph(dbuf,"xtc")
user.prompt()
eigraph(dbuf,"xt")
user.prompt()
eigraph(dbuf,"xgraph")
user.prompt()
eigraph(dbuf,"xgraphc")
user.prompt()
eigraph(dbuf,"goodman")
user.prompt()
eigraph(dbuf,"profile")
user.prompt()
eigraph(dbuf,"profileR")
user.prompt()
eigraph(dbuf,"postb")
user.prompt()
eigraph(dbuf,"postW")
user.prompt()
eigraph(dbuf,"post")
user.prompt()
message("Running beta with kern=E")
eigraph(dbuf,"betaB")
user.prompt()
eigraph(dbuf,"betaW")
user.prompt()
eigraph(dbuf,"beta")
user.prompt()
message("Running beta with kern=TN")
eigraph(dbuf,"beta",kern="TN")
user.prompt()
eigraph(dbuf,"results", kern="E")
user.prompt()
eigraph(dbuf,"movie")
user.prompt()
eigraph(dbuf,"movied")
user.prompt()
eigraph(dbuf,"lines")
user.prompt()
eigraph(dbuf,"bivar")
user.prompt()
eigraph(dbuf,"betabw")
user.prompt()
eigraph(dbuf,"biasb")
user.prompt()
eigraph(dbuf,"biasw")
user.prompt()
eigraph(dbuf,"bias")
user.prompt()
eigraph(dbuf,"boundxb")
user.prompt()
eigraph(dbuf,"boundxw")
user.prompt()
eigraph(dbuf,"boundx")
user.prompt()
message("Addition: dependences of beta's vs T,N")
eigraph(dbuf,"betast")
user.prompt()
eigraph(dbuf,"betasn")
user.prompt()
message("Addition: three-dimensional dependences of beta's vs X,T,N")
eigraph(dbuf,"betaxn")
user.prompt()
eigraph(dbuf,"betatn")
user.prompt()
}
###DESCRIPTION to run demos for parametric estimations
###
###INPUT: dat the data as data.frame or mat
###       tind,xind,nind numeric indeces to identify x,t,n in dat
###       invn, boolean to invert dat[[nidn]]
### Run the demo
###
eidemopar <- function(dat,tind=3,xind=1,nind=5,invn=FALSE,...){
res <- dat
t <- res[[3]]
x <- res[[1]]
invtvap <- res[[5]]
if(invn)
  tvap <- 1/(invtvap +.Machine$double.eps)
else
  tvap <- invtvap
xind <- which(x <= 0 | x >= 1)
tind <- which(t <= 0 | t >= 1)
nind <- which(tvap<=0)
ind <- unique.default(c(xind,tind,nind))
if(length(ind)) {
  x <- x[-ind]
  t <- t[-ind]
  tvap <- tvap[-ind]
}
n <- round(tvap)
message("Running default parametric estimation")
###user.prompt()
dbuf <- ei(t,x,n,1,1,EdoML=1,dbug=FALSE,...)
print(names(dbuf))
message("Running graphics:") 
eigraph(dbuf,"tomog")
user.prompt()
eigraph(dbuf,"tomogp")
user.prompt()
eigraph(dbuf,"tomoge")
user.prompt()
eigraph(dbuf,"tomogci")
user.prompt()
eigraph(dbuf,"tomogci95")
user.prompt()
eigraph(dbuf,"estsims")
user.prompt()
eigraph(dbuf,"tomogs")
user.prompt()
eigraph(dbuf,"nonpar")
user.prompt()
eigraph(dbuf,"xtc")
user.prompt()
eigraph(dbuf,"xt")
user.prompt()
eigraph(dbuf,"xgraph")
user.prompt()
eigraph(dbuf,"xgraphc")
user.prompt()
eigraph(dbuf,"goodman")
user.prompt()
eigraph(dbuf,"xtfit")
user.prompt()
eigraph(dbuf,"xtfitg")
user.prompt()
eigraph(dbuf,"fit")
user.prompt()
eigraph(dbuf,"profile")
user.prompt()
eigraph(dbuf,"profileR")
user.prompt()
eigraph(dbuf,"postb")
user.prompt()
eigraph(dbuf,"postW")
user.prompt()
eigraph(dbuf,"post")
user.prompt()
message("Running beta with kern=E")
eigraph(dbuf,"betaB")
user.prompt()
eigraph(dbuf,"betaW")
user.prompt()
eigraph(dbuf,"beta")
user.prompt()
message("Running beta with kern=TN")
eigraph(dbuf,"beta",kern="TN")
user.prompt()
eigraph(dbuf,"results",kern="E")
user.prompt()
eigraph(dbuf,"movie")
user.prompt()
eigraph(dbuf,"movied")
user.prompt()
eigraph(dbuf,"lines")
user.prompt()
eigraph(dbuf,"bivar")
user.prompt()
eigraph(dbuf,"betabw")
user.prompt()
eigraph(dbuf,"biasb")
user.prompt()
eigraph(dbuf,"biasw")
user.prompt()
eigraph(dbuf,"bias")
user.prompt()
eigraph(dbuf,"boundxb")
user.prompt()
eigraph(dbuf,"boundxw")
user.prompt()
eigraph(dbuf,"boundx")
user.prompt()
message("Addition: dependences of beta's vs T,N")
eigraph(dbuf,"betast")
user.prompt()
eigraph(dbuf,"betasn")
user.prompt()
message("Addition: three-dimensional dependences of beta's vs X,T,N")
eigraph(dbuf,"betaxn")
user.prompt()
eigraph(dbuf,"betatn")
}













