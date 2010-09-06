#Maximum likelihood estimation
library(msm)
library(mvtnorm)
library(tmvtnorm)
library(ucminf)
#ei1 <- ei(t,x,n,Zb,Zw,erho=.5,esigma=.5,ebeta=0,ealphab=NA,ealphaw=NA, truth=NA)
ei <- function(t,x,n,Zb,Zw, erho=.5, esigma=.5, ebeta=0, ealphab=NA, ealphaw=NA, truth=NA){
Zb <- as.matrix(Zb)
Zw <- as.matrix(Zw)
if (dim(Zb)[1]==1 & Zb==1) Zb <- as.matrix(rep(1,length(x)))
if (dim(Zw)[1]==1 & Zw==1) Zw <- as.matrix(rep(1,length(x)))
numb <- dim(Zb)[2]
numw <- dim(Zw)[2]
start <- c(0,0,-1.2,-1.2, 0, rep(0, numb+numw))
message("Maximizing likelihood")
solution <- ucminf(start, like, y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb, erho=erho, esigma=esigma, ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw, hessian=3)
print(solution$par)
print(solution$convergence) 
#solution <- genoud(like, y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb, erho=erho, esigma=esigma, #ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw, nvars=5, starting.values=start)
#solution <- maxLik(like, y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb, erho=erho, esigma=esigma, #ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw,start=start)
#solution <- subplex(start, like, y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb, erho=erho, esigma=esigma, #ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw)
#solution <- nlminb(start, like,y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb, erho=erho, #esigma=esigma, #ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw, hessian=T)

covs <- as.logical(ifelse(diag(solution$hessian)==0,0,1))
varcv <- solution$hessian[covs,covs]
#varcv <- varcv

message("Importance Sampling..")
keep <- matrix(data=NA, ncol=(length(solution$par)))
resamp <- 0
while(dim(keep)[1]<100){
	keep <- samp(t,x,n, Zb, Zw, solution$par, varcv, 1000, keep, numb=numb, covs, erho, esigma, ebeta, ealphab, ealphaw)
	resamp = resamp + 1
	}


keep <- keep[2:100,]
mu <- keep[,1:2]
sd <- keep[,3:4]
rho <- keep[,5]
Bb0v <- keep[,6:(5+numb)]
Bw0v <- keep[,(6+numb):length(solution$par)]
sd[,1] <- exp(sd[,1])
sd[,2] <- exp(sd[,2])
Zb <- as.matrix(Zb)
Zw <- as.matrix(Zw)
mu1 <- mu[,1]*(.25 + sd[,1]^2) + .5 + t(as.matrix(apply(Zb,2, function (x) x - mean(x)))%*%t(as.matrix(Bb0v)))
mu2 <- mu[,2]*(.25 + sd[,2]^2) + .5 + t(as.matrix(apply(Zw,2, function (x) x - mean(x)))%*%t(Bw0v))
#phin <- dmvnorm(psi, par, log=T)
rho <- (exp(2*rho)-1)/(exp(2*rho) +1)
psi <- cbind(mu1, mu2, sd, rho)
bb <- psi[,1:length(x)]	
bw <- psi[,(length(x)+1):(length(x)*2)]
sb <- psi[,(length(x)*2+1)]
sw <- psi[,(length(x)*2+2)]
rho <- psi[,(length(x)*2+3)]
omx <- 1-x
sbw <- rho*sb*sw
betab <- matrix(nrow=length(x),ncol=dim(keep)[1])
betaw <- matrix(nrow=length(x),ncol=dim(keep)[1])

homoindx <- ifelse(x==0, 1, 0)
homoindx <- ifelse(x==1, 2, homoindx)
enumtol=.0001
cT0 <- t<enumtol & homoindx==0
cT1 <- t>(1-enumtol) & homoindx==0
ok <- ifelse(homoindx==0 & cT0==0 & cT1==0,T, F)
wh <- homoindx==1
bl <- homoindx==2


for (i in 1:dim(keep)[1]){
sig2 <- sb[i]^2*x^2 + sw[i]^2*omx^2 + sbw[i]*2*x*omx
omega <- sb[i]^2*x + sbw[i]*omx
eps <- t - (bb[i,])*x - (bw[i,])*omx
mbb <- bb[i,] + omega/sig2*eps
vbb <- sb[i]^2 - (omega^2)/sig2
s <- sqrt(vbb)
bounds <- bounds1(x,t,n)
out<- NULL
for(j in 1:length(x[ok])){
out[ok][j] <- rtnorm(1, mean=mbb[ok][j], sd=s[ok][j], lower=bounds[ok,][j,1], upper=bounds[ok,][j,2])
}
out[wh] <- NA
out[bl] <- t[bl]
out[cT1] <- bounds[cT1,1]
out[cT0] <- bounds[cT0,1]

betab[,i] = out
}

omx <- 1 - x
for (j in 1:length(x[ok])){
	betabs <- betab[ok,][j,]
	betaw[ok,][j,] <- t[ok][j]/omx[ok][j]-betabs*x[ok][j]/omx[ok][j]
	betaw[wh,] <- rep(1, dim(keep)[1])*t[wh]
	betaw[bl,] <- rep(NA, dim(keep)[1])
	betaw[cT1,] <- rep(bounds[cT1,3],dim(keep)[1])
	betaw[cT0,] <- rep(bounds[cT0,3], dim(keep)[1])
	}

mbetab <- apply(betab,1,mean)
mbetaw <- apply(betaw,1,mean)
sdbetab <- apply(betab,1,sd)
sdbetaw <- apply(betaw,1,sd)
output <- list(solution$par, solution$hessian, psi, mbetab, mbetaw, sdbetab, sdbetaw, betab, betaw, resamp, erho, esigma, ebeta, ealphab, ealphaw, numb, x, t, n, Zb, Zw, truth)
names(output) <- c("phi", "hessian", "psi", "betab", "betaw", "sbetab", "sbetaw", "betabs", "betaws", "resamp", "erho", "esigma", "ebeta", "ealphab", "ealphaw", "numb", "x", "t", "n", "Zb", "Zw", "truth")
class(output) <- "ei"
return(output)
}


samp <- function(t,x,n, Zb, Zw, par, varcv, nsims, keep, numb, covs, erho, esigma, ebeta, ealphab, ealphaw){
import1 <- NULL
varcv2 <- solve(varcv)/4
draw <- rmvnorm(nsims, par[covs], varcv2)
varcv3 <- solve(varcv2)
phiv <- dmvnorm(draw, par[covs], varcv3, log=T)
zbmiss <- ifelse(covs[6:(5+numb)]==FALSE,TRUE,FALSE)
zwmiss <- ifelse(covs[(6+numb):length(covs)]==FALSE, TRUE, FALSE)
if(zbmiss==TRUE&zwmiss==FALSE){
	draw <- cbind(draw[,1:5], rep(1,nsims), draw[,(5+numb):sum(covs)])
	}
if(zbmiss==FALSE&zwmiss==TRUE){
	draw <- cbind(draw, rep(1,nsims))
	}
if(zbmiss==TRUE&zwmiss==TRUE){
	draw <- cbind(draw, rep(1,nsims), rep(1,nsims))
	}
for(i in 1:nsims){
import1[i] <- -like(as.vector(draw[i,]), t, x, n, Zb, Zw, numb=numb, erho, esigma, ebeta, ealphab, ealphaw) - phiv[i]
}
lnir <- import1-max(import1[1:nsims])
ir <- exp(lnir)
print(mean(ir))
tst <- ir[1:nsims]>runif(nsims,0,1)
keep <- rbind(keep, draw[tst,])
return(keep)
}


print.ei <- function(ei.object, ...){
message("resamp")
message(ei.object$resamp)
message("Maximum likelihood results in scale of estimation")
print(ei.object$phi)
message("Use summary() and eiread() for more quantities of interest.")
}
