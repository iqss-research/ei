#Maximum likelihood estimation
library(msm)
library(mvtnorm)
library(tmvtnorm)

samp <- function(t,x,n, Zb, Zw, par, varcv, nsims, keep, numb, covs, erho, esigma, ebeta, ealphab, ealphaw){
#a <- t(chol(varcv))
#mat <- matrix(rnorm(5*nsims), nrow=5)
#draw <- t(par + a%*%mat)
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
#ymu <- draw[i,] - par
#phiv = -(t(as.matrix(ymu))%*%solve(varcv)%*%as.matrix(ymu))/2
import1[i] <- like(as.vector(draw[i,]), t, x, n, Zb, Zw, numb=numb, erho, esigma, ebeta, ealphab, ealphaw) - phiv[i]
}
lnir <- import1-max(import1[1:nsims])
ir <- exp(lnir)
tst <- ir[1:nsims]>runif(nsims,0,1)
keep <- rbind(keep, draw[tst,])
return(keep)
}

## commented out
## Zb <- as.matrix(data$x)
## Zw <- as.matrix(data$x)

#ei1 <- ei(t,x,n,Zb,Zw,erho=.5,esigma=.5,ebeta=0,ealphab=NA,ealphaw=NA)
ei <- function(t,x,n,Zb,Zw, erho, esigma, ebeta, ealphab, ealphaw){
Zb <- as.matrix(Zb)
Zw <- as.matrix(Zw)
numb <- dim(Zb)[2]
numw <- dim(Zw)[2]
start <- c(0,0,-1.2,-1.2, 0, rep(0, numb+numw))
solution <- optim(start, like, y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb, erho=erho, esigma=esigma, ebeta=ebeta, ealphab =ealphab, ealphaw=ealphaw, method="BFGS", control=list(fnscale=-1), hessian=T) 
covs <- as.logical(ifelse(diag(solution$hessian)==0,0,1))
varcv <- -solution$hessian[covs,covs]
#varcv <- varcv


keep <- matrix(data=NA, ncol=(length(solution$par)))
resamp <- 0
while(dim(keep)[1]<100){
	keep <- samp(t,x,n, Zb, Zw, solution$par, varcv, 1000, keep, numb=numb, covs, erho, esigma, ebeta, ealphab, ealphaw)
	resamp = resamp + 1
	}

#for (i in 1:dim(keep)[1]){
#	keeplik <- like(as.vector(keep[i,]), data$t, data$x, data$n, #data$x, data$x)
#	}

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
for (i in 1:dim(keep)[1]){
sig2 <- sb[i]^2*x^2 + sw[i]^2*omx^2 + sbw[i]*2*x*omx
omega <- sb[i]^2*x + sbw[i]*omx
eps <- t - (bb[i,])*x - (bw[i,])*omx
mbb <- bb[i,] + omega/sig2*eps
vbb <- sb[i]^2 - (omega^2)/sig2
s <- sqrt(vbb)
bounds <- bounds1(x,t,n)
out<- NULL
for(j in 1:length(x)){
out[j] <- rtnorm(1, mean=mbb[j], sd=s[j], lower=bounds[j,1], upper=bounds[j,2])
}
betab[,i] = out
}
omx <- 1 - x
for (j in 1:length(x)){
	betabs <- betab[j,]
	betaw[j,] <- t[j]/omx[j]-betabs*x[j]/omx[j]
	}

mbetab <- apply(betab,1,mean)
mbetaw <- apply(betaw,1,mean)
sdbetab <- apply(betab,1,sd)
sdbetaw <- apply(betaw,1,sd)
output <- list(solution$par, solution$hessian, psi, mbetab, mbetaw, sdbetab, sdbetaw, betab, betaw, resamp, erho, esigma, ebeta, ealphab, ealphaw, numb, x, t, n)
names(output) <- c("phi", "hessian", "psi", "betab", "betaw", "sbetab", "sbetaw", "betabs", "betaws", "resamp", "erho", "esigma", "ebeta", "ealphab", "ealphaw", "numb", "x", "t", "n")
class(output) <- "ei"
return(output)
}

summary.ei <- function(ei1){
#Calculate maximum likelihood results in the scale of estimation
numb <- ei1$numb
covs <- as.logical(ifelse(diag(ei1$hessian)==0,0,1))
sdphi <- diag(sqrt(solve(-ei1$hessian[covs,covs])))
zbmiss <- ifelse(covs[6:(5+numb)]==FALSE,TRUE,FALSE)
zwmiss <- ifelse(covs[(6+numb):length(covs)]==FALSE, TRUE, FALSE)
names <- c("Bb0", "Bw0", "sigB", "sigW", "rho")
if(zbmiss==TRUE&zwmiss==FALSE){
	sdphi <- c(sdphi[1:5], 0, sdphi[(5+numb):sum(covs)])
	numw <- length(ei1$phi) - (5+numb)
	wname <- NULL
	for (i in 1:numw){
	wname[i] = paste("Zw",(i-1), sep="")
	}
	names <- c(names, "Zb0", wname)
	}
if(zbmiss==FALSE&zwmiss==TRUE){
	sdphi <- c(sdphi, 0)
	bname <- NULL
	for (i in 1:numb){
	bname[i] = paste("Zb",(i-1), sep="")
	}
	names <- c(names, bname, "Zw0")
	}
if(zbmiss==TRUE&zwmiss==TRUE){
	sdphi <- c(sdphi, 0, 0)
	names <- c(names, "Zb0", "Zw0")
	}
mle <- as.data.frame(rbind(ei1$phi, sdphi))
names(mle) <- names


#Calculate untruncated psis
n <- length(ei1$betab)
BB <- mean(ei1$psi[,1:n])
BW <- mean(ei1$psi[,(n+1):(2*n)])
SB <- mean(ei1$psi[,((2*n)+1)])
SW <- mean(ei1$psi[,((2*n)+2)])
RHO <- mean(ei1$psi[,((2*n)+3)])
psiu <- as.data.frame(t(c(BB,BW,SB,SW,RHO)))
names(psiu) <- c("BB","BW", "SB", "SW", "RHO")

#Calculate truncated psis
BB <- mean(ei1$betabs)
BW <- mean(ei1$betaws)
SB <- sd(as.vector(ei1$betabs))
SW <- sd(as.vector(ei1$betaws))
RHO <- cor(as.vector(ei1$betabs), as.vector(ei1$betaws))
psit <- as.data.frame(t(c(BB,BW,SB,SW,RHO)))
names(psit) <- c("BB","BW", "SB", "SW", "RHO")

#Aggregate Bounds
ab <- matrix(abounds(ei1$x, ei1$t, ei1$n), nrow=2)


#Estimates of Aggregate Quantities of Interest
magg <- matrix(maggs(ei1$x,ei1$t,ei1$n, ei1$betabs, ei1$betaws), nrow=2)

output <- list(ei1$erho, ei1$esigma, ei1$ebeta, n, ei1$resamp, mle, psiu, psit, ab, magg)
names(output) <- c("Erho", "Esigma", "Ebeta", "N", "Resamp", "MLE", "psiu", "psit", "ab", "magg")
return(output)
}
