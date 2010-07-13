#Maximum likelihood estimation
library(msm)
library(mvtnorm)
library(tmvtnorm)

samp <- function(t,x,n, Zb, Zw, par, varcv, nsims, keep, numb, covs){
#a <- t(chol(varcv))
#mat <- matrix(rnorm(5*nsims), nrow=5)
#draw <- t(par + a%*%mat)
import1 <- NULL
draw <- rmvnorm(nsims, par[covs], varcv)
varcv2 <- solve(varcv)
phiv <- dmvnorm(draw, par[covs], varcv2, log=T)
zbmiss <- ifelse(covs[6:(5+numb)]==FALSE,TRUE,FALSE)
zwmiss <- ifelse(covs[(6+numb):length(covs)]==FALSE, TRUE, FALSE)
if(zbmiss==TRUE&zwmiss==FALSE){
	draw <- cbind(draw[,1:5], rep(1,nsims), draw[,(6+numb):length(par)])
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
import1[i] <- like(as.vector(draw[i,]), t, x, n, Zb, Zw, numb=1) - phiv[i]
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


ei <- function(t,x,n,Zb,Zw){
numb <- dim(Zb)[2]
numw <- dim(Zw)[2]
start <- c(0,0,-1.2,-1.2, 0, rep(0, numb+numw))
solution <- optim(start, like, y=t, x=x, n=n, Zb=Zb, Zw=Zw,numb=numb, method="BFGS", control=list(fnscale=-1), hessian=T) 

covs <- as.logical(ifelse(diag(solution$hessian)==0,0,1))
varcv <- solve(-solution$hessian[covs,covs])
varcv <- varcv/4


keep <- matrix(data=NA, ncol=(length(solution$par)))
while(dim(keep)[1]<100){
	keep <- samp(t,x,n, Zb, Zw, solution$par, varcv, 1000, keep, numb=numb, covs)
	print(dim(keep)[1])
	}

#for (i in 1:dim(keep)[1]){
#	keeplik <- like(as.vector(keep[i,]), data$t, data$x, data$n, #data$x, data$x)
#	}


keep <- keep[2:dim(keep)[1],]
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
bw <- psi[,(length(x)+1):length(data$x)*2]
sb <- psi[,(length(x)*2+1)]
sw <- psi[,(length(x)*2+2)]
rho <- psi[,(length(x)*2+3)]
omx <- 1-data$x
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
output <- list(solution$par, solution$hessian, psi, mbetab, mbetaw, sdbetab, sdbetaw)
names(output) <- c("phi", "hessian", "psi", "betab", "betaw", "sbetab", "sbetaw")
return(output)
}
