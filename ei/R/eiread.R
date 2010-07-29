#Quantities of Interest

abounds <- function(x,t,n){
	bounds <- bounds1(x,t,n)
	omx <- 1-x
	Nb <- n*x
	Nw <- n*omx
	LAbetaB <- weighted.mean(bounds[,1],Nb)
	UAbetaB <- weighted.mean(bounds[,2], Nb)
	LAbetaW <- weighted.mean(bounds[,3], Nw)
	UAbetaW <- weighted.mean(bounds[,4], Nw)
	return(c(LAbetaB, UAbetaB, LAbetaW,UAbetaW))	
	}
	
aggs <- function(x,t,n, betab, betaw){
	betab <- as.matrix(betab)
	betaw <- as.matrix(betaw)
	omx <- 1-x
	Nb <- n*x
	Nw <- n*omx
	Bbgg <- vector(mode="numeric", length=dim(betab)[2])
	for (i in 1:dim(betab)[2]){
		Bbgg[i] <- weighted.mean(betab[,i], Nb)
		}
	Bwgg <- vector(mode="numeric", length=dim(betaw)[2])
	for (i in 1:dim(betaw)[2]){
		Bwgg[i] <- weighted.mean(betaw[,i], Nw)
		}
	return(cbind(Bbgg, Bwgg))
	}
	
maggs <- function(x,t,n, betab, betaw){
	omx <- 1-x
	Nb <- n*x
	Nw <- n*omx
	betab <- as.matrix(betab)
	betaw <- as.matrix(betaw)
	Bbgg <- vector(mode="numeric", length=dim(betab)[2])
	for (i in 1:dim(betab)[2]){
		Bbgg[i] <- weighted.mean(betab[,i], Nb)
		}
	Bwgg <- vector(mode="numeric", length=dim(betaw)[2])
	for (i in 1:dim(betaw)[2]){
		Bwgg[i] <- weighted.mean(betaw[,i], Nw)
		}
	return(c(mean(Bbgg), mean(Bwgg), sd(Bbgg), sd(Bwgg)))
	}
	
VCaggs <- function(x,t,n,betab, betaw){
	betab <- as.matrix(betab)
	betaw <- as.matrix(betaw)
	omx <- 1-x
	Nb <- n*x
	Nw <- n*omx
	Bbgg <- vector(mode="numeric", length=dim(betab)[2])
	for (i in 1:dim(betab)[2]){
		Bbgg[i] <- weighted.mean(betab[,i], Nb)
		}
	Bwgg <- vector(mode="numeric", length=dim(betaw)[2])
	for (i in 1:dim(betaw)[2]){
		Bwgg[i] <- weighted.mean(betaw[,i], Nw)
		}
	vc <- matrix(c(var(Bbgg), cov(Bbgg, Bwgg), cov(Bbgg, Bwgg), var(Bwgg, Bwgg)), nrow=2)
	return(vc)
	}
	
CI80b <- function(betab){
	lwr <- apply(betab, 1, function(x) quantile(x, probs=c(.1)))
	upr <-apply(betab, 1, function(x) quantile(x, probs=c(.9)))
	return(cbind(lwr,upr))
	}
	
CI80w <- function(betaw){
	lwr <- apply(betaw, 1, function(x) quantile(x, probs=c(.1)))
	upr <-apply(betaw, 1, function(x) quantile(x, probs=c(.9)))
	return(cbind(lwr,upr))
	}
	
eaggbias <- function(x, mbetab, mbetaw){
	lm.b <- lm(mbetab ~ x)
	lm.w <- lm(mbetaw ~ x)
	return(rbind(summary(lm.b)$coefficients, summary(lm.w)$coefficients))
	}
	
goodman <- function(t, x){
	lm.g <- lm(t ~ x)
	w <- 1-x
	lm.w <- lm(t ~ w)
	BetaW <- summary(lm.g)$coefficients[1,]
	BetaB <- summary(lm.w)$coefficients[1,]
	coefs <- rbind(BetaB, BetaW)
	return(coefs)
	}
