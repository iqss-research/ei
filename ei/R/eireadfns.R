#Quantities of Interest

betaB <- function(ei.object){
ei.object$betab
}
betaW <- function(ei.object){
ei.object$betaw
}

sbetab <- function(ei.object){
ei.object$sbetab
}

sbetaw <- function(ei.object){
ei.object$sbetaw
}
phi <- function(ei.object){
ei.object$phi
}

psisims <- function(ei.object){
ei.object$psi
}

bounds <- function(ei.object){
bounds1(ei.object$x, ei.object$t, ei.object$n)
}

abounds <- function(ei.object){
	x <- ei.object$x
	t <- ei.object$t
	n <- ei.object$n
	bounds <- bounds1(x,t,n)
	omx <- 1-x
	Nb <- n*x
	Nw <- n*omx
	LAbetaB <- weighted.mean(bounds[,1],Nb)
	UAbetaB <- weighted.mean(bounds[,2], Nb)
	LAbetaW <- weighted.mean(bounds[,3], Nw)
	UAbetaW <- weighted.mean(bounds[,4], Nw)
	return(matrix(c(LAbetaB, UAbetaB, LAbetaW,UAbetaW), nrow=2))	
	}
	
aggs <- function(ei.object){
	ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
	x <- ei.object$x[ok]
	t <- ei.object$t[ok]
	n <- ei.object$n[ok]
	betab <- ei.object$betabs[ok,]
	betaw <- ei.object$betaws[ok,]
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
	
maggs <- function(ei.object){
	ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
	x <- ei.object$x[ok]
	t <- ei.object$t[ok]
	n <- ei.object$n[ok]
	betab <- ei.object$betabs[ok,]
	betaw <- ei.object$betaws[ok,]
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
	return(c(mean(Bbgg), mean(Bwgg), sd(Bbgg), sd(Bwgg)))
	}
	
VCaggs <- function(ei.object){
	ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
	x <- ei.object$x[ok]
	t <- ei.object$t[ok]
	n <- ei.object$n[ok]
	betab <- ei.object$betabs[ok,]
	betaw <- ei.object$betaws[ok,]
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
	
CI80b <- function(ei.object){
	ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
	betab <- ei.object$betabs[ok,]
	lwr <- vector(mode="numeric",length=length(x))
	upr <- vector(mode="numeric",length=length(x))
	lwr[ok] <- apply(betab, 1, function(x) quantile(x, probs=c(.1)))
	lwr[!ok] <- NA
	upr[ok] <-apply(betab, 1, function(x) quantile(x, probs=c(.9)))
	upr[!ok] <- NA
	return(cbind(lwr,upr))
	}
	
CI80w <- function(ei.object){
	ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
	betaw <- ei.object$betaws[ok,]
	lwr <- vector(mode="numeric",length=length(x))
	upr <- vector(mode="numeric",length=length(x))
	lwr[ok] <- apply(betaw, 1, function(x) quantile(x, probs=c(.1)))
	lwr[!ok] <- NA
	upr[ok] <-apply(betaw, 1, function(x) quantile(x, probs=c(.9)))
	upr[!ok] <- NA
	return(cbind(lwr,upr))
	}
	
eaggbias <- function(ei.object){
	x <- ei.object$x
	mbetab <- ei.object$betab
	mbetaw <- ei.object$betaw
	lm.b <- lm(mbetab ~ x)
	lm.w <- lm(mbetaw ~ x)
	output <- list(summary(lm.b)$coefficients,summary(lm.w)$coefficients)
	names(output) <- c("betaB", "betaW")
	return(output)
	}
	
goodman <- function(ei.object){
	x <- ei.object$x
	t <- ei.object$t
	lm.g <- lm(t ~ x)
	w <- 1-x
	lm.w <- lm(t ~ w)
	BetaW <- summary(lm.g)$coefficients[1,]
	BetaB <- summary(lm.w)$coefficients[1,]
	coefs <- rbind(BetaB, BetaW)
	return(coefs)
	}
