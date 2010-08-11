summary.ei <- function(ei1){
#Calculate maximum likelihood results in the scale of estimation
numb <- ei1$numb
covs <- as.logical(ifelse(diag(ei1$hessian)==0,0,1))
sdphi <- sqrt(diag((solve(ei1$hessian[covs,covs]))))
zbmiss <- ifelse(covs[6:(5+numb)]==FALSE,TRUE,FALSE)
zwmiss <- ifelse(covs[(6+numb):length(covs)]==FALSE, TRUE, FALSE)
names <- c("Bb0", "Bw0", "sigB", "sigW", "rho")
if(zbmiss==FALSE&zwmiss==FALSE){
sdphi <- c(sdphi[1:5], 0,0)
names <- c(names, "Zb", "Zw")
}
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
mle <- rbind(ei1$phi, sdphi)
colnames(mle) <- names
rownames(mle) <- c("","")

#Calculate untruncated psis
n <- length(ei1$betab)
BB <- mean(ei1$psi[,1:n])
BW <- mean(ei1$psi[,(n+1):(2*n)])
SB <- mean(ei1$psi[,((2*n)+1)])
SW <- mean(ei1$psi[,((2*n)+2)])
RHO <- mean(ei1$psi[,((2*n)+3)])
psiu <- t(matrix(c(BB,BW,SB,SW,RHO)))
colnames(psiu) <- c("BB","BW", "SB", "SW", "RHO")
rownames(psiu) <- c("")

#Calculate truncated psis
BB <- mean(ei1$betabs)
BW <- mean(ei1$betaws)
SB <- sd(as.vector(ei1$betabs))
SW <- sd(as.vector(ei1$betaws))
RHO <- cor(as.vector(ei1$betabs), as.vector(ei1$betaws))
psit <- t(matrix(c(BB,BW,SB,SW,RHO)))
colnames(psit) <- c("BB","BW", "SB", "SW", "RHO")
rownames(psit) <- c("")

#Aggregate Bounds
ab <- matrix(abounds(ei1), nrow=2)
rownames(ab) <- c("lower", "upper")
colnames(ab) <- c("betab", "betaw")

#Estimates of Aggregate Quantities of Interest
magg <- matrix(maggs(ei1), nrow=2)
rownames(magg) <- c("Bb", "Bw")
colnames(magg) <- c("mean", "sd")

output <- list(ei1$erho, ei1$esigma, ei1$ebeta, n, ei1$resamp, mle, psiu, psit, ab, magg)
names(output) <- c("Erho", "Esigma", "Ebeta", "N", "Resamp", "Maximum likelihood results in scale of estimation (and se's)", "Untruncated psi's", "Truncated psi's (ultimate scale)", "Aggregate Bounds", "Estimates of Aggregate Quantities of Interest")
class(output) <- "summary"
return(output)
}

print.summary <- function(sum.object, ...){
for (key in names(sum.object)){
	val <- sum.object[[key]]
	if ((is.character(val) || is.numeric(val)) && length(val) < 2){
		message(cat(key, "=", val))
	}
	else{
		message()
		message(key)
		print(val)
}
}
}
