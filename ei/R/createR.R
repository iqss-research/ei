createR <- function(sub, Rfun, bb, bw, sb,sw, rho){
out <- NULL
lower = cbind(-bb[sub]/sb, -bw[sub]/sw)
upper = cbind(-bb[sub]/sb+1/sb, -bw[sub]/sw+1/sw)
mean=c(0,0)
corr=matrix(c(1,rho,rho,1),nrow=2)
#pmvnorm

if (Rfun==1){
	for(i in 1:length(x[sub])){
		qi <- pmvnorm(lower=lower[i,], upper=upper[i,], mean=mean, corr=corr)
		qi <- ifelse(qi<0|qi==0, 1*10^-322,qi)
		out[i] <- log(qi)
		if(is.na(out[i])|abs(out[i]==Inf)) print("R not real")
		out[i] <- ifelse(is.na(out[i])|abs(out[i]==Inf), 999, out[i])
		}
	return(out)
	}

if (Rfun==2){
	for(i in 1:length(x[sub])){
		qi <- sadmvn(lower=lower[i,], upper=upper[i,], mean=mean, varcov=corr)
		qi <- ifelse(qi<0|qi==0, 1*10^-322,qi)
		out[i] <- log(qi)
		if(is.na(out[i])|abs(out[i]==Inf)) print("R not real")
		out[i] <- ifelse(is.na(out[i])|abs(out[i]==Inf), 999, out[i])
		}
	return(out)
	}
if (Rfun==3){
	fun <- function(x) dmvnorm(x,mean,corr)
	for (i in 1:length(x[sub])){
		qi <- adaptIntegrate(fun, lower=lower[i,], upper=upper[i,])$integral
		out[i] <- log(qi)
		if(is.na(out[i])|abs(out[i]==Inf)) print("R not real")
		out[i] <- ifelse(is.na(out[i])|abs(out[i]==Inf), 999, out[i])
		}
	return(out)
	}
if (Rfun==4){
	for(i in 1:length(x[sub])){
		qi <- mvnprd(A=upper[i,], B=lower[i,], BPD=c(rho,rho),INF=rep(2,2))$PROB
		qi <- ifelse(qi<0|qi==0, 1*10^-322,qi)
		out[i] <- log(qi)
		if(is.na(out[i])|abs(out[i]==Inf)) print("R not real")
		out[i] <- ifelse(is.na(out[i])|abs(out[i]==Inf), 999, out[i])
		}
	return(out)
	}
}	

