library(msm)
library(mvtnorm)
library(tmvtnorm)

## repar <- function(Bb0, Bw0, sb0, sw0, rho0){
## 	sb=exp(sb0)
## 	sw=exp(sw0)
## 	bb=Bb0*(.25+sb^2) + .5
## 	bw=Bw0*(.25+sw^2) + .5
## 	rho=(exp(2*rho0)-1)/(exp(2*rho0) +1)
## 	return(c(bb,bw,sb,sw,rho))
## 	}

repar <- function(Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, Zw){
	sb=exp(sb0)
	sw=exp(sw0)
	bb=Bb0*(.25+sb^2) + .5 + as.matrix(apply(Zb,2, function(x) x - mean(x)))%*%as.matrix(Bb0v)
	bw=Bw0*(.25+sw^2) + .5 + as.matrix(apply(Zw,2, function(x) x - mean(x)))%*%as.matrix(Bw0v)
	rho=(exp(2*rho0)-1)/(exp(2*rho0) +1)
	return(c(t(bb), t(bw), sb, sw, rho))
	}
	
## #exvar <- function(x,y,n,bb, bw, sb,sw, rho){
## 	sigb2 <- sb^2
## 	sigw2 <- sw^2
## 	sigbw = rho*sb*sw
## 	omx <- 1-x
## 	mu = bb*x + bw*omx
## 	epsilon = y - mu
## 	#s2 = sigw2 + (2*sigbw-2*sigw2)*x + (sigb2 + sigw2 - 2*sigbw)*(x^2)
## 	s2 = sigb2*(x^2) + sigw2*(omx^2) + 2*sigbw*x*omx
## 	omega = sigb2*x + sigbw*omx
## 	ebb = bb + (omega/s2)*epsilon
## 	vbb = sigb2 - (omega^2)/s2
## 	vbb = ifelse(vbb<1*10^-322, 1*10^-322, vbb)
## 	return(cbind(mu, s2, epsilon, omega, ebb,vbb))
## 	}

#numb <- 2
like <- function(param, y, x, n, Zb, Zw, numb, erho, esigma, ebeta, ealphab, ealphaw){
	Bb0 <- param[1]
	Bw0 <- param[2]
	sb0 <- param[3]
	sw0 <- param[4]
	rho0 <- param[5]
	Bb0v <- param[6:(5+numb)]
	Bw0v <- param[(numb+6):length(param)]
	sb=exp(sb0)
	sw=exp(sw0)
	Zb <- as.matrix(Zb)
	Zw <- as.matrix(Zw)
	bb=Bb0*(.25+sb^2) + .5 + as.matrix(apply(Zb,2, function(x) x - mean(x)))%*%as.matrix(Bb0v)
	bw=Bw0*(.25+sw^2) + .5 + as.matrix(apply(Zw,2, function(x) x - mean(x)))%*%as.matrix(Bw0v)
	rho=(exp(2*rho0)-1)/(exp(2*rho0) +1)
	sigb2 <- sb^2
	sigw2 <- sw^2
	sigbw = rho*sb*sw
print(c(bb[1],bw[1],sb[1],sw[1],rho[1]))
#Compute likelihood for different categories
homoindx <- ifelse(x==0, 1, 0)
homoindx <- ifelse(x==1, 2, homoindx)
enumtol=.0001
cT0 <- y<enumtol & homoindx==0
cT1 <- y>(1-enumtol) & homoindx==0
ok <- ifelse(homoindx==0 & cT0==0 & cT1==0,T, F)

#0<T<1, 0<X<1

	omx <- 1-x
	mu = bb*x + bw*omx
	epsilon = y - mu
	s2 = sigb2*(x^2) + sigw2*(omx^2) + 2*sigbw*x*omx
	omega = sigb2*x + sigbw*omx
	ebb = bb + (omega/s2)*epsilon
	vbb = sigb2 - (omega^2)/s2
	vbb = ifelse(vbb<1*10^-322, .0001, vbb)
	bounds <- bounds1(x, y, n)
	s <- sqrt(vbb)
	res <- NULL
	#res[ok] <- log(pnorm(bounds[ok,2], mean=ebb[ok], sd=s[ok]) - #pnorm(bounds[ok,1], mean=ebb[ok], sd=s[ok]))
b.s = (bounds[ok,][,2]-ebb[ok])/s[ok]
as = (bounds[ok,][,1]-ebb[ok])/s[ok]
res[ok] <- log(pnorm(as, lower.tail=F) - pnorm(b.s, lower.tail=F))

	R <- NULL
	bs <- as.matrix(cbind(bb, bw))
for (i in 1:length(x[ok])){
R[ok][i] <-log(pmvnorm(lower=c(-bb[ok][i]/sb,-bw[ok][i]/sw), upper=c(-bb[ok][i]/sb+1/sb,-bw[ok][i]/sw+1/sw), mean=c(0, 0), corr=matrix(c(1,rho,rho,1), nrow=2)))
}
	#R[ok] <- apply(bs[ok,], 1, function (x) log(pmvnorm(lower=
#c(0,0), upper=c(1,1), mean=as.vector(x), sigma=matrix(c(sigb2, #sigbw,sigbw,sigw2),nrow=2), algorithm=Miwa())))
	#for(i in 1:length(x)){
	#R[i] <- log(pmvnorm(lower=c(0,0), upper=c(1,1), #mean=c(bb[i], bw[i]), sigma=matrix(c(sigb2,sigbw,sigbw,sigw2), nrow=2)))
	#}
	

	lliki  <- -.5*(log(s2[ok])+epsilon[ok]^2/s2[ok])
	llik.het <- -.5*sum((log(s2[ok])+(epsilon[ok]^2)/(s2[ok]))) 
	llik.het <- llik.het + sum(res[ok]) - sum(R[ok])

#Homogenous precincts

#X=0
wh <- homoindx==1
llik.wh=0
if (sum(wh)>0){
epsilon= y[wh]-bw[wh]
llik.whi = -.5*(log(sigw2)+(epsilon^2)/(sigw2))
llik.wh=-.5*sum((log(sigw2)+(epsilon^2)/(sigw2)))
bnds=cbind(rep(0,sum(wh)),rep(1,sum(wh)))
Ebb = bb[wh]+rho*(sb/sw)*epsilon
Vbb = sigb2*(1-rho^2)
s <- sqrt(Vbb)
res <- log(pnorm(bnds[,2], mean=Ebb, sd=s) - pnorm(bnds[,1], mean=Ebb, sd=s))
if(sum(wh)==1){
R <- log(pmvnorm(lower=c(0,0), upper=c(1,1), mean=bs[wh,], sigma=matrix(c(sigb2, sigbw,sigbw,sigw2),nrow=2)))
}
if(sum(wh)>1){
R <- apply(bs[wh,], 1, function (x) log(pmvnorm(lower=c(0,0), upper=c(1,1), mean=as.vector(x), sigma=matrix(c(sigb2, sigbw,sigbw,sigw2),nrow=2))))
}
llik.wh = llik.wh + sum(res)-sum(R)

}

#X=1
bl <- homoindx==2
llik.bl=0
if (sum(bl)>0){
epsilon = y[bl] - bb[bl]
llik.bl = -.5*sum((log(sigb2)+(epsilon^2)/(sigb2)))
bnds = cbind(rep(0, sum(bl)),rep(1,sum(bl)))
Ebb=bw[bl] + rho*(sw/sb)*epsilon
Vbb=sigw2*(1-rho^2)
s <- sqrt(Vbb)
res <- log(pnorm(bnds[,2], mean=Ebb, sd=s) - pnorm(bnds[,1], mean=Ebb, sd=s))
if(sum(bl)==1){
R <- log(pmvnorm(lower=c(0,0), upper=c(1,1), mean=bs[bl,], sigma=matrix(c(sigb2, sigbw,sigbw,sigw2),nrow=2)))
}
if(sum(bl)>1){
R <- apply(bs[bl,], 1, function (x) log(pmvnorm(lower=c(0,0), upper=c(1,1), mean=as.vector(x), sigma=matrix(c(sigb2, sigbw,sigbw,sigw2),nrow=2))))
}
llik.bl = llik.bl + sum(res)-sum(R)
}

#T=0, 0<X<1
llik.cT0=0
if(sum(cT0)>0){
if(sum(cT0)==1){
first = log(dmvnorm(c(0,0),mean=bs[cT0,], sigma=matrix(c(sigb2,sigbw,sigbw,sigw2),nrow=2)))

second = log(pmvnorm(lower=c(0,0), upper=c(1,1), mean=bs[cT0,], sigma=matrix(c(sigb2, sigbw,sigbw,sigw2),nrow=2)))
llik.cT0=sum(first)-sum(second)
}
else{
first = apply(bs[cT0,], 1, function (x) log(dmvnorm(c(0,0),mean=as.vector(x), sigma=matrix(c(sigb2,sigbw,sigbw,sigw2),nrow=2))))

second = apply(bs[cT0,], 1, function (x) log(pmvnorm(lower=c(0,0), upper=c(1,1), mean=as.vector(x), sigma=matrix(c(sigb2, sigbw,sigbw,sigw2),nrow=2))))
llik.cT0=sum(first)-sum(second)
}
}

#T=1, 0<X<1
llik.cT1=0
if(sum(cT1)>0){
if(sum(cT1)==1){
first = log(dmvnorm(c(1,1),mean=bs[cT1,], sigma=matrix(c(sigb2,sigbw,sigbw,sigw2),nrow=2)))
bb.cT1 = bs[cT1,][1]
bw.cT1 = bs[cT1,][2]
second = log(pmvnorm(lower=c(-bb.cT1/sb,-bw.cT1/sw), upper=c(-bb.cT1/sb + 1/sb, -bw.cT1/sw + 1/sw), mean=c(0,0), corr=matrix(c(1,rho,rho,1), nrow=2)))
llik.cT1=sum(first)-sum(second)
}
if(sum(cT1)>1){
first = apply(as.matrix(bs[cT1,]), 1, function (x) log(dmvnorm(c(1,1),mean=as.vector(x), sigma=matrix(c(sigb2,sigbw,sigbw,sigw2),nrow=2))))

second = apply(as.matrix(bs[cT1,]), 1, function (x) log(pmvnorm(lower=c(0,0), upper=c(1,1), mean=as.vector(x), sigma=matrix(c(sigb2, sigbw,sigbw,sigw2),nrow=2))))
llik.cT1=sum(first)-sum(second)
}
}
llik=llik.het + llik.bl + llik.wh + llik.cT0 + llik.cT1


#Priors
	prior=0
	lpdfnorm = log(dnorm(rho0,0,sd=erho))
	if (esigma>0) prior = prior-(1/(2*esigma^2))*(sigb2+sigw2);
	if(erho>0) prior = prior +lpdfnorm;
	if(ebeta>0 & (mean(bb)<0)) prior = prior -.5*((mean(bb)^2)/ebeta);
	if(ebeta>0 & mean(bb)>1) prior = prior -.5*((mean(bb)-1)^2/ebeta);
	if(ebeta>0 & mean(bw)<0) prior = prior -.5*((mean(bw)^2)/ebeta);
	if(ebeta>0 & mean(bw)>1) prior = prior -.5*((mean(bw)-1)^2/ebeta);
	if(sum(is.na(ealphab))==0) prior=prior + sum(dmvnorm(Bb0v, ealphab[,1], sigma=diag(ealphab[,2]^2), log=T));
	if(sum(is.na(ealphaw))==0) prior=prior + sum(dmvnorm(Bw0v, ealphaw[,1], sigma=diag(ealphaw[,2]), log=T));
	llik = llik + prior
	return(-llik)
      }


#Alternative res
#s = sqrt(vbb)
#b.s = (bounds[,2]-ebb)/s
#as = (bounds[,1]-ebb)/s
#log(pnorm(b.s) - pnorm(as))

#for (i in 1:75) {
#res[i] <- log(pmvnorm(lower=c(as[i]), upper=c(b.s[i]), mean=c(0), #sigma=as.matrix(1)))
#}

#for (i in 1:75) {
#res[i] <- log(pmvnorm(lower=c(bounds[,1][i]), upper=c(bounds[,2][i]), #mean=ebb[i], sigma=as.matrix(vbb[i])))
#}
#Alternative R

#R <- ifelse(rho==0,
		#log(pmvnorm(lower=c(0), upper=c(1), mean=bb, #sigma=as.matrix(sigb2))) + log(pmvnorm(lower=c(0), upper=c(1), #mean=bw, sigma=as.matrix(sigw2))), log(pmvnorm(lower=c(0,0), #upper=c(1,1), mean=c(bb, bw), sigma=matrix(c(sigb2,sigbw,sigbw,sigw2), nrow=2))))
		
