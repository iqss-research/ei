#Graphs
library(ellipse)
library(plotrix)
library(MASS)

#Tomography plot

tomog <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
bounds <- bounds1(x,t,n)
bbounds <- cbind(bounds[,1	], bounds[,2])
wbounds <- cbind(bounds[,4], bounds[,3])
n <- dim(bounds)[1]
plot(c(100,200), xlim=c(0,1), ylim=c(0,1), col="white", ylab="betaW", xlab="betaB", xaxs="i",yaxs="i", main="Tomography Plot")
for(i in 1:n){
	lines(bbounds[i,], wbounds[i,], col="yellow")
	}
}

tomogd <- function(x,t,n,title){
bounds <- bounds1(x,t,n)
bbounds <- cbind(bounds[,1	], bounds[,2])
wbounds <- cbind(bounds[,4], bounds[,3])
n <- dim(bounds)[1]
plot(c(100,200), xlim=c(0,1), ylim=c(0,1), col="white", ylab="betaW", xlab="betaB", xaxs="i",yaxs="i", main=title)
for(i in 1:n){
	lines(bbounds[i,], wbounds[i,], col="yellow")
	}
}

#Tomography plot with ML contours
tomogl <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
Zb <- ei.object$Zb
Zw <- ei.object$Zw
phi <- ei.object$phi
tomogd(x,t,n, "Tomography Plot with ML Contours")
numb <- dim(Zb)[2]
numw <- dim(Zw)[2]
Bb0 <- phi[1]
Bw0 <- phi[2]
sb0 <- phi[3]
sw0 <- phi[4]
rho0 <- phi[5]
Bb0v <- phi[6:(5+numb)]
Bw0v <- phi[(6+numb):length(phi)]
vars <- repar(Bb0,Bw0,sb0,sw0,rho0, Bb0v, Bw0v, Zb, Zw)
bb <- vars[1:length(x)]
bw <- vars[(length(x)+1):(2*length(x))]
sb <- vars[2*length(x)+1]
sw <- vars[2*length(x)+2]
rho <- vars[2*length(x)+3]
tomog3 <- function(bb,bw,sb,sw,rho){
	lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.914), col="blue",lwd=4)
	lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.35), col="red",lwd=4)
	points(mean(bb),mean(bw),col="pink",  pch=15)
	}
#tomog4 <- function(bb,bw,sb,sw,rho){
#	lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), #scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.914), #col="blue",lwd=1, lty=3)
#	lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), #scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.35), #col="red",lwd=1, lty=3)
#	}
#or
#for (i in 1:length(x)){
	#points(mean(bb[i]), mean(bw[i]), col="red", #pch=19, cex=.1)
#	tomog3(bb[i], bw[i], sb, sw, rho)
	#}
tomog3(bb,bw,sb,sw,rho)
}

tomog3 <- function(bb,bw,sb,sw,rho){
	lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.914), col="blue",lwd=4)
	lines(ellipse(matrix(c(1,rho,rho,1),nrow=2), scale=c(sb,sw),centre=c(mean(bb),mean(bw)),level=.35), col="red",lwd=4)
	points(mean(bb),mean(bw),col="pink",  pch=15)
	}


#Tomography plot with 80% CIs
tomog80CI <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
betabs <- ei.object$betabs
betaws <- ei.object$betaws
tomogd(x,t,n,"Tomography Plot with 80% CIs")
betabcd <- apply(betabs,1,function(x) quantile(x, probs=c(.1,.9)))
betawcd <- apply(betaws,1,function (x) quantile(x,probs=c(.1,.9)))
n <- dim(betabcd)[2]
for(i in 1:n){
	lines(betabcd[,i], sort(betawcd[,i],decreasing=T), col="red", lwd=3)
	}
}
#Tomography plot with 95% CIs
tomog95CI <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
betabs <- ei.object$betabs
betaws <- ei.object$betaws
tomogd(x,t,n,"Tomography Plot with 95% CIs")
betabcd <- apply(betabs,1,function(x) quantile(x, probs=c(.025,.975)))
betawcd <- apply(betaws,1,function (x) quantile(x,probs=c(.025,.975)))
n <- dim(betabcd)[2]
for(i in 1:n){
	lines(betabcd[,i], sort(betawcd[,i], decreasing=T), col="red", lwd=3)
	}
}

#TomogE -- Tomography plot with mean posterior betabs and betaws
tomogE <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
betabs <- ei.object$betabs
betaws <- ei.object$betaws
tomogd(x,t,n,"Tomography Plot with Mean Posterior Betabs and Betaws")
betabm <- apply(betabs,1,mean)
betawm <- apply(betaws,1,mean)
points(betabm, betawm, col="red", pch=19)
}

#TomogP -- Tomography plot with contours based on mean posterior psi
#tomogd(data$x,data$t, data$n, "Tomography with simulated contours based on mean posterior")
#vars <- apply(psi,2,mean)
#bb <- vars[1]
#bw <- vars[2]
#sb <- vars[3]
#sw <- vars[4]
#rho <- vars[5]
#tomog(bb,bw,sb,sw,rho)
#dev.off()


tomogP2 <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
psi <- ei.object$psi
tomogd(x,t,n,"Tomography Plot with Contours Based on Posterior")
bbp <- psi[,1:length(x)]
bwp <- psi[,(length(x)+1):(2*length(x))]
sbp <- psi[,2*length(x)+1]
swp <- psi[,2*length(x)+2]
rhop <- psi[,2*length(x)+3]
points(mean(bbp), mean(bwp), col="red", pch=19)
tomog3(bbp,bwp,mean(sbp),mean(swp),mean(rhop))
#or
#points(mean(bbp), mean(bwp), col="red", pch=19)
#for (i in 1:length(x)){
#	points(mean(bbp[,i]), mean(bwp[,i]), col="red", #pch=19, cex=.1)
#	tomog4(bbp[,i], bwp[,i], mean(sbp), mean(swp), mean #(rhop))
#	}
}


#Density plots
#tomogd(data$x,data$t,data$n, "Tomography with contours #from posterior")
#bivn.kde <- kde2d(psi[,1], psi[,2], n=100)
#contour(bivn.kde, add=T, nlevels=7, drawlabels=F)
#points(mean(psi[,1]), mean(psi[,2]), col="red", pch=15)

#Density of betab
betabd <- function(ei.object){
betabs <- ei.object$betabs
betabm <- apply(betabs,1,mean)
plot(density(betabm), xlim=c(0,1),  col="green", xlab="betaB", ylab="density across precincts, f(betaB)", main="Density of betaB")
vb <- as.vector(betabm)
for (i in 1:length(vb)){
	lines(c(vb[i], vb[i]), c(0,.25))
	}
}

#Density of betaw
betawd <- function(ei.object){
betaws <- ei.object$betaws
betawm <- apply(betaws,1,mean)
plot(density(betawm), xlim=c(0,1), col="green", xlab="betaW", ylab="density across precincts, f(betaW)", main="Density of betaW")
vw <- as.vector(betawm)
for (i in 1:length(vw)){
	lines(c(vw[i], vw[i]), c(0,.25))
	}
}


#XT plot
xt <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T Scatterplot", ylab="T", xlab="X", pch=20)
}

#XTc plot

xtc <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
circ <- .04
plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T Scatterplot with Population Density Circles", ylab="T", xlab="X", pch=20)
minn <- min(n)
maxn <- max(n)
for (i in 1:length(x)){
	radius = (n[i]-minn+1)/(1+maxn-minn)
	draw.circle(x[i], t[i], radius*circ)
	}
	}
#xtc(data$x,data$t,data$n,.04, "X and T Scatterplot")


#XTfit plot

xtfit <- function(ei.object){
betabs <- ei.object$betabs
betaws <- ei.object$betaws
low <- .1
up <- .9
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
circ <- .04
plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T Scatterplot with E(T|X) and 80% CIs", ylab="T", xlab="X", pch=20)
minn <- min(n)
maxn <- max(n)
for (i in 1:length(x)){
	radius = (n[i]-minn+1)/(1+maxn-minn)
	draw.circle(x[i], t[i], radius*circ)
	}
x <- seq(0,1,by=.01)
betabs <- as.vector(betabs)
betaws <- as.vector(betaws)
t <- matrix(ncol=length(x), nrow=length(betabs))
for(i in 1:length(x)){
t[,i] <- betabs*x[i] + betaws*(1-x[i])
	}
et <- apply(t,2,mean)
lines(x,et, col="yellow")
lwr <- apply(t,2,function (x) quantile(x, probs=c(low)))
upr <- apply(t,2,function (x) quantile(x, probs=c(up)))
lines(x, lwr, col="red")
lines(x, upr, col="red")
}

#xtfit(x,t,n,.04,"",ei.1$betabs, ei.1$betaws, .2,.8)


#XTfitg plot
xtfitg <- function(ei.object){
betabs <- ei.object$betabs
betaws <- ei.object$betaws
low <- .1
up <- .9
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
circ <- .04
plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T Scatterplot with E(T|X), 80% CIs, and Goodman", ylab="T", xlab="X", pch=20)
minn <- min(n)
maxn <- max(n)
for (i in 1:length(x)){
	radius = (n[i]-minn+1)/(1+maxn-minn)
	draw.circle(x[i], t[i], radius*circ)
	}
x <- seq(0,1,by=.01)
betabs <- as.vector(betabs)
betaws <- as.vector(betaws)
t <- matrix(ncol=length(x), nrow=length(betabs))
for(i in 1:length(x)){
t[,i] <- betabs*x[i] + betaws*(1-x[i])
	}
et <- apply(t,2,mean)
lines(x,et, col="yellow")
lwr <- apply(t,2,function (x) quantile(x, probs=c(low)))
upr <- apply(t,2,function (x) quantile(x, probs=c(up)))
lines(x, lwr, col="red")
lines(x, upr, col="red")
t <- ei.object$t
x <-ei.object$x
lm.fit <- lm(t ~ x)
abline(lm.fit, col="green")
}

#Goodman plot
goodman <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
plot(x, t, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="X and T Scatterplot with Goodman", ylab="T", xlab="X", pch=20)
lm.fit <- lm(t ~ x)
abline(lm.fit, col="red")
}


#Estsims plot

estsims <- function(ei.object){
betabs <- as.vector(ei.object$betabs)
betaws <- as.vector(ei.object$betaws)
colors = runif(length(betabs),26,51)
plot(betabs, betaws, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="Simulations of betaW and betaB", ylab="betaW simulations", xlab="betaB simulations", pch=20, col=colors, lty=2, cex=.25)
}

#boundXB


boundXb <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
truebb <- ei.object$truth[,1]
bounds <- bounds1(x, t, n)
plot(x, truebb, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="Aggregation Bias for betaB", ylab="True betab", xlab="X", pch=20)
for (i in 1:length(x)){
	lines(c(x[i], x[i]), c(bounds[,1][i], bounds[,2][i]))
	}
lm.xb <- lm(truebb ~ x)
abline(lm.xb, lty=2)
}


boundXw <- function(ei.object){
x <- ei.object$x
t <- ei.object$t
n <- ei.object$n
truebw <- ei.object$truth[,2]
bounds <- bounds1(x, t, n)
plot(x, truebw, xlim=c(0,1), ylim=c(0,1),xaxs="i",yaxs="i", main="Aggregation Bias for betaW", ylab="True betaw", xlab="X", pch=20)
for (i in 1:length(x)){
	lines(c(x[i], x[i]), c(bounds[,3][i], bounds[,4][i]))
	}
lm.xw <- lm(truebw ~ x)
abline(lm.xw, lty=2)
}


#truth
truthfn <- function(ei.object){
n <- ei.object$n
truebb <- ei.object$truth[,1]
truebw <- ei.object$truth[,2]
betabs <- ei.object$betabs
betaws <- ei.object$betaws
truthbb <- sum(truebb*n)/sum(n)
truthbw <- sum(truebw*n)/sum(n)
par(mfrow=c(1,2))
ag <- aggs(x,t,n,betabs,betaws)
plot(density(ag[,1]), xlim=c(0,1),ylim=c(0,max(density(ag[,1])$y)+1), yaxs="i",xaxs="i", main="Density of Bb Posterior & Truth", xlab="Bb",ylab="Density")
lines(c(truthbb, truthbb), c(0,.25*(max(density(ag[,1])$y)+1)), lwd=3)
plot(density(ag[,2]), xlim=c(0,1),ylim=c(0,max(density(ag[,2])$y)+1), yaxs="i", xaxs="i", main="Density of Bw Posterior & Truth", xlab="Bw",ylab="Density")
lines(c(truthbw, truthbw), c(0,.25*(max(density(ag[,2])$y)+1)), lwd=3)
}
