#setwd("/Users/mollyroberts/Winter2010/GaryResearch")
#data <- read.table("sample.csv", header=T)
#names(data) <- tolower(names(data))
bounds1 <- function(x,t,n){
	homindx <- NULL; tx <- NULL; tomx <- NULL; LbetaB <- NULL; UbetaB <- NULL; LbetaW <- NULL; UbetaW <- NULL;
	omx = 1-x
	Nb = x*n
	Nw =omx*n
	p = length(x)
	homoindx <- ifelse(x==0, 1, 0)
	homoindx <- ifelse(x==1, 2, homoindx)

#Heterogenous precincts
tx <- as.matrix(t/x)
tomx = as.matrix(t/omx)
tomxx <- as.matrix(tx-(omx/x))
txx <- as.matrix(tomx-x/(1-x))
LbetaB <- apply(tomxx, 1, function (x) max(0,x))
UbetaB <- apply(tx, 1, function (x) min(x,1))
LbetaW <- apply(txx, 1, function (x) max(0, x))
UbetaW <- apply(tomx, 1, function (x) min(x, 1))

bl <- homoindx==2
#Homogenously black	
		LbetaB[bl] = t[bl]
		UbetaB[bl] = t[bl]
		LbetaW[bl] = NA
		UbetaW[bl] = NA
	

#Homogenously white
wh <- homoindx==1
		LbetaB[wh] = NA
		UbetaB[wh] = NA
		LbetaW[wh] = t[wh]
		UbetaW[wh] = t[wh]
		
return(cbind(LbetaB, UbetaB, LbetaW, UbetaW))
}
	
