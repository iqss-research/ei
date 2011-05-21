#Plot to help visualize multiple dimensions.

eiRxCplot <- function(ei.object, random =FALSE, black=TRUE, hispanic=TRUE,white=TRUE,informative=FALSE, threshold, title,xaxis, yaxis, percent, data, estimates=TRUE, legendpos = c(.69,1)){
  #ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  bounds <- bounds1(x,t,n)
  bbounds <- cbind(bounds[,1], bounds[,2])
  wbounds <- cbind(bounds[,4], bounds[,3])
  bounds3 <- na.omit(bounds)
  unan <- sum(bounds3[,1]==bounds3[,2] & bounds3[,1]==bounds3[,3] & bounds3[,1]==bounds3[,4] & bounds3[,1]==1)
  n <- dim(bounds)[1]
 
  plot(c(100,200), xlim=c(0,1), ylim=c(0,1),
       col="white", xaxs="i",
       yaxs="i", main=title, xlab=xaxis, ylab=yaxis)
ok2 <- NULL
 rand <- rep(TRUE, length(x))
  if (random ==TRUE) {
       rand <- sample(c(TRUE, FALSE),length(x),replace=T, prob=c(percent,1-percent))
}
  if (black == TRUE) {
    bl <- data$black.ei >threshold}
  if (hispanic==TRUE){
 his <- data$hisp.ei > threshold}
  if (white==TRUE){
 whit <- data$whit.ei > threshold
 }
 none <- !bl & !his & !whit
#print(sum(ok2))
  if (informative==TRUE) {
     ok2 <- bounds[,2]-bounds[,1] <.15 | bounds[,4]-bounds[,3] < .15
}
  for(i in 1:n){
    if (bl[i] == TRUE & rand[i]==TRUE)  lines(bbounds[i,], wbounds[i,], col="red", lwd=1)
    if (whit[i] == TRUE & rand[i]==TRUE)  lines(bbounds[i,], wbounds[i,], col="green", lwd=1)
    if (his[i] == TRUE & rand[i]==TRUE)  lines(bbounds[i,], wbounds[i,], col="blue", lwd=1)
  if (none[i] == TRUE & rand[i]==TRUE)  lines(bbounds[i,], wbounds[i,], col="black", lwd=1)
  }
if (estimates==TRUE){
  ok <- !is.na(ei.object$betab)&!is.na(ei.object$betaw)
  betabs <- ei.object$betabs[ok,]
  betaws <- ei.object$betaws[ok,]
  betabcd <- apply(betabs,1,function(x) quantile(x, probs=c(.1,.9)))
  betabm <- apply(betabs,1, mean)
  betawcd <- apply(betaws,1,function (x) quantile(x,probs=c(.1,.9)))
  betawm <- apply(betaws,1, mean)
  #n <- dim(betabcd)[2]
  for(i in 1:sum(ok)){
    if (rand[i]==TRUE) {lines(betabcd[,i], sort(betawcd[,i],decreasing=T), col="yellow",lwd=1.5)
    	}
    if(random==FALSE){
    points(betabm, betawm, col="yellow", cex=1, pch=16)
    }
  }
}
  #if (random==TRUE) text(.85,.97,paste("Unanimous",unan))
  legend(legendpos[1], legendpos[2],1,c("White", "Black", "Hispanic", "Mixed", "CI of Estimates"), col=c("green", "red", "blue", "black", "yellow"), lwd=1)
}

