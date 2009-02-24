##/*
##** {a1,b1}=listwis2(a,b);
##**
##** INPUTS:  a,b have same number of rows
##**
##** a1,b2, equal a,b except rows with any missing values are deleted.
##*/
listwis2 <- function(a,b){
  if(rows(a) != rows(b)) {
    message("listwis2: applies only to matrices same number of rows")
    return(NULL)
  }
  ###vector indicating if rows contain NA's 
  ind <- rowSums(is.na(cbind(a,b))) >=1
  a1 <- subset(a,!ind)
  b1 <- subset(b, !ind)
  lst <- c(list(a1), list(b1))
  names(lst) <- c("a", "b")
  return(lst)
}
###/*
###** {a1,b1,c1}=listwis3(a,b,c);
##**
##** INPUTS:  a,b,c have same number of rows
##**
##** a1,b1,c1 equal a,b,c except rows with any missing values are deleted.
###*/
listwis3 <- function(a,b,c){
 if(rows(a) != rows(b) || rows(a) != rows(c) || rows(b) != rows(c)) {
    message("listwis3: applies only to matrices same number of rows")
    return(NULL)
  }
  ###vector indicating if rows contain NA's 
  ind <- rowSums(is.na(cbind(a,b,c))) >=1
  a1 <- subset(a,!ind)
  b1 <- subset(b, !ind)
  c1 <- subset(c,!ind)
  lst <- c(list(a1), list(b1), list(c1))
  names(lst) <- c("a", "b","c")
  return(lst)
}
##/*
##**  (C) Copyright 1999 Gary King
##**  All Rights Reserved.
##**  http://GKing.Harvard.Edu, King@Harvard.Edu
##**  Department of Government, Harvard University
##**
##** triple scatter plot
##**
##** Usage:  call triple(x,y,z);
##**
##**         x is plotted horizontally,
##**         y is plotted vertically,
##**         z is the size of the circle to be plotted at the x,y coordinate.
##**           it should be linearly scaled to be from about 0.5 to about 10.
##**
##**         xtics and ytics should be set before running this proc.
##*/
##external matrix _pxscale;
##external matrix _pyscale;
##declare matrix _psym != 0;

triple <- function(x,y,z, xlabel="", ylabel="",
                   title="",pxscale=0,pyscale=0,pzscale=1,evbase=NULL){
 
  z <- as.matrix(z)
 
   z <- as.matrix(z)*pzscale
  if( rows(z)==1){  ###plots circles
    os <- matrix(1, nrow=rows(x),ncol=1)
    z <- os%*%z
  }
    lst <- listwis3(x,y,z)
    x <- lst[[1]]
    y <- lst[[2]]
    z <- lst[[3]]
    os <- matrix(1, nrow=rows(x),ncol=1)
    minx <- miny <- 0
    maxx <- maxy <- 1
    if(pxscale==0){
        minx <- floor(minc(x)*10)/10
        maxx <- ceiling(maxc(x)*10)/10
      }
 
    if(pyscale==0){
        miny <- floor(minc(y)*10)/10
        maxy <- ceiling(maxc(y)*10)/10;

      }
  plot(x,y,cex = z, xlim=c(minx,maxx), ylim=c(miny,maxy),
       xlab=xlabel, ylab=ylabel,main=title)
}
