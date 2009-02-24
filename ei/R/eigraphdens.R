##/* KERNEL DENSITY ESTIMATION
##**
##**   Format:    {px,py}=dens(y);
##**
##**   input:     y = Nx1 vector
##**
##**   output:   plot of y (horizontally) by its est'd density f(y) (vertically).
##**             Whiskers are also plotted along horizontal axis for each datum
##**
##**              px = x coordinates for the plot (_pts x 1 vector)
##**              py = y coordinates for the plot (_pts x 1 vector)
##**
##**    Global values are as follows, and may be changed:
##**   _smth=0;         smoothing parameter>0. the larger, the more smoothing.
##**                    =0 for automatic calculation of density.
##**                    =-1 for automatic calculation for description.
##**
##**   _strt=0;         x axis minimum. (if _strt=_endd, set at min and max).
##**   _endd=0;         x axis maximum.
##**   _pts=100;        number of points to plot
##**   _whiskr=-1;      -1 = draw whiskers at 1/10th size;
##**                    >0 = draw whiskers at size _whiskr;
##**                     0 = do not draw whiskers
##**                    Whiskers are always drawn at y==0, so you may have to
##**                    use ytics to scale the y-axis.
##**   _jitter=0;       0=nothing extra. #=add # jitter to each
##**                      whisker if _whiskr/=0;.
##**
##**   let _kern=E;     kernel type. N=normal, E=Epanechnikov,
##**                                 B=biweight, T=triangular, R=rectangular
##**                                 TN=Doubly Truncated Normal
##**  for option _kern=="TN",
##**      _Tleft=_strt;   truncate distribution on left at _Tleft
##**      _Tright=_endd;  truncate distribution on right at _Tright
##**
##**   _output = 1;     1 = print density plot (default); 0 = do not print.
##*/
##@ Kernels @
##@ arguments of kernels: z=input value, m=mean, h=smoothing parameter _smth @
kerneln <- kernelN <- function(z,m,h){
###  pi <- 3.1415927
  z <- as.matrix(z) ###scalar
  if(length(z)<= 1) z <- z[1,1]
  
  m <- as.matrix(m) ###matrix Esims rows, 1column
  if(length(m) <= 1) m <- m[1,1] ### scalar
  
  h <- as.matrix(h) ###scalar
  if(length(h)<=1) h <- h[1,1]
  
  z <- (z-m)%dot/%h ###      @ NORMAL kernel @
  res <- (1/sqrt(2*pi))*exp(-(1/2)*(z^2))
  return(res)
}
kernelTN <- kerneltn <- kernelTn <- function(z,m,h,Tleft,Tright){
  z <- as.matrix(z) ###scalar
  if(length(z)<= 1) z <- z[1,1]
  
  m <- as.matrix(m) ###matrix Esims rows, 1column
  if(length(m) <= 1) m <- m[1,1] ### scalar
  
  h <- as.matrix(h) ###scalar
  if(length(h)<=1) h <- h[1,1]
  
  Tleft <- as.matrix(Tleft) ###scalar
  if(length(Tleft)<=1) Tleft <- Tleft[1,1]
  
  Tright <- as.matrix(Tright) ###scalar
  if(length(Tright)<=1) Tright <- Tright[1,1]
  
  zz <- (z-m)%dot/%h    ###      @ TRUNCATED NORMAL kernel @
    tl <- (Tleft-m)%dot/%h
    tr <- (Tright-m)%dot/%h
    zvec <- as.vector(z)
    t <- as.numeric(((zvec >Tleft) & (zvec <Tright)))
     
    res <- (t%dot*%dnorm(zz))%dot/%(1-pnorm(tl)-(1-pnorm(tr)))
    return(res)
  }
kernele <- kernelE <- function(z,m,h){
   z <- as.matrix(z) ###scalar
  if(length(z)<= 1) z <- z[1,1]
  
  m <- as.matrix(m) ###matrix Esims rows, 1column
  if(length(m) <= 1) m <- m[1,1] ### scalar
  
  h <- as.matrix(h) ###scalar
  if(length(h)<=1) h <- h[1,1]
  
  z <- (z-m)%dot/%h ###  @ EPANECHNIKOV kernel @
  t <- (abs(z) <sqrt(5))
   a <- code(t,rbind(sqrt(5),1))
  res <- t%dot*% ((3/4)*(1-(1/5)*(z^2))%dot/%a)
  return(res)
}
kernelb <- kernelB <- function(z,m,h){
  z <- as.matrix(z) ###scalar
  if(length(z)<= 1) z <- z[1,1]
  
  m <- as.matrix(m) ###matrix Esims rows, 1column
  if(length(m) <= 1) m <- m[1,1] ### scalar
  
  h <- as.matrix(h) ###scalar
  if(length(h)<=1) h <- h[1,1]
  
  z <- (z-m)%dot/%h ###  @ BIWEIGHT kernel @
  t <- as.numeric(abs(z)<1)
  res <- t%dot*%((15/16)*((1-(z^2))^2))
  return(res)
}
kernelt <- kernelT <- function(z,m,h){
  z <- as.matrix(z) ###scalar
  if(length(z)<= 1) z <- z[1,1]
  
  m <- as.matrix(m) ###matrix Esims rows, 1column
  if(length(m) <= 1) m <- m[1,1] ### scalar
  
  h <- as.matrix(h) ###scalar
  if(length(h)<=1) h <- h[1,1]
  
  z <- (z-m)%dot/%h ###  @ TRIANGULAR kernel @
  t <- as.numeric(abs(z)<1)
  res <- t%dot*%(1-abs(z))
  return(res)
}
kernelr <- kernelR <- function(z,m,h){
  z <- as.matrix(z) ###scalar
  if(length(z)<= 1) z <- z[1,1]
  
  m <- as.matrix(m) ###matrix Esims rows, 1column
  if(length(m) <= 1) m <- m[1,1] ### scalar
  
  h <- as.matrix(h) ###scalar
  if(length(h)<=1) h <- h[1,1]
  
  z <- (z-m)%dot/%h ###  @ RECTANGULAR kernel @
  t <- as.numeric(abs(z)<1)
  res <- t*0.5
  return(res)
}

###@ creates two nx1 vectors to plot @
density <- function(y,strt,endd,pts,h,kerna,Tleft,Tright){
   
    kerna <- as.vector(toupper(kerna))
    px <- as.matrix(seqas(strt,endd,pts))   
    py <- px
    
  ###  format /rdn 4,0;
    message("Kernel Density Estimation")
    message("Calculating ", pts," Points: ")
    for (i in 1:pts){
    
        if(identical(kerna,"N"))
          t <- (kernelN(py[i+0,1],y,h))%dot/%h
        else if( identical(kerna,"E")){
         
          t <- (kernelE(py[i+0,1],y,h))%dot/%h
        }else if( identical(kerna,"B"))
          t <- (kernelB(py[i+0,1],y,h))%dot/%h
        else if(identical(kerna,"T"))
          t <- (kernelT(py[i+0,1],y,h))%dot/%h
        else if(identical(kerna,"R"))
          t <- (kernelR(py[i+0,1],y,h))%dot/%h
        else if(identical(kerna,"TN"))
          t <- (kernelTN(py[i+0,1],y,h,Tleft,Tright))%dot/%h
        else
          stop( "kernel specified incorrectly")
       
            
        py[i+0,1] <- colSums(as.matrix( t ))/rows(y);
        
        if( (i+0)==10*floor(i/10))
          message("Index is ",i+0)
      }
    lst <- c(list(px),list(py))
    names(lst) <- c("px", "py")
    return(lst)
  }

###@ Kernel Density Estimate; sets defaults, calls density, then xy @
dens <- function(y,evbase=NULL){
    
    y <- as.matrix(y)
    if( cols(y)!=1)
      stop("Argument must be a column vector")
    strt <- get("strt", env=evbase)
    endd <- get("endd",env=evbase)
    if(strt>endd)
      stop("error: _strt>_endd")
    else if( strt==endd){
      strt <- minc(y)
      endd <- maxc(y)
    }
    pts <- floor(get("pts",env=evbase))
    if (pts<=2)
      stop("pts must be greater than 2.  Try pts=100;")
    smth <- get("smth", env=evbase)
    if ((smth<0)& (smth!=-1))
      stop("smth must be -1 or > than 0")
    else if( smth==0){
      py <- sortc(y,1)
      std <- minc(cbind((py[floor(3*rows(py)/4)]-py[floor(rows(py)/4)])/1.34,stdc(py)))
      smth <- 0.9*std*(rows(py)^(-0.2))
    }else if (smth==-1){
      py <- sortc(y,1)
      std <- minc(cbind((py[floor(3*rows(py)/4)]-py[floor(rows(py)/4)])/1.34,stdc(py)))
      smth <- 0.9*std*(rows(py)^(-0.2))/2
    }
    kern <-get("kern", env=evbase)
    
    Tleft <- get("Tleft", env=evbase)
    Tright <- get("Tright", env=evbase)
   
    if (identical(kern,"TN")){
      bool <- as.vector(Tleft==Tright)
      Tleft <- ifelse(bool, strt,Tleft)
      Tright <- ifelse(bool, endd, Tright)
   
    }
    jitter <- get("jitter", env=evbase)
    if(length(jitter)<=1) jitter <- jitter[1,1]
    lst <- density(y,strt,endd,pts,smth,kern,Tleft,Tright);
    px <- lst[["px"]]
    py <- lst[["py"]]
    
    std <- rndu(rows(y),1)*jitter-(jitter/2)
    y <- y+std
    whiskr <- get("whiskr", env=evbase)
    output <- get("output", env=evbase)
    
    if(output==1){
    if(whiskr==-1){
       os <- matrix(1,nrow=rows(y),ncol=1)
##	_pline=os~ (os*6)~ y~ (os*0)~ y~ (os*(maxc(py)/15))~ os~ (os*15)~(os*0);
       
    }else if(whiskr>0){
        os=matrix(1, nrow=rows(y),ncol=1)
  ###      _pline=os~ (os*6)~ y~ (os*0)~ y~ (os*_whiskr)~ os~ (os*15)~(os*0);
        
      }
        plot(px,py,type="l")
    ###    format/rd 5,4;
        if(rows(smth)==1)
          message("? Smoothing Parameter: smth=",smth)
  }
    lst <- c(px=list(px),py=list(py))
    return(lst)
  }
