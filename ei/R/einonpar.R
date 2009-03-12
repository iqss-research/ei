##/*
##** betabs = einonp(t,x);
##**
##** INPUTS: t,x = (p X 1) turnout, race
##**
##** OUTPUTS: betabs = (p x _Esims) simulations of beta^b, 
##**                   computed nonparametrically
##
## FUNCS: getEnvVar,bounds1, nonbiv 
##*/
einonp <- function(t,x, evbase=parent.frame()){
  t <- matrix(t,ncol=1)
  x <- matrix(x,ncol=1)
  evei <- evloc <- getEnvVar(evbase, environment())###, vecvar=c("eimetar"))
  if(exists("Eprt")) Eprt <- get("Eprt", env=evei)
  if(exists("EnonEval")) EnonEval <- get("EnonEval", env=evei)
  if(exists("Esims")) Esims <- get("Esims", env=evei)
  nobs <- rows(t);
  
  if(Eprt>=2)
    message("Nonparametric Density Estimation...")
 
  Enumtol <- EnumTol <- as.vector(EnumTol)
  x <- recode(x,cbind(x<Enumtol,x>(1-Enumtol)),rbind(Enumtol, (1-Enumtol)));
  t <- recode(t,cbind(t<Enumtol, t>(1-Enumtol)),rbind(Enumtol, (1-Enumtol)));
  lst <- bounds1(t,x,1,Enumtol);
  
  bnds <- lst[[1]]
  tt <- lst[[2]]
  betab <- matrix(0, nrow=EnonEval,ncol=nobs);

  for (i in 1:nobs)
    betab[,i+0] <- seqase(bnds[i+0,1],bnds[i+0,2],EnonEval)
 
  x1 <- 1-x;
 
  betaw1 <- betaw <- betab
  for(n in 1:cols(betab))
    betaw1[,n] <- (x/x1)[n]* betab[,n]
  for(n in 1:cols(betaw1))
    betaw[,n] <- (t/x1)[n] - betaw1[,n]
 
  
  if (Eprt>=2)
    message("Computing conditional density for each line...");
  
  pz <- nonbiv(t,x,betab,betaw, evbase);
  
  if(Eprt>=2)
    message("Drawing simulations...");
  
  ###/* compute interpolated trapazoidal areas */
  areas1 <- (betab[2,]-betab[1,])
 
  areas2 <- na.omit(pz+0.5*abs(pz-lag(pz)))
  
  areas <- areas1 %dot*% areas2
 
 
  colS <- colSums(areas)
  areas <- areas %dot/% colS
  
  areas <- cumsum(as.data.frame(areas));
  areas <- as.matrix(areas)
  
  betabs <- matrix(0,nrow=nobs,ncol=Esims);
  for (i in 1:nobs){
  
 ###   /* count within categories defined by areas of interpolated trapazoids */
    cnts <- counts(matrix(runif(Esims), nrow=Esims, ncol=1),areas[,i+0]);
   
 ###   /* use inverse CDF to draw from trapazoidal dist  within categories */
    fa <- trimr(lag(pz[,i+0]),1,0);
    fb <- trimr(pz[,i+0],1,0);
    a <- trimr(lag(betab[,i+0]),1,0);
    b <- trimr(betab[,i+0],1,0);
    D <- fa*b-fb*a;
    E <- fb-fa
   
    E[abs(E)<0.000001] <- 0.000001; ###@ fix for uniform distribution @
    C <- 1/( (b-a)*D+0.5*(b^2-a^2)*E );
    CE <- C*E;
    CD <- C*D;
    
    res <- 0;
    for (j in 1:(EnonEval-1)){
      if (!is.na(cnts[j+0]) &&  cnts[j+0]>0 ){
        res <- rbind(res, (-CD[j+0]+sqrt((C[j+0]^2)*(D[j+0]^2)+ 2*CE[j+0]*(CD[j+0]*a[j+0]+0.5*CE[j+0]*(a[j+0]^2)+as.matrix(runif(cnts[j+0])))))/CE[j])
      }
    }

    res <- as.vector(trimr(res,1,0))
    betabs[i+0,] <- res[order(runif(Esims))];

  }
  if(dim(betabs)[[1]] != length(x) ||  dim(betabs)[[2]] != Esims)
    stop("Dimension betabs do not agree with spec")
 
  return(betabs);
}
##/*
##** pz = nonbiv(t,x,px,py);
##**
##** INPUTS: t, x, = (px1) turnout, race
##**         px,py = (MxQ) coordinates to evaluate nonparametric biv density at
##**
##** OUTPUT: pz = (MxQ) height of the nonparametric bivariate density at px,py,
##**               which is f(px,py), using a sheet-normal kernel.
##
###FUNCS: unitarea, perpdist
##*
##** If memory is available, use vec(px) and vec(py) and reshape pz upon output
##** to make this proc run faster.
##*/
nonbiv <- function(t,x,px,py, evbase=parent.frame(),eigraph.bvsmth=NULL,Enumtol=NULL){
  if(length(evbase)){
    eigraph.bvsmth <- get("eigraph.bvsmth",env=evbase)
    Enumtol <- as.vector(get("EnumTol", env=evbase))
  }
  t <- as.matrix(t)
  x <- as.matrix(x)
  px <- as.matrix(px)
  py <- as.matrix(py)
  col <- cols(px);
  pz <- matrix(0, nrow=rows(px),ncol=col)
  c0 <- unitarea(t,x,evbase);###         @ scale factor to divide by @
 
  r <- rows(x);
  pz <- sapply(as.list(1:col), function(i,t,x,px,py,Enumtol, eigraph.bvsmth,c0,r){
  
   
    d <- perpdist(t,x,px[,i+0],py[,i+0],Enumtol) ### @ perpendicular distance to line @
  ###  /* sheet-normal kernel */
    z <- d / eigraph.bvsmth;
    ln <- cols(z)
    c <- matrix(c0, nrow=length(c0), ncol=ncol(z))
    mat <- exp(-0.5*(z*z))/c
    csum <- colSums(mat)
    return(pz  <- csum/r/sqrt(2*pi))},t,x,px,py,Enumtol, eigraph.bvsmth,c0,r)
  
  pz <- unlist(pz)
  
  pz <- pz/eigraph.bvsmth^2;
  if(any(dim(pz)!= dim(px)) ||any(dim(pz)!= dim(py)))
    stop("Wrong dimensions for pz <- nonbiv()")
  return(pz)
}
##/*  area = unitarea(t,x);
##**  
##**  t,x = (px1) turnout, race
##**
##**  area = (px1) area within unit square of truncated nonparametric 
##**         normal-sheet kernel for each tomography line
##
## FUNCS:bounds1, maxr, perpdist, seqas 
##*/
unitarea <- function(t,x, evbase=parent.frame(),EnonNumInt=NULL, Enumtol=NULL,eigraph.bvsmth=NULL ){
  if(length(evbase)){
    EnonNumInt <- get("EnonNumInt", env=evbase)
    Enumtol <- as.vector(get("EnumTol", env=evbase))
    eigraph.bvsmth <- get("eigraph.bvsmth", env=evbase)
  }
  x <- recode(x,cbind(x<Enumtol,x>(1-Enumtol)),rbind(Enumtol, (1-Enumtol)))
  lst <- bounds1(t,x,1,Enumtol);
 
  bnds <- lst[[1]]
  tt <- lst[[2]]
  lb <- lB <- bnds[,1];
  ub <- uB <- bnds[,2];
  lw <- lW <- bnds[,3];
  uw <- uW <- bnds[,4];
  x1 <- 1-x;
  xx1 <- x/x1;
  x1x <- x1/x;
  nobs <- rows(x);
  area <- matrix(0, nrow=nobs,ncol=1);
  var <- eigraph.bvsmth^2;
 
    
   area <- lapply(as.list(1:nobs), function(i,uB,lW,lB,uW,t,x,x1,x1x,xx1,var,Enumtol,EnonNumInt ){
     lb <- lB
     ub <- uB
    c <- cbind(maxr(1-uB[i+0],lW[i+0]),maxr(lB[i+0],1-uW[i+0]))
   
    c <- substute(c,c<0.00001,c*0+0.00001);
   
    d <- perpdist(t[i+0],x[i+0],c(1,0),c(0,1), Enumtol);
 
    a <- sqrt(c^2-d^2);
    k <- a*cos(asin(a/c));

    g <- sqrt(a^2-k^2);
    
  ###  /* coordinates of tomography line extended so that the end points
  ##  are perpendicular to the 1,0 0,1 coordinates of the unit square:
  ##  (uB+k[.,1])~(lW-g[.,1])~(lB-k[.,2])~(uW+g[.,2]);  */
    
  ###  /* points to evaluate on the tomography line */
    betabS <- seqas(lb[i+0]-k[,2],ub[i+0]+k[,1],EnonNumInt);
    betawS <- (t[i+0]/x1[i+0])-(x[i+0]/x1[i+0])*betabS;
    z <- matrix(0,nrow=EnonNumInt,ncol=1);
    o <- matrix(1, nrow=EnonNumInt,ncol=1);
    
  ###  /* lengths of perpendicular lines within the unit square */
    minbetaw <- maxr(z,betawS-x1x[i+0]*betabS);
    maxbetaw <- minr(o,betawS+x1x[i+0]*(1-betabS));
    minbetab <- maxr(z,betabS-xx1[i+0]*betawS);
    maxbetab <- minr(o,betabS+xx1[i+0]*(1-betawS));
    a <- sqrt((maxbetaw-betawS)^2+(maxbetab-betabS)^2);
    b <- sqrt((betabS-minbetab)^2+(betawS-minbetaw)^2);
    S <- (betabS>=minbetab)&(betabS<=maxbetab);
    
   ### /* area within square for each perpendicular line */
   
    return(area <- colMeans(cdfnorm(maxr(a,b),0,var)+((S-0.5)*2)*cdfnorm(minr(a,b)-S,0,var)))},uB,lW,lB,uW,t,x,x1,x1x,xx1,var,Enumtol,EnonNumInt )
     
    
  area <- unlist(area)
   if(any(dim(area)!= dim(t)))
      stop("area dimension do not agrre with x and t")
   
  return(area);
}


##/* dist = perpdist(t,x,px,py);
##**
##** INPUTS:  t,x   = (px1) turnout, race
##**          px,py = (Mx1) coordinates of points to compute distance
##**
##** OUTPUT: dist = (PxM) perpendicular distance from each point in px,py
##**                 to the tomography line betaW=(t./(1-x))-(x./(1-x))*betaB
##
## FUNCS: bounds1
##**/
perpdist <- function(t,x,px,py,tol){
 
  t <- matrix(t, ncol=1)
  x <- matrix(x,ncol=1)
  px <- matrix(px,ncol=1)
  py <- matrix(py,ncol=1)

  lst <- bounds1(t,x,1,tol);
  bnds <- lst[[1]]
 
  tt <- lst[[2]]

  lB <- as.matrix(bnds[,1]);
  uB <- as.matrix(bnds[,2]);
  lW <- as.matrix(bnds[,3]);
  uW <- as.matrix(bnds[,4]);
 
  A <- sqrt((uW-lW)^2+(uB-lB)^2);
 
  a <- A <- recode(A,A<=0.0001,0.0001);

  B <- sqrt((uW %-% py)^2+(lB %-% px)^2);
  b <- B <- substute(B,B<0.00001,B*0+0.00001);
  c <- C <- sqrt((lW %-% py)^2 + (uB %-% px)^2); 
  a <- A <- matrix(A, nrow=rows(A), ncol=cols(B))

  acarg <- (a^2+b^2-c^2)/(2 *a *b);
  tst1 <- (acarg>1);
  tst2 <- (acarg < -1);
  acarg <- substute(acarg,(tst1 |tst2),0*acarg+tst1-tst2);
###  acarg <- na.omit(acarg)
 
   
  D <- B*sin(acos(acarg));
  if(dim(D)[[1]] != dim(x) || dim(D)[[2]] != dim(px))
    stop("dimension perpdist do not agree with dim(x) X dim(px)")
  
  return(D)
}
 
