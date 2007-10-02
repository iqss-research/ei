##/*
##** betabs = einonp(t,x);
##**
##** INPUTS: t,x = (p X 1) turnout, race
##**
##** OUTPUTS: betabs = (p x _Esims) simulations of beta^b, 
##**                   computed nonparametrically
##*/
einonp <- function(t,x, evbase=parent.frame()){

  evloc <- getEnvVar(evbase, environment())###, vecvar=c("eimetar"))
  
  nobs <- rows(t);
  
  if(Eprt>=2)
    message("Nonparametric Density Estimation...")
 
  Enumtol <- EnumTol <- as.vector(EnumTol)
  xr <- recode(x,cbind(x<Enumtol,x>(1-Enumtol)),rbind(Enumtol, (1-Enumtol)));
  tr <- recode(t,cbind(t<Enumtol, t>(1-Enumtol)),rbind(Enumtol, (1-Enumtol)));
  x[x<Enumtol] <- EnumTol
  t[t<Enumtol] <- EnumTol
  x[x>(1-Enumtol)] <- 1-EnumTol
  t[t>(1-Enumtol)] <- 1-EnumTol
  if(any(x != xr) || any(t != tr))
    stop("Recode does not work")
  lst <- bounds1(t,x,1);
  
  bnds <- lst[[1]]
  tt <- lst[[2]]
  betab <- matrix(0, nrow=EnonEval,ncol=nobs);

  for (i in 1:nobs)
    betab[,i+0] <- seqase(bnds[i+0,1],bnds[i+0,2],EnonEval)
 
  x1 <- 1-x;
  betaw <- as.vector((t/x1))-as.vector((x/x1))*betab;
  
  if (Eprt>=2)
    message("Computing conditional density for each line...");
  
  pz <- nonbiv(t,x,betab,betaw, evbase);
  
  if(Eprt>=2)
    message("Drawing simulations...");
  
  ###/* compute interpolated trapazoidal areas */
  areas <- (betab[2,]-betab[1,])*na.omit(pz+0.5*abs(pz-lag(pz)));
  nr <- rows(areas)
  areas <- na.omit(areas)
  dc <- nr - rows(areas)
  areas <- cumsum(as.data.frame(areas/as.vector(colSums(areas))));
   
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
  message("betabs <- einonp()")
  print(betabs)
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
##*
##** If memory is available, use vec(px) and vec(py) and reshape pz upon output
##** to make this proc run faster.
##*/
nonbiv <- function(t,x,px,py, evbase=parent.frame()){

  eigraph.bvsmth <- get("eigraph.bvsmth",env=evbase)
  col <- cols(px);
  pz <- matrix(0, nrow=rows(px),ncol=col)

  c <- unitarea(t,x,evbase);###         @ scale factor to divide by @

  r <- rows(x);
 
  for (i in 1:col){
    d <- perpdist(t,x,px[,i+0],py[,i+0],evbase) ### @ perpendicular distance to line @

  ###  /* sheet-normal kernel */
    z <- d / eigraph.bvsmth;
    ln <- cols(z)
    mat <- matrix(0, nrow=rows(c), ncol=cols(z))
    for(n in 1:length(c))
      mat[n, ] <- exp(-0.5*(z*z))/c[n]
    csum <- colSums(mat)
    pz[,i+0]  <- 1 / sqrt(2*pi) * csum/r
  
  }
  

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
##*/
unitarea <- function(t,x, evbase=parent.frame()){

  EnonNumInt <- get("EnonNumInt", env=evbase)
  Enumtol <- as.vector(get("EnumTol", env=evbase))
  eigraph.bvsmth <- get("eigraph.bvsmth", env=evbase)
  xr <- recode(x,cbind(x<Enumtol,x>(1-Enumtol)),rbind(Enumtol, (1-Enumtol)))
  x[x < Enumtol] <- Enumtol
  x[x > (1-Enumtol)] <- (1-Enumtol)
  if(any(x != xr))
    stop("Recode does not work")
 
  lst <- bounds1(t,x,1);
 
  bnds <- lst[[1]]
  tt <- lst[[2]]
  lb <- lB <- na.omit(bnds[,1]);
  ub <- uB <- na.omit(bnds[,2]);
  lw <- lW <- na.omit(bnds[,3]);
  uw <- uW <- na.omit(bnds[,4]);
  x1 <- 1-x;
  xx1 <- x/x1;
  x1x <- x1/x;
  nobs <- rows(x);
  area <- matrix(0, nrow=nobs,ncol=1);
  var <- eigraph.bvsmth^2;

  for (i in 1:nobs){
    
    c <- maxr(cbind(1-uB[i+0],lW[i+0]),maxr(lB[i+0],1-uW[i+0]))
   
    c <- substute(c,c<0.00001,c*0+0.00001);
    ###c[ c<0.00001] <- c*0+0.00001
    ###if(any(cs != c))
    ###  stop("substute not working")

    d <- perpdist(t[i+0],x[i+0],c(1,0),c(0,1));

    mx <- max(c(length(c), length(d)))
    if(!is.na(mx))
      {
        if(length(c) < mx && length(c))
          c <- rep(c, mx/length(c))
        if(length(d) < mx && length(d))
          d <- rep(d, mx/length(d))
      }
  
    a <- sqrt(c^2-d^2);
    k <- a*cos(asin(a/c));
    mx <- max(c(length(a), length(k)))
     if(!is.na(mx))
      {
        if(length(a) < mx && length(a))
          a <- rep(a, mx/length(a))
        if(length(k) < mx && length(k))
          k <- rep(k, mx/length(k))
      }
     
    g <- sqrt(a^2-k^2);
  
  ###  /* coordinates of tomography line extended so that the end points
  ##  are perpendicular to the 1,0 0,1 coordinates of the unit square:
  ##  (uB+k[.,1])~(lW-g[.,1])~(lB-k[.,2])~(uW+g[.,2]);  */
    
  ###  /* points to evaluate on the tomography line */
   
    
    betabS <- seqase(lb[i+0]-k[,2],ub[i+0]+k[,1],EnonNumInt);
     
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
   ##  print(a)
   ## print(b)
   ## print(var)
   ## print(S)
    area[i+0] <- colMeans(cdfnorm(maxr(a,b),0,var)+((S-0.5)*2)*cdfnorm(minr(a,b)-S,0,var));
    
  }
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
##**/
perpdist <- function(t,x,px,py,evbase=parent.frame()){
 
  t <- matrix(t, ncol=length(t))
  x <- matrix(x, ncol=length(x))
  px <- as.matrix(px)
  py <- as.matrix(py)

  lst <- bounds1(t,x,1);
  bnds <- lst[[1]]
 
  tt <- lst[[2]]

  lB <- bnds[,1];
  uB <- bnds[,2];
  lW <- bnds[,3];
  uW <- bnds[,4];
 
  A <- sqrt((uW-lW)^2+(uB-lB)^2);
  a <- A[A <=0.0001] <- 0.0001
  aR <- AR <- recode(A,A<=0.0001,0.0001);
  if(any(AR != A))
    stop("Recode is not working properly")
  B <- sqrt((uW-t(py))^2+(lB-t(px))^2);
  
  b <- B <- substute(B,B<0.00001,B*0+0.00001);
 ### b <- B[B<0.00001] <- B*0+0.00001
 ### if(any(B != Bs))
 ###   stop("substute not working properly")
  c <- C <- sqrt((lW-t(py))^2+(uB-t(px))^2);
  mx <- max(c(length(a), length(b), length(c)))
    if(!is.na(mx))
      {
        if(length(a) < mx)
          a <- A <- rep(a, mx/length(a))
    
        if(length(b) < mx)
          b <- B <- rep(b, mx/length(b))
  
        if(length(c) < mx)
          c <- C <- rep(c, mx/length(c))
      }
  
  acarg <- (a^2+b^2-c^2)/(2 *a *b);
  tst1 <- (acarg>1);
  tst2 <- (acarg < -1);
  acarg <- substute(acarg,rbind(tst1,tst2),0*acarg+tst1-tst2);
  acarg <- na.omit(acarg)
  
  D <- B*sin(acos(acarg));
  if(dim(D)[[1]] != dim(x) || dim(D)[[2]] != dim(px))
    stop("dimension perpdist do not agree with dim(x) X dim(px)")
  return(D)
}
