##/*
####**  This archive is part of the program EI
####**  (C) Copyright 1995-2001 Gary King
##**  All Rights Reserved.
##*/
##/*
##** call eigraph(dbuf,"name");
##**
##** graphs results taken from from EI output data buffer in format of "name",
##** with options listed below.  eiread() is used to extract information from 
##** dbuf.  Some options require extra information in dbuf (such as "truth")
##**
##** When used with ei2 data buffers, eigraph uses the mean posterior estimate
##**
##** INPUTS:
##** dbuf     a data buffer, possibly created by ei()
##** name     a string with the name of the graph to produce
##**          chosen from this list below
##**
##** OUTPUT:##** numeric code indicating the success or failure of the call
##**          0  call succeeded
##**         -1  failed because needed truth and truth was not present
##**         -2  failed because inputs were bad or because no such graph
##**  
##** OPTIONS for name (*=combined graphs)
##** tomogD    tomography plot with data only
##** tomog     tomography plot with ML contours
##** tomogP    tomography plot with mean posterior contours 
##** tomogE    tomography plot with estimated betaB,betaW points
##** tomogT    tomography plot with true betaB,betaW points
##** tomogCI   tomography plot with 80% confidence intervals
##** tomogCI95 tomography plot with 95% confidence intervals
##** tomogS   *tomog,tomogp,tomogCI,Tbivar(or estsims if truth isn't available)
##** nonpar    tomogD & nonparametric density est with contours and surface plot
##** xt        basic x by t graph
##** xtC       basic x by t graph with circles sized proportional to N
##**           or some other variable defined in _eigraph_circ
##** Xgraph    an X-graph with data plotted
##** XgraphC   an X-graph with data plotted with size prop'l to N
##** goodman   x by t plot with goodman's regression line plotted
##** xtfit     x by t plot with E(T|X) and cond'l 80% confidence intervals 
##** xtfitg    xtfit with goodman's regression line superimposed
##** fit      *xtfit and tomogP
##** fitT     *xtfit and tomogT
##** postB     density estimate of (weighted) district aggregate B^b
##** postW     density estimate of (weighted) district aggregate B^w
##** post     *postB and postW
##** betaB     density estimate of est'd betaB's with whiskers at point ests
##** betaW     density estimate of est'd betaW's with whiskers at point ests
##** beta     *betaB & betaW
##** results  *postB, postW, betaB, betaW
##** movie     posterior of betaW_i and betaW_i for each i; hit key to cont.
##** movieD    highlight individual lines on a tomography plot; hit key to cont. 
##** prectB    plot of estimated betaB_i by true betaB_i
##** prectW    plot of estimated betaW_i by true betaW_i
##** truth    *post, prectB, prectW (cf truth to ests at dist and pcnt level)
##** lines     xt plot with one ESTIMATED line per precinct
##** Tlines    xt plot with one TRUE line per precinct
##** bivar     ESTimated betaB by betaW
##** Tbivar    TRUE betaB by betaW
##** betabw   *lines,Tlines,bivar,Tbivar (Tlines,Tbivar if truth is available)
##** profile   profile posterior plot of elements of phi
##** profileR  profile plot of elements of phi for R function (based on cdfbvn)
##** ptileB    true percentile at which betaB falls by est'd betaB
##** ptileW    true percentile at which betaW falls by est'd betaW
##** ptile    *ptileB & ptileW
##** simsB     sims of betaB by true betaB
##** simsW     sims of betaW by true betaW
##** sims     *simsB & simsW
##** estsims   simulated betaB's by simulated betaW's
##** biasB     X by EST'd betaB
##** biasW     X by EST'd betaW
##** TbiasB    X by TRUE'd betaB
##** TbiasW    X by TRUE'd betaW
##** bias     *biasB, biasW, TbiasB, TbiasW
##** boundXB   X by bounds on betaB
##** boundXW   X by bounds on betaW
##** boundX   *boundXW and boundXW
##**
### EV: I have added
###     betaXN thre dimensional plots fore betab and betaw as functions of X and N (or tvap)
###     betaTN thre dimensional plots fore betab and betaw as functions of T and N (or tvap)
#######
##** GLOBALS:
##** _eigraphC  = multiply by circle size to change sizes (>0)
##** _eigraph_Xlo = 0 low end of X graphs
##** _eigraph_Xhi = 1 high end of X graphs
##** _eigraph_Tlo = 0 low end of T graphs
##** _eigraph_Thi = 1 high end of T graphs
##** _eigraph_x  = "X"  xlabel for xt plots
##** _eigraph_t  = "T"  ylabel for xt plots
##** _eigraph_bb = "betaB" xlabel for bbXbw plots
##** _eigraph_bw = "betaW" ylabel for bbXbw plots
##** _eigraph_bblo = 0 low end for betaB plots
##** _eigraph_bbhi = 1 high end for betaB plots
##** _eigraph_bwlo = 0 low end for betaW plots
##** _eigraph_bwhi = 1 high end for betaW plots
##** _eigraph_loess = 1 show simulated and fitted loess for xtfit; 
##**                  0 fit only (default)
##** _eigraph_thick = add to line thickness parameter (default=1) 
##** _eigraph_bvsmth = bivariate density smoothing parameter (default=0.08)
##**             (this is stored in dbuf as 'bvsmth', but that value is not used
##**             to draw the contours under nonpar.)
##** _eigraph_eval = number of points to evaluate each side of nonparametric
##**               contour and surface plots (default=31)
##** _eigraph_dbuf = 0 don't do anything (default); 1 save inputs to contour
##** plot routine as elements in a data buffer called _eigraph_dbuf
##**               (px, py, pz are the elements)
##** _eigraph_smpl = 1.  Set this to (0,1] to randomly select this fraction of
##**               observations to use in tomography plots.  This is useful if 
##**               p is so large that it is difficult to see patterns. 
##** _tomogClr =  { 12, 9, 10, 11, 13, 5 } colors for each contour drawn 
##** _tomogPct = {.5, .95}    percentage values at which to draw contours
##**           rows(_tomogClr) must be >= rows(_tomogPct)
##** _EIMetaR = If dbuf is a meta-data buffer (_ei2_mta output from ei2),
##**          this global denotes which of the imputed data buffers stored
##**          in dbuf should be accessed when running this procedure (default=1).
##*/


###  require("lattice")
eigraph <- function(dbuf, str,psiu=NA,...){
 
  nm <- names(dbuf)
  nm <- sapply(nm,tolower)
  names(dbuf) <- nm
  if("evbase" %in% names(dbuf))
    evbase <-ev <- dbuf[["evbase"]]
  else if(exists("evbase"))
    ev <- evbase
  else
    ev <- evbase <- eiset(...)
 
 #### assigning the '..." to this environmnet
  drvdot <- match.call(expand.dots=TRUE)
  drv  <-  match.call(expand.dots=FALSE)
  entries <- expanddots(drvdot,drv,evbase)

  param  <- ls(env=evbase)

### copy variables from evbase to local environment 
  evei <- getEnvVar(evbase, environment())  ##environment
### we need these assigments so that R CMD ceck does not complaint
### The code works without them as getEnvVar is assigning them in this environment anyway
  if(exists("tomogPct")) tomogPct <- get("tomogPct", env=evei)
  if(exists("eigraph.bvsmth")) eigraph.bvsmth <- get("eigraph.bvsmth", env=evei)
  if(exists("eigraph.t")) eigraph.t <- get("eigraph.t", env=evei)
  if(exists("eigraph.loess")) eigraph.loess <- get("eigraph.loess", env=evei)
  if(exists("eigraph.bblo")) eigraph.bblo <- get("eigraph.bblo", env=evei)
  if(exists("eigraph.bbhi")) eigraph.bbhi <- get("eigraph.bbhi", env=evei)
  if(exists("eigraph.bwlo")) eigraph.bwlo <- get("eigraph.bwlo", env=evei)
  if(exists("eigraph.bwhi")) eigraph.bwhi <- get("eigraph.bwhi", env=evei)
  if(exists("eigraph.bb")) eigraph.bb <- get("eigraph.bb", env=evei)
  if(exists("eigraph.bw")) eigraph.bw <- get("eigraph.bw", env=evei)
  if(exists("eigraph.tit")) eigraph.tit <- get("eigraph.tit", env=evei)
####end of R CMD complaints
  if(any(!is.na(psiu))) eigraph.psiu <- psiu
 
  if(exists("EnumTol")) tol <- EnumTol <- get("EnumTol", env=evei)
  else tol <- EnumTol <- get("EnumTol", env=evbase)
  if(class(tol)=="try-error") tol <- 0.0001
  ### fill the titles
  str <- tolower(str)
  #### f=checkginputs(dbuf);
  if(identical(str,"tomogd")){###	@ tomograph with data only  @
   
    if(!exists("eigraph.psiu")|| !exists("psiu"))
      psiu <- eigraph.psiu <- get("eigraph.psiu", env=evbase)

    t <- eiread(dbuf,"t")
    x <- eiread(dbuf,"x")
    f <- rows(x)
    i <- order(runif(f))
    i <- i[1:floor(f* eigraph.smpl)]
    tst <- eiread(dbuf, "Eselect")
    if(length(tst) <= 0 || is.na(tst)) tst <- get("Eselect", env=ev)
    if(scalone(tst)) tst <- matrix(1,nrow=f,ncol=1)

    tomog(x[i],t[i],eigraph.psiu,tst[i],tol)
   
    return("tomogd")
  }

  
    if(identical(str,"tomog")){	### @ tomography with ML contours  @
      if(exists("EnonPar")) EnonPar <- get("EnonPar", env=evei)
     
      if(as.logical(EnonPar)|| EnonPar >= 1)
        stop("ei: Graph 'tomog' set for parametric estimations only")     
      if ("psiu" %in% nm) 
        psiu <- eigraph.psiu <- dbuf[["psiu"]]
      else
        psiu <- eigraph.psiu <- eiread(dbuf,"psiu")
      
      eigraph(dbuf,"tomogd",psiu=psiu,eigraph.tit="Tomography with ML contours",...)
           
      return("tomog")
    }
   if(identical(str,"tomogp")) { ###		@ tomography with mean post contours @
   
     if ("mppsiu" %in% nm) 
       mppsiu <- eigraph.psiu <- dbuf[["mppsiu"]]
     else
       mppsiu <- eigraph.psiu <- eiread(dbuf,"mppsiu")
     if(all(is.na(mppsiu))) message("ei: Graph 'tomogp' set for parametric estimations") 
     eigraph(dbuf,"tomogd",psiu=mppsiu,eigraph.tit="Tomography with mean post contours",...)
     return("tomogp")
   }
  if(identical(str,"tomoge")) { ### 	@ tomography w/ Estimated betab,betaw @
     if ("betab" %in% nm) 
       betab <- dbuf[["betab"]]
     else
       betab <- eiread(dbuf,"betab")
     if ("betaw" %in% nm) 
       betaw <- dbuf[["betaw"]]
     else
       betaw <- eiread(dbuf,"betaw")
       
     eigraph(dbuf,"tomogd",eigraph.tit="Tomography with data and estimated betab, betaw",...)
     points(as.vector(betab), as.vector(betaw),pch=21)
       if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct) 
     return("tomoge")
   }

  if(identical(str,"tomogt")) { ### 	@ tomography w/ true betab,betaw @
    ### I do not know what is truth check with Gary
     if (!("truth" %in% nm)){
       warning( "ei: 'truth' must be stored first")
       return(NA)
      
     }
       betab <- eiread(dbuf,"truthb")
       betaw <- eiread(dbuf,"truthw")
           
     eigraph(dbuf,"tomogd",eigraph.tit="Tomography with data and estimated betab, betaw",...)
     points(as.vector(betab), as.vector(betaw),pch=21)
     if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct)
     return("tomogt")
   }
   if(identical(str,"tomogci")) {    ###	@ tomography with 80% conf intervals @
     b <- eiread(dbuf,"ci80bw")
     if(scalmiss(b)) {
       eigraph(dbuf,"tomogd",...)
       warning("No confidence intervals")
       return(NULL) 
     }
     tit <- "Tomography with 80% CI's"
###the confidence interval lines
     conf.intv(b,tit,npts=500,dbuf)
    
     if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct)  
     return("tomogci")
   }
    if(identical(str,"tomogci95")) {    ###	@ tomography with 95% conf intervals @
     b95 <- eiread(dbuf,"ci95bw")
     if(scalmiss(b95)) {
       eigraph(dbuf,"tomogd",...)
       warning("No confidence intervals 95")
       return(NULL) 
     }
     tit <- "Tomography with 95% CI's"
###the confidence interval lines
     conf.intv(b95,tit,npts=500,dbuf)
     if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct)
     return("tomogci95")
   }
 if(identical(str,"tomogs")){ ####  @ tomog,tomogp,tomogCI,tomog(or tomogE)@
    
   op <- par(no.readonly=TRUE)
   if(scalzero(eiread(dbuf,"Enonpar"))) par(mfrow=c(2,2))##2 rows and 2 columns
   else  par(mfrow=c(1,2))   ##1 row and 2 columns
 
   if(scalzero(eiread(dbuf,"Enonpar"))){
     tit <- "ML Contours" 
     eigraph(dbuf,"tomog",...)
     tit <- "Mean Posterior Contours"
     eigraph(dbuf,"tomogp",...)
   }
   tit <- "80% Confidence Intervals"
   eigraph(dbuf,"tomogci")
   if( vin(dbuf,"truth")){
     tit <- "True coordinates"
     eigraph(dbuf,"Tbivar",eigraph.tit=tit,...)
   }else{
     tit <- "Point estimates"
     eigraph(dbuf,"estsims",eigraph.tit=tit,...)
   }
   
   par(op) ### reset to default values
   assign("eigraph.tit","", env=evbase)
   return("tomogs")
    }
  if(identical(str,"nonpar")){ ###	@ tomogd & nonparametric contour & surface  @
 op <- par(no.readonly=TRUE)    
    Enumtol <- as.vector(get("EnumTol", env=evbase))
    x <- eiread(dbuf,"x")
    t <- eiread(dbuf,"t")
    if(vin(dbuf,"Eselect")){
      Eselect <- eiread(dbuf,"Eselect")
      if (rows(Eselect)==rows(x)){
        x <- subset(x,as.logical(Eselect))
        t <- subset(t,as.logical(Eselect))
      }
    }
    e <- cbind((x < Enumtol), (x>(1-Enumtol)))
    v <- rbind(Enumtol,(1-Enumtol))
    x <- recode(x,e,v)
    e <- cbind((t < Enumtol), (t>(1-Enumtol)))
    t <- recode(t,e,v)
    eigraph.eval <- get("eigraph.eval", env=evbase)
    tt <- as.matrix(seqase(0,1,eigraph.eval))
    o <- matrix(1, nrow=eigraph.eval,ncol=1)
    ka = kronecker(tt,o)
    kb = kronecker(o,tt)
    a <-t(reshape(ka,eigraph.eval,eigraph.eval))  ###   @ px @
    b <- t(reshape(kb,eigraph.eval,eigraph.eval))  ### @ py @
   
    z <- nonbiv(t,x,a,b,evbase=evbase,
                eigraph.bvsmth=eigraph.bvsmth,Enumtol=Enumtol)  ###                              @ pz @
    op <- par(no.readonly=TRUE)
  

###  x11()--OR--- par()
    par(mfrow=c(1,2)) 
    z <- z%dot/%maxc(as.vector(z))
    if(scalone(eigraph.dbuf)){
      eigraph.dbuf <- as.list(eigraph.dbuf)
      eigraph.dbuf <- c(eigraph.dbuf,px=list(t(tt)))
      eigraph.dbuf <- vput(eigraph.dbuf,tt,"py")
      eigraph.dbuf <- vput(eigraph.dbuf,z,"pz")
    }
    dat <- data.frame(dataz=list(z),datax=list(tt),datay=list(tt))
    dataz <- z
   
    contour(tt,tt,z, main= "Nonparametric Contours", xlab="betab",ylab="betaw",
           col=rainbow(length(tt)))
###  lattice: contourplot(dataz,main=tit,xlab="betab",ylab="betaw",
###          col=rainbow(length(tt)))
    tit <- "Nonparametric Surface"
#### lattice: same as wireframe(z,tt,tt)
    par(mfrow=c(1,2),new=TRUE, fin=c(6,6),plt=c(0.465,1,0,1))
    persp(tt,tt,z,theta=-30,phi=30,col=rainbow(length(tt)),
          main="Nonparametric Surface",ticktype="detailed",xlab="betab", ylab="",
          zlab="",nticks=5,pty="m")
###    mtext("betaw",side=4)
    mtext(tit,side=3,cex=1.25,font=2)
    mtext("betaw", side=2)
    par(op) ### reset to default values
    return("nonpar")
  }
  if (identical(str,"xtc")||identical(str,"xtbubble")){ ###			@ basic x by t with sized circ's @ 
    x <- eiread(dbuf,"x")
    t <- eiread(dbuf,"t")
    prop <- 1
    eigraph.circ0 <- eigraph.circ
    if (eigraph.circ==0){
      tvap <- prop*eiread(dbuf,"n")
      if(all(is.na(tvap))) tvap <- eiread(dbuf,"tvap")
      eigraph.circ <- 20*(tvap-minc(tvap)+1)%dot/%(1+maxc(tvap)-minc(tvap)) 
      if (stdc(tvap)<50)
	  eigraph.circ <- 0.5*eigraph.circ 
    }
    y <- t
    if(identical(eigraph.tit, "") )
      eigraph.tit <- "Basic X by T with sized circles"
   
   
    xst <- ifelse(min(x) > 0, 0,min(x))
    yst <- ifelse(min(y) > 0, 0,min(y))
    xend <- ifelse(max(x) <1, 1,max(x))
    yend <- ifelse(max(y) <1, 1, max(y))
    plot(x,t,cex=eigraph.circ*eigraphC,xlab=eigraph.x, 
         ylab=eigraph.t, main=eigraph.tit,xlim=c(xst,xend), ylim=c(yst,yend))
    assign("eigraph.tit","",env=evbase)
    assign("eigraph.circ",eigraph.circ0,env=evbase)
   return("xtbubble")
  }
  if(identical(str,"xt")){
    eigraph.circ0 <- eigraph.circ
    eigraph(dbuf,"xtc",eigraph.circ=0.5,eigraph.tit="Basic X by T",...)
 ###   mtext("Basic X by T", side=3)
    assign("eigraph.tit","",env=evbase)
    assign("eigraph.circ",eigraph.circ0,env=evbase)
     
    return("xt")
  }
   if (identical(str,"xgraph")){###		@ X-graph  @
    eigraph.circ0 <- eigraph.circ
    eigraph(dbuf,"xtc",eigraph.circ=0.5,eigraph.tit="X-graph by T",...)
    abline(a=0.,b=1)   ###   pline(0,0,1,1);
    abline(a=1, b=-1)  ###   pline(0,1,1,0);
    assign("eigraph.tit","",env=evbase)
    assign("eigraph.circ",eigraph.circ0,env=evbase)
    return("xgraph")
  }
   if( identical(str,"xgraphc")||  identical(str,"xgraphbubble")){###		@ X-graph with sized circles  @
     eigraph.circ0 <- eigraph.circ
     eigraph(dbuf,"xtc",eigraph.tit="X-graph with sized circles",...)
     abline(a=0.,b=1) ### pline(0,0,1,1);
     abline(a=1.,b=-1) ###pline(0,1,1,0);
     assign("eigraph.tit","",env=evbase)
     assign("eigraph.circ",eigraph.circ0,env=evbase)
   return("xgraphbubble")
  }
  if(identical(str,"goodman")){	### @ x by t plot w/ goodman's line @
    eigraph.circ0 <- eigraph.circ
    eigraph(dbuf,"xtc",eigraph.circ=0.5,eigraph.tit="X by T plot w/ Goodman's line",...)
    b <- eiread(dbuf,"goodman")
    pline(0,b[2,1],1,b[1,1], new =TRUE,color="red")
    assign("eigraph.tit","",env=evbase)
    assign("eigraph.circ",eigraph.circ0,env=evbase)
    return("goodman")
  }
  if(identical(str,"xtfit")){ ###			@ x by t with Exp and 80%CI's  @
    b <- eiread(dbuf,"expvarcis") 
    if(scalmiss(b)) {
      warning("ei: Graph 'xtfit' is not set for non-parametric estimations")
      return(NA)
    }
### @ CI's drawn only if no Z's are specified @
    
    if(identical(eigraph.tit,""))
      eigraph(dbuf,"xtc",eigraph.tit="X by T with Exp and 80%CI's",...)
    else
       eigraph(dbuf,"xtc")
    
    plines(b[,1],b[,3])
    plines(b[,1],b[,2],color="blue")
    plines(b[,1],b[,4],color="blue")
 
    if(eigraph.loess){
 	b <- eiread(dbuf,"expvarci")

	plines(b[,1],b[,3])
     plines(b[,1],b[,2],color="red")
	plines(b[,1],b[,4],color="red")
   }
    assign("eigraph.tit","",env=evbase)
    return("xtfit")
  }
  if (identical(str,"xtfitg")){###		@ xtfit with goodman's regression @
    b <- eiread(dbuf,"goodman")
  ###  assign(eigraph.tit,"X by T with Exp & 80%CI's and Goodman",env=evbase)
    y <- eigraph(dbuf,"xtfit",eigraph.tit="X by T with Exp & 80%CI's and Goodman",...)
    if(is.na(y)){
       warning("ei: Graph 'xtfit' is not set for non-parametric estimations")
      return(NA)
    }
    pline(0,b[2,1],1,b[1,1],color="red")
    assign("eigraph.tit","",env=evbase)
    return("xtfitg")
  }
  if(identical( str,"fit")){ ###		        @ tomogp and xtfit @
    op <-  par(no.readonly=TRUE)

    if(vin(dbuf,"truth")) par(mfrow=c(2,1))
     
    y <- eigraph(dbuf,"xtfit",...)
    if(is.na(y)) {
      warning("ei: Graph 'fit' is not set for non-parametric estimations")
      return(NA)
    }
    if(vin(dbuf,"truth") )eigraph(dbuf,"tomogp",...)
    else
      warning("truth is not stored")
    par(op)
    return("fit")
  }
   if(identical(str,"fitt")){	###	        @ tomog and mlfit  @
     op <-  par(no.readonly=TRUE)
     if(vin(dbuf,"truth")) par(mfrow=c(2,1))
   
     y <- eigraph(dbuf,"xtfit",...)
     if(is.na(y)) {
       warning("ei: Graph 'fit' is not set for non-parametric estimations")
       return(NA)
     }
   
    if (vin(dbuf,"truth"))
   
      eigraph(dbuf,"tomogt",...)
     else
       warning("truth is not stored")
    par(op)
    return("fitt")
  }

  if(identical(str,"profile")){ ###               @ profile plot of phi @
    op <-  par(no.readonly=TRUE)
    profileit(dbuf,0)###, eigraph.pro=c(list(c(rep(-1,5),0,0)*15), list(c(10,20,rep(10,3),0,0))))
    par(op)
    return("profile")
  }
 
  if(identical(str,"profiler")){ ###               @ profile plot of R (cdfbvn) @
    op <-  par(no.readonly=TRUE)
    profileit(dbuf,1)###,eigraph.pro=c(list(c(rep(-1,5),0,0)*5), list(c(rep(1,2)*5,rep(1,3),0,0))))
    par(op)
    return("profileR")
  }
 
 
   if( (postb <- identical(str,"postb")) ||(postw <- identical(str,"postw")) ){ ### @ posterior of dist agg B^b and B^w  @
     a <- eiread(dbuf,"abounds")
    strt <- ifelse(postb, a[1,1],a[1,2])
    endd <- ifelse(postb,a[2,1],a[2,2])
   
    assign("strt",strt,env=evbase)
    assign("endd",endd,env=evbase)
    kern <- get("kern",env=evbase) 
    kern <- as.vector(kern)
    whiskr <- 0
    if (vin(dbuf,"truth")){
      b <- eiread(dbuf,"aggtruth")
      if(postb)
        pline(b[1],0,b[1],3)
      else
        pline(b[2],0,b[2],3) 

    }

    output <- 0
    assign("output", output,env=evbase)
    b <- eiread(dbuf,"aggs") 
    if(postb)
      lst <- dens(b[,1],evbase)
    else
      lst <- dens(b[,2],evbase)
    a <- lst[[1]]
    b <- lst[[2]]
    a <- rbind(a[1],a,a[rows(a)])
    b <- rbind(0,b,0)

###   ylabel("density, f(B^w)");
###   xlabel(_eigraph_bw);
    if(postb){
      eigraph.x <-get("eigraph.bb", env=evbase)
      tit <- paste("Posterior dist aggregate B^b; kern=",kern)
    }else{
      eigraph.x <- get("eigraph.bw", env=evbase)
      tit <- paste("Posterior dist aggregate B^w; kern=", kern)
    }
    if(output==1)  par(new=TRUE)
    lst <-listwis2(a,b) 
    plot(lst[[1]],lst[[2]],type="l",xlab=eigraph.x,ylab="Density",
         main=tit)
 
   if(postb) return("postb")
   else return("postw")
  }
  if(identical(str,"post")){ ### @ posterior of dist aggs B^b B^w @
    op <-  par(no.readonly=TRUE)
    par(mfrow=c(1,2))
    eigraph(dbuf,"postb",...)
    eigraph(dbuf,"postw",...)
    par(op)
    return("post")
  }
 if((identical( str,"betab")) ||
    (identical(str,"betaw"))){ ###		@ density est of beta^b  @
 
   betb <- identical( str,"betab")
   strt <- 0
   endd <- 1
   kern <- get("kern",env=evbase)
   kern <- as.vector(kern)
   assign("strt",strt,env=evbase)
   assign("endd",endd,env=evbase)
   
   if(betb)
     bet <- na.omit(eiread(dbuf,"betaB"))
   else
     bet <- na.omit(eiread(dbuf,"betaW"))
    output <- 0
    assign("output", 0, env=evbase)
    lst <- dens(bet,evbase) 
    a <- lst[[1]]
    b <- lst[[2]]
    o <- matrix(1, nrow=rows(bet),ncol=1)
    a <- rbind(a[1],a,a[rows(a)])
    b <- rbind(0,b,0)

###    _pline=o~ (o*6)~ betaB~ (o*0)~ betaB~ (o*(maxc(b)/15))~ o~ (o*15)~(o*0);  
###   ylabel("density across precincts, f(betaB)");
    ylabel <- "Density across precints" 
##    xlabel(_eigraph_bb) 
    xlabel <- ifelse(betb, get("eigraph.bb", env=evbase),
                     get("eigraph.bw", env=evbase))
    tit1 <- paste("Density est of beta^b; kern=",kern)
    tit2 <- paste("Density est of beta^w; kern=",kern) 
    tit <- ifelse(betb, tit1,tit2)

###    xtics(_eigraph_bblo,_eigraph_bbhi,(_eigraph_bbhi-_eigraph_bblo)/4,5);

###    gosub setup;
       
    lst <-listwis2(a,b) 
    plot(lst[[1]],lst[[2]],type="l",xlab=xlabel,ylab=ylabel,
         main=tit)
 
   if(betb)
     return("betab")
   else
     return("betaw")

}
  if(identical(str,"beta")){###	@ density est of betab and betaw @
    op <-  par(no.readonly=TRUE)
    par(mfrow=c(1,2))
    eigraph(dbuf,"betab",...)
    eigraph(dbuf,"betaw",...)
    par(op)
    return("beta")  
}
  if(identical(str,"results")){###   @ posterior of dist aggs B^b B^w, hist of beta's @
    op <-  par(no.readonly=TRUE)
    par(mfrow=c(2,2))
    eigraph(dbuf,"postb",...)
    eigraph(dbuf,"postw",...)
    eigraph(dbuf,"betab",...)
    eigraph(dbuf,"betaw",...)
    par(op)
    return("results")
  }

  if(identical( str,"movie")){ ###	     @ moving picture of betaB and betaW densities  @
      op <-  par(no.readonly=TRUE)
      
    Esims <- eiread(dbuf,"Esims")
    assign("strt", 0, env=evbase)
    strt <- 0
    assign("endd", 1, env=evbase)
    endd <- 1
    kern <- get("kern",env=evbase)
    kern <- as.vector(kern)
    output <- 0
    assign("output", 0, env=evbase)

    x <- eiread(dbuf,"x")
    t <- eiread(dbuf,"t")
    n <- eiread(dbuf,"n")
    p <- rows(x)
    betaB <- eiread(dbuf,"betaBs")
    betaW <- eiread(dbuf,"betaWs")
    i <- 1

    while( i<=p){
      cc <- betaB[i,]
      d <- betaW[i,]

      par(mfrow=c(2,2))
    
 ####fst plot:betaB     
	if( !(ismiss(cc))){

       lst <- dens(cc, evbase)
       a <- lst[[1]]
       b <- lst[[2]]
	  o <- matrix(1, nrow=Esims,ncol=1)
	  a <- as.matrix(c(a[1],a,a[rows(a)]))
	  b <- matrix(c(0,b,0))
###	  _pline=o~ (o*6)~ c~ (o*0)~ c~ (o*(maxc(b)/15))~ o~ (o*15)~(o*0);

	  ylabel <- "density, f(betaB)"
	  xlabel <- get("eigraph.bb",env=evbase)

       lst <- listwis2(a,b)
       plot(lst[[1]], lst[[2]], xlab=xlabel, ylab=ylabel,type="h",
            main="Moving picture of betaB density") 
}
      
###scnd plot:betaW     
	if (!(ismiss(d))){

       lst <- dens(d,evbase)
       a <- lst[[1]]
       b <- lst[[2]]
	  o <- matrix(1,nrow=Esims,ncol=1)
	  a <- rbind(a[1],a,a[rows(a)])

	  b <- as.matrix(c(0,b,0))

       ##	  _pline=o~ (o*6)~ d~ (o*0)~ d~ (o*(maxc(b)/15))~ o~ (o*15)~(o*0);

	  ylabel <- "density, f(betaW)"

	  xlabel <- get("eigraph.bw", env=evbase)

      lst <-listwis2(a,b) 
	 plot(lst[[1]], lst[[2]], type="h", xlab = xlabel, ylab=ylabel,
           main="Moving picture of betaW density")
      }


####thrd plot:simulations
      if (!(ismiss(c(cc,d)))){

	  xlabel <- get("eigraph.bb",env=evbase)
	  ylabel <- get("eigraph.bw",env=evbase)
	  title <- "Simulations"
       lst <- listwis2(cc,d)
 
	  plot(lst[[1]], lst[[2]], type="l", xlab=xlabel, ylab=ylabel, main=title)

     }
      
      plot(-1:1,-1:1, xlab="", ylab="", type="n",
           axes=FALSE, xaxt="n",yaxt="n")
      text(-0.25,0.25,paste("Observation =",i))
      text(0.19, 0,   paste("X =", x[i]),adj=c(0.5,0))
      text(0.2, -0.25,   paste("T =", t[i]),adj=c(0.5,0))
      text(0.21, -0.5,   paste("N =", n[i]),adj=c(0.5,0))        

      cat("Observation = ",i,"\nX = ", x[i],"\nT = ",t[i],"\nN = ", n[i])
      
      answer <- user.prompt(1) 
      answer <- trim.blanks(answer)
      bool <- FALSE
      if(!identical(answer,""))
        bool <- chk.numeric(answer)
      ob <- 0
      if(bool) ob <- as.numeric(answer)
      if(ob>=1 && ob <= p) i <- ob
      else  i <- i+1

    }
      par(op)
   return("movie")
  }
    

  if(identical(str,"movied")){
    ###@ pick out individual lines on tomography plot @
  op <-  par(no.readonly=TRUE)
    output <- 0
    assign("output", 0, env=evbase)
    x <- eiread(dbuf,"x") 
    t <- eiread(dbuf,"t") 
    n <- eiread(dbuf,"n") 
    p <- rows(x) 
    i <- 1

    while( i<=p){

      title <- paste("Dark line is observation",i)
      xlabel <- get("eigraph.bb", env=evbase)
      ylabel <- get("eigraph.bw", env=evbase)
      f <- sortind(rndu(p,1))
      eigraph.smpl <- get("eigraph.smpl", env=evbase)
      f <- f[1:floor(p*eigraph.smpl)]
      tst <- eiread(dbuf,"Eselect") 

      if (scalone(Eselect))
        tst <- matrix(1,nrow=p,ncol=1)

      tst[i] <- 2
      eigraph.psiu <- get("eigraph.psiu", env=evbase)
     
      tomog(x[f],t[f],eigraph.psiu,tst[f]) 
      mtext(paste("Observation:", i,"; X=",x[i],"; T=",t[i],"; N=", n[i]), side=3)
      answer <- user.prompt(1) 
      answer <- trim.blanks(answer)
      bool <- FALSE
      if(!identical(answer,""))
        bool <- chk.numeric(answer)
      ob <- 0
      if(bool) ob <- as.numeric(answer)
      if(ob>=1 && ob <= p) i <- ob
      else  i <- i+1

    }
  par(op)
  return("movied")
  }

 if( prectb <- identical(str,"prectb")){ ###	@ est by true of betaB_i  @
    if(! vin(dbuf,"truth")){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }
    
    b <- eiread(dbuf,"truthB")
    a <- eiread(dbuf,"betaB")
    xlabel <- paste("ESTIMATED",get("eigraph.bb", env=evbase))
    ylabel <- paste("TRUE",get("eigraph.bb", env=evbase))
    cc <- eiread(dbuf,"bvap") 
    evscale <- 0.25
    cc <- evscale*get("eigraphC", env=evbase)*24.5*cc/maxc(cc) 
    triple(a,b,cc, xlabel=xlabel, ylabel=ylabel,
           title="Est by true of betaB_i")
    bblo <- get("eigraph.bblo", env=evbase)
    bbhi <- get("eigraph.bbhi",env=evbase)
   
    d <- eiread(dbuf,"ci80b") 
    e <- meanc(na.omit(abs(d-a)))

    pline(bblo,bblo,bbhi,bbhi,negslope=FALSE,col="red")
    pline(bblo,bblo+e[1],bbhi,bbhi+e[1],col="blue", negslope=FALSE)
    pline(bblo,bblo-e[2],bbhi,bbhi-e[2], negslope=FALSE,col="blue")
    return("prectb")
  }
  if(identical(str,"prectw")){ ###	@ est by true of betaW_i  @
    if(! vin(dbuf,"truth")){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }
    
    a <- eiread(dbuf,"betaW")
    b <- eiread(dbuf,"truthW")
    xlabel <- paste("ESTIMATED",get("eigraph.bw", env=evbase))
    ylabel <- paste("TRUE",get("eigraph.bw", env=evbase))
    cc <- eiread(dbuf,"wvap")
    evscale <- 0.25
    cc <- evscale*get("eigraphC", env=evbase)*24.5*cc/maxc(cc) 
    triple(a,b,cc, xlabel=xlabel, ylabel=ylabel,
                   title="Est by true of betaW_i")
    bwlo <- get("eigraph.bwlo", env=evbase)
    bwhi <- get("eigraph.bwhi",env=evbase)

    d <- eiread(dbuf,"ci80w") 
    e <- meanc(na.omit(abs(d-a)))

    pline(bwlo,bwlo,bwhi,bwhi,negslope=FALSE,col="red")
    pline(bwlo,bwlo+e[1],bwhi,bwhi+e[1],col="blue", negslope=FALSE)
    pline(bblo,bblo-e[2],bbhi,bbhi-e[2], negslope=FALSE,col="blue")
   return("prectw")
  }
 if(identical(str,"truth")){ ###			@ compare truth to estimates @
   op <-  par(no.readonly=TRUE)
    if( !(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }

    par(mfrow=c(2,2))
    eigraph(dbuf,"postb")
    eigraph(dbuf,"postw")
    eigraph(dbuf,"prectb")
    eigraph(dbuf,"prectw")
    par(op)
    return("truth")
  }
  if(identical(str,"lines")){ ###			@ xt with one EST'd per precinct  @

    betab <- eiread(dbuf,"betab")
    betaw <- eiread(dbuf,"betaw")
    lst <- listwis2(betab,betaw)
    betab <- lst[[1]]
    betaw <- lst[[2]]
    a <- rows(betab)
    z <- matrix(0,nrow=a,ncol=1)
    o <- matrix(1,nrow=a,ncol=1)
    
    eigraph(dbuf,"xtc",eigraph.tit="XT with one EST'd per precinct",...)
    pline(z,betab,o,betaw)
    assign("eigraph.tit","",env=evbase)
    return("lines")
  }


  if(identical(str,"tlines")){###		@ xt with one TRUE line per precinct @

    if (!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }
     
    betab <- eiread(dbuf,"truthb")
    betaw <- eiread(dbuf,"truthw") 
    lst <- listwis2(betab,betaw)
    betab <- lst[[1]]
    betaw <- lst[[2]]
    a <- rows(betab) 
    z <- matrix(0,nrow=a,ncol=1)
    o <- matrix(1, nrow=a,ncol=1)
    eigraph.circ0 <- eigraph.circ
   
    eigraph(dbuf,"xtc",eigraph.tit="XT with one TRUE line per precinct",
            eigraph.circ=1,...)
    pline(z,betab,o,betaw)
    assign("eigraph.tit","",env=evbase)
    assign("eigraph.circ",eigraph.circ0,env=evbase)
    ###Gauss  pline(z,betaw,o,betab,negslope=FALSE)
    return("tlines")
  }
  
  if(identical(str,"bivar")){	###		@ EST'd betab by betaw @
    
    betaB <- eiread(dbuf,"betab") 
    betaW <- eiread(dbuf,"betaw") 
    xlabel <- get("eigraph.bb", env=evbase)
    ylabel <- get("eigraph.bw", env=evbase)
    lst <- listwis2(betaB,betaW)
    betaB <- lst[[1]]
    betaW <- lst[[2]]
    psymsiz <- .5*get("eigraphC",env=evbase)
   
    plot(betaB, betaW, xlab=xlabel, ylab=ylabel,xlim=c(0,1), ylim=c(0,1),
         main="EST'd betab by betaw", type="p", cex=psymsiz)
 
    pline(0,0,1,1,negslope=FALSE)
    if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct) 

    return("bivar")
}
  
  if (identical(str,"tbivar")){ ###		@ TRUE betab by betaw @

    if(!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }
    betaB <- eiread(dbuf,"truthb")
    betaW <- eiread(dbuf,"truthw")
    lst <- listwis2(betaB,betaW)
    betaB <- lst[[1]]
    betaW <- lst[[2]]
    xlabel <- paste("TRUE ", get("eigraph.bb", env=evbase))
    ylabel <- paste("TRUE ",get("eigraph.bw",env=evbase))
    psymsiz=.5*get("eigraphC", env=evbase)
   
    plot(betaB, betaW, xlab=xlabel, ylab=ylabel,xlim=c(0,1), ylim=c(0,1),
         main="TRUE betab by betaw", type="p", cex=psymsiz)
 
    pline(0,0,1,1,negslope=FALSE) 
    if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct)  
    return("tbivar")
  }

 if(identical( str,"betabw")){###		@ lines,tlines,bivar,tbivar @
      op <-  par(no.readonly=TRUE)
     if(!(vin(dbuf,"truth")))
       par(mfrow=c(1,2))
     else
       par(mfrow=c(2,2))

      eigraph(dbuf,"lines",...)
      eigraph(dbuf,"bivar",...)
     if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct)  

     if(vin(dbuf,"truth")){

       eigraph(dbuf,"tlines",...)
       
       eigraph(dbuf,"tbivar",...)
       if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct)
   }
      par(op)
     return("betabw")
   }

  if(identical(str,"ptileb")){	###	@ true percentile for betaB @

    if (!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }

    xlabel <- get("eigraph.bb", env=evbase)
    ylabel <- "True Percentile"
    betaB <- eiread(dbuf,"betab")
    a <- eiread(dbuf,"truPtile")
    bvap <- eiread(dbuf,"bvap")
    evscale <- 0.3
    cc <- evscale*get("eigraphC", env=evbase)*24.5*bvap/maxc(bvap)
    
    triple(betaB,a[,1],cc,xlabe=xlabel,ylabel=ylabel,title="True percentile for betaB",pxscale=0,pyscale=0)

    pline(0,0.25,1,0.25,negslope=FALSE)
    pline(0,0.5,1,0.5,negslope=FALSE)
    pline(0,0.75,1,0.75,negslope=FALSE)
    return("ptileb")
  }

  if(identical(str,"ptilew")){###		@ true percentile for betaW @

    if(!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }
    xlabel <- get("eigraph.bw", env=evbase)
    ylabel <- "True Percentile"

    betaW <- eiread(dbuf,"betaw")
    a <- eiread(dbuf,"truPtile")

    wvap <- eiread(dbuf,"wvap")
    evscale <- 0.3
    cc <- evscale*get("eigraphC", env=evbase)*24.5*wvap/maxc(wvap)

 
    triple(betaW,a[,1],cc,xlabe=xlabel,ylabel=ylabel,title="True percentile for betaW",pxscale=0,pyscale=0)
    pline(0,0.25,1,0.25,negslope=FALSE)
    pline(0,0.5,1,0.5,negslope=FALSE)
    pline(0,0.75,1,0.75,negslope=FALSE)
    return(,"ptilew")
  }

  

  if(identical(str,"ptile")){###			@ ptileb & ptilew @

    if(!(vin(dbuf,"truth"))){
      stop("ei: 'truth' must be stored first")
      return(NA)
    }
    op <-  par(no.readonly=TRUE)
    par(mfrow=c(1,2))
    eigraph(dbuf,"ptileb",...)
    eigraph(dbuf,"ptilew",...)
    par(op)
    return("ptile")
  }

  if(simb <- identical(str,"simsb")){ ###			@ sims of betab by betab  @

    if(!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }

    psymsiz <- 0.1*get("eigraphC", env=evbase)
    xlabel <- paste(get("eigraph.bb", env=evbase),"simulations")
    ylabel <- paste("True", get("eigraph.bb",env=evbase))

    betab <- eiread(dbuf,"betabs")
    a <- eiread(dbuf,"truthb") 
    lst <- listwis2(betab,a)
    betab <- lst[[1]]
    dm <- dim(betab)
    a <- lst[[2]]
  ###  aa <- matrix(a,nrow=dm[1],ncol=dm[2])
    mxb <- max(betab)
    mnb <- min(betab)
    phcev <- 0:18
    cl <- ifelse(dim(a)[2]>=dim(betab)[2],dim(a)[2],dim(betab)[2])
    phcev <- rep(phcev, floor(cl/length(phcev))+1)
   
    matplot(betab, a, xlab=xlabel, ylab=ylabel,pch=phcev,
         main="Sims of betab by betab true", type="p", cex=psymsiz)
 
    pline(0,0,1,1, negslope=FALSE)
    return("simsb")
  }
  
  if(simw <- identical(str,"simsw")){ ###			@ sims of betaw by betaw  @

    if(!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }

    psymsiz <- 0.1*get("eigraphC", env=evbase)
    xlabel <- paste(get("eigraph.bw", env=evbase),"simulations")
    ylabel <- paste("True", get("eigraph.bw",env=evbase))

    betaw <- eiread(dbuf,"betaws")
    a <- eiread(dbuf,"truthw") 
    lst <- listwis2(betaw,a)
    betaw <- lst[[1]]
    a <- lst[[2]]
  ###  aa <- matrix(a,nrow=dm[1],ncol=dm[2])
    phcev <- 0:18
    cl <- ifelse(dim(a)[2]>=dim(betaw)[2],dim(a)[2], dim(betab)[2])
    phcev <- rep(phcev, floor(cl/length(phcev))+1)
      
    matplot(betaw, a, xlab=xlabel, ylab=ylabel,pch=phcev,
         main="Sims of betaw by betaw true", type="p", cex=psymsiz)
 
    pline(0,0,1,1, negslope=FALSE)
    return(,"simsw")
  }

    if(identical(str,"sims")){ ###			@ simsB & simsW @

    if (!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }
    op <-  par(no.readonly=TRUE)
    par(mfrow=c(1,2))
    
    eigraph(dbuf,"simsb",...)

    eigraph(dbuf,"simsw",...)
    par(op)
    return("sims")
  }
  if(identical(str,"estsims")){ ###	@ sim'd est betab by est betab  @
          
      tit <- "Sim\'d est betaw by est\'d betab"
      betab <- eiread(dbuf,"betabs")
      betaw <- eiread(dbuf,"betaws")
      lst <- listwis2(betab,betaw)
### same as matplot((lst[[2]],lst[[1]],...) 
      plot(lst[[2]]~lst[[1]], type="p",col="blue", 
           xlim=c(eigraph.bblo,eigraph.bbhi), ylim=c(eigraph.bwlo,eigraph.bwhi),
           xlab= eigraph.bb,ylab= eigraph.bw, main=tit)
      
      if ("betab" %in% nm) 
        betab <- dbuf[["betab"]]
      else
        betab <- eiread(dbuf,"betab")
      if ("betaw" %in% nm) 
        betaw <- dbuf[["betaw"]]
      else
        betaw <- eiread(dbuf,"betaw")
      
      points(as.vector(betab), as.vector(betaw),col="red",pch=20,cex=1.35)
      
      if ("psiu" %in% nm) 
        psiu <- eigraph.psiu <- dbuf[["psiu"]]
      else
        psiu <- eigraph.psiu <- eiread(dbuf,"psiu")
      
      if(!is.na(psiu)&& length(psiu) >1) MLKcontours(psiu,tomogPct)
     
      return("estsims")
    }
  
  if(identical(str,"biasb")){###			@ x by betaB @

    psymsiz <- 0.3*get("eigraphC", env=evbase)
    xlabel <- get("eigraph.x",env=evbase)
    ylabel <- paste("EST'd",get("eigraph.bb",env=evbase))

    x <- eiread(dbuf,"x")

    betab <- eiread(dbuf,"betab")
    a <- eiread(dbuf,"eaggbias") 
    lst <- listwis2(x,betab)
    x <- as.vector(lst[[1]])
    betab <- lst[[2]]
    plot(x,betab, xlab=xlabel, ylab=ylabel,
         main="X by betaB", type="p", cex=psymsiz)
    pline(0,a[1,1],1,a[1,1]+a[2,1],negslope=FALSE)
    return("biasb")
  }


  if(identical( str,"biasw")){###			@ x by betaW @

     psymsiz=0.3*get("eigraphc", env=evbase)

    xlabel <- get("eigraph.x", env=evbase)
    ylabel <- paste("EST'd ",get("eigraph.bw", env=evbase))

     x <- eiread(dbuf,"x")
     betaw <- eiread(dbuf,"betaw")
     a <- eiread(dbuf,"eaggbias") 
     lst <- listwis2(x,betaw)
     x <- as.vector(lst[[1]])
     betaw <- lst[[2]]
     plot(x,betaw, xlab=xlabel, ylab=ylabel,
          main="X by betaW", type="p", cex=psymsiz)
     pline(0,a[3,1],1,a[3,1]+a[4,1],negslope=FALSE)
     return("biasw")
   }
   if(identical(str,"tbiasb")){	###	@ x by true betaB @

    if(!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }

    psymsiz <- 0.3*get("eigraphC", env=evbase)

    xlabel <- get("eigraph.x",env=evbase)

    ylabel <- paste("True ", get("eigraph.bb",env=evbase))

    x <- eiread(dbuf,"x")
    betab <- eiread(dbuf,"truthb")
    a <- eiread(dbuf,"aggbias")
    lst <- listwis2(x,betab)
    x <- lst[[1]]
    betab <- lst[[2]]
    plot(x,betab, xlab=xlabel, ylab=ylabel,
         main="X by true betaB", type="p", cex=psymsiz)
 
    pline(0,a[1,1],1,a[1,1]+a[2,1],negslope=FALSE)
    return("tbiasb")
}  
 

  if(identical(str,"tbiasw")){###		@ x by true betaW @

    if(!(vin(dbuf,"truth"))){
      warning("ei: 'truth' must be stored first")
      return(NA)
    }

    psymsiz <- 0.3*get("eigraphC", env=evbase)
    xlabel <- get("eigraph.x", env=evbase)
    ylabel <- paste("True ",get("eigraph.bw", env=evbase))

    x <- eiread(dbuf,"x")
    betaw <- eiread(dbuf,"truthw")
    a <- eiread(dbuf,"aggbias")
    lst <- listwis2(x,betaw)
    x <- lst[[1]]
    betaw <- lst[[2]]
    plot(x,betaw, xlab=xlabel, ylab=ylabel,
         main="X by true betaW", type="p", cex=psymsiz)
    
    pline(0,a[3,1],1,a[3,1]+a[4,1],negslope=FALSE)
    return("tbiasw")
  }

   if(identical(str,"bias")){	###		@ biasB,biasW,TbiasB,TbiasW @

     op <-  par(no.readonly=TRUE)
    if(!vin(dbuf,"truth"))
      par(mfrow=c(1,2))
    else
      par(mfrow=c(2,2)) 
   
    eigraph(dbuf,"biasB",...)
    eigraph(dbuf,"biasW",...)
    if(vin(dbuf,"truth")){
      eigraph(dbuf,"TbiasB",...)
      eigraph(dbuf,"TbiasW",...)
    }else
    warning("ei: Truth is not stored")
     par(op) 
    return("bias")
  }
   if(identical(str,"boundxb")){###		@ x by bounds on betaB @

    psymsiz <- 0.3*get("eigraphC", env=evbase)
    xlabel <- get("eigraph.x",env=evbase)
    if(vin(dbuf,"truth"))
      ylabel <- paste("True ",get("eigraph.bb",env=evbase))
    else
      ylabel <- paste("Estimated ",get("eigraph.bb",env=evbase))

    x <- eiread(dbuf,"x")
    b <- eiread(dbuf,"bounds")
    minx <- ifelse((t <- min(x))>0,0,t)
    maxx <- ifelse((t <- max(x))<1,1,t)
    miny <- ifelse((t <- min(b))>0,0,t)
    maxy <- ifelse((t <- max(b))<1, 1, t)

    pline.bounds(x,b[,1],x,b[,2],negslope=FALSE,color=rainbow(length(x)),
         xylabels=c(xlabel,ylabel),title="X by bounds on betaB")
  
    if(!vin(dbuf,"truth")){
      a <- eiread(dbuf,"eaggbias")
    }else{
      a <- eiread(dbuf,"aggbias")
    }
    pline(0,a[1,1],1,a[1,1]+a[2,1],negslope=FALSE)
    return("boundxb")
    
  }

 if(identical(str,"boundxw")){###		@ x by bounds on betaW @

   psymsiz <- 0.3*get("eigraphC", env=evbase)
   xlabel <- get("eigraph.x",env=evbase)
   if(vin(dbuf,"truth"))
     ylabel <- paste("True ",get("eigraph.bw",env=evbase))
   else
     ylabel <- paste("Estimated ",get("eigraph.bw",env=evbase))

    x <- eiread(dbuf,"x")
    b <- eiread(dbuf,"bounds")
    minx <- ifelse((t <- min(x))>0,0,t)
    maxx <- ifelse((t <- max(x))<1,1,t)
    miny <- ifelse((t <- min(b))>0,0,t)
    maxy <- ifelse((t <- max(b))<1, 1, t)

    pline.bounds(x,b[,3],x,b[,4],negslope=FALSE,color=rainbow(length(x)),
         xylabels=c(xlabel,ylabel),title="X by bounds on betaW")
  
    if(!vin(dbuf,"truth")){
      a <- eiread(dbuf,"eaggbias")
    }else{
      a <- eiread(dbuf,"aggbias")
    }
   pline(0,a[3,1],1,a[3,1]+a[4,1],negslope=FALSE)
   return("boundxw") 
  }  
 if(identical(str,"boundx")){###			@ boundXB boundXW @
    op <-  par(no.readonly=TRUE)
    if(!vin(dbuf,"truth")){
      par(mfrow=c(1,2))
      eigraph(dbuf,"boundXB",...)
      eigraph(dbuf,"boundXW",...)
   
    }else{
      par(mfrow=c(2,2))
      eigraph(dbuf,"boundXB",...)
      eigraphC <- 7 
      assign("eigraphC", 7 , env=evbase)
      eigraph(dbuf,"TbiasB",...)
      eigraph(dbuf,"boundXW",...)
      assign("eigraphC", 7 , env=evbase)
      eigraph(dbuf,"TbiasW",...)
    }
    par(op)
    return("boundx")
  }
  ### plots the betab and betaw for each precinct as functions of X and N(or tvap)
  if(identical(str,"betaxn")){
    op <-  par(no.readonly=TRUE)
    betab <- eiread(dbuf,"betab")
    betaw <- eiread(dbuf,"betaw")
  
    X <- eiread(dbuf,"x")
    T <- eiread(dbuf,"t")
    N <- eiread(dbuf,"tvap")
   
    par(mfrow=c(1,2))
    eiscatterplot3d(betab,betaw,X,highlight.3d=TRUE,pch=20)
    eiscatterplot3d(betab,betaw,N,highlight.3d=TRUE,pch=20)
     par(op)
    return("betaxn")
  }
   ### plots the betab and betaw for each precinct as functions of T and N(or tvap)
  if(identical(str,"betatn")){
    op <-  par(no.readonly=TRUE)
    betab <- eiread(dbuf,"betab")
    betaw <- eiread(dbuf,"betaw")
  
    X <- eiread(dbuf,"x")
    T <- eiread(dbuf,"t")
    N <- eiread(dbuf,"tvap")
   
    par(mfrow=c(1,2))
    eiscatterplot3d(betab,betaw,T,highlight.3d=TRUE,pch=20)
    eiscatterplot3d(betab,betaw,N,highlight.3d=TRUE,pch=20)
    par(op)
 ###   scatterplot3d(betab,betaw,T,highlight.3d=TRUE,pch=20)
    return("betatn")
  }
  if(identical(str,"betast")){
    op <-  par(no.readonly=TRUE)
    betaB <- eiread(dbuf,"betab")
    betaW <- eiread(dbuf,"betaw")
    T <- eiread(dbuf,"t")
    par(mfrow=c(1,2))
    plot(betaB, T,main="Est'd betaB vs T")
    plot(betaW,T,main="Est'd betaW vs T")
    par(op)
    return("betast")
    
   }
   if(identical(str,"betasx")){
     op <-  par(no.readonly=TRUE)
     betaB <- eiread(dbuf,"betab")
     betaW <- eiread(dbuf,"betaw")
     X <- eiread(dbuf,"x")
     par(mfrow=c(1,2))
     plot(betaB, X,main="Est'd betaB vs X")
     plot(betaW,X,main="Est'd betaW vs X")
     par(op)
     return("betasx")
    
   }
   if(identical(str,"betasn")){
     op <-  par(no.readonly=TRUE)
     betaB <- eiread(dbuf,"betab")
     betaW <- eiread(dbuf,"betaw")
     N <- eiread(dbuf,"tvap")
     par(mfrow=c(1,2))
     plot(betaB, N,main="Est'd betaB vs N")
     plot(betaW,N,main="Est'd betaW vs N")
     par(op)
     return("betasn")
    
   }
   if(identical(str,"betastxn")){
   op <-  par(no.readonly=TRUE)
     betaB <- eiread(dbuf,"betab")
     betaW <- eiread(dbuf,"betaw")
     X <- eiread(dbuf,"x")
     T <- eiread(dbuf,"t")
     N <- eiread(dbuf,"tvap")
    
     par(mfrow=c(3,2))
     plot(betaB, T,main="Est'd betaB vs T")
     plot(betaW,T,main="Est'd betaW vs T")
     plot(betaB, X,main="Est'd betaB vs X" )
     plot(betaW,X,main="Est'd betaW vs X")
     plot(betaB, N,main="Est'd betaB vs N")
     plot(betaW,N,main="Est'd betaW vs N")
   par(op)
     return("betastxn")
    
   }
  message("The graphics is not supported")    

  
}
##/*
##    tomog(x,t,psiu,sel);
##**
##**  Ecological inference "tomography plot" 
##**  support proc for eigraph()
##**
##** INPUTS:
##** x = explanatory variable
##** t = outcome variable
##** psiu = 5x1 parameters of the untruncated bivariate normal, 
##**      = mean1|mean2|StanDev1|StanDev2|correlation;
##**        or 0 or missing to skip contours
##**      mean1 and mean1 can each be px1 vectors, in which case they will
##**      each be averaged before drawing ONE set of contours
##** sel = rows(x)x1 1 draw tomog line solid yellow, 0 to use dashed blue
##**                 2  draw very thick red
##**
##** OUTPUT:
##** draws lines defined by b1=(t/(1-x))-(x/(1-x))*b2, where b1,b2 are the
##**       the variables for the 2 axes and each \in[0,1]
##** also draws contour lines if psiu/=0 of the untruncated bivariate normal
##**       with parameters psiu, truncated to the unit square.  Contours are
##**       drawn so that 50% and 95% of the area within the unit square
##**       falls within the respective contours
##**
##** GLOBALS:
##** _tomogClr =  { 12, 9, 10, 11, 13, 5 }  colors for each contour.  
##** _tomogPct = {.5, .95}    percentage values at which to draw contours
##**           rows(_tomogClr) must be >= rows(_tomogPct)
##**
##** OUTPUT GLOBAL:
##** _tomogVals = values of parameter that gives 50% and 95% contours
##*/
tomog <- function(x,t,psiu,sel,choice=1,tol=0.0001,npts=100){
  evbase <- parent.frame()
  if(class(evbase)=="try-error")
    evbase <- eiset()
  ev <- evbase
  
  eigraph.bblo <- get("eigraph.bblo", env=evbase)
  eigraph.bbhi<- get("eigraph.bbhi", env=evbase)
  eigraph.bwlo<- get("eigraph.bwlo", env=evbase)
  eigraph.bwhi<- get("eigraph.bwhi", env=evbase)
  tomogClr <- get("tomogClr", env = evbase)
  eigraph.bb<- get("eigraph.bb", env=evbase)
  eigraph.bw<- get("eigraph.bw", env=evbase)
  tomogPct <- get("tomogPct", env=evbase)
  tit <-get("eigraph.tit", env=evbase)
  
  eigraph.psiu0 <- try(get("eigraph.psiu", env=evbase), silent=TRUE)
  psiu0 <- try(get("psiu", env=evbase), silent=TRUE)
  if(identical(tit, "")) 
    tit <- "Tomography lines"
     
  if(scalone(sel)) sel <- x*0+1
  plotsiz <- matrix(5,nrow=2, ncol=1)
  plctrl <- -1 
  ###  /* draw lines */
  lst <- bounds1(t,x,1,tol) 
  bnds <- lst[[1]]
  indx <- seqa(1,1,rows(t))
  j <- selif(indx,sel==1)
  b1 <- missrv(bnds[j,1],0)
  b2 <- missrv(bnds[j,2],1)
  b3 <- missrv(bnds[j,3],1)
  b4 <- missrv(bnds[j,4],0)
  lst <- drawlines(b1,b2,b3,b4,npts=npts)
  matx <- lst[["matx"]]
  maty <- lst[["maty"]]
  phcev <- 0:18
  cl <- ifelse(dim(matx)[2]>=dim(maty)[2], dim(matx)[2],dim(maty)[2])
  phcev <- rep(phcev, floor(cl/length(phcev))+1)
  matplot(matx,maty,type="l", pch=phcev, 
  xlim=c(eigraph.bblo,eigraph.bbhi), ylim=c(eigraph.bwlo,eigraph.bwhi),
      xlab= eigraph.bb, ,ylab= eigraph.bw, main=tit)
  if(any(sel>1)) {
     ch <- sel[sel>1]
     ind <- unlist(sapply(as.character(ch), FUN=grep, sel))
     mx <- max(ind)
     if(length(ind) && mx<= dim(matx)[2] && mx <= dim(maty)[2]) 
       matlines(matx[,ind], maty[,ind],col="black", lwd=4)
  }
 ###  /* draw ellipses */
   if(!scalzero(psiu) && !scalmiss(psiu)) 
      MLKcontours(psiu,tomogPct)
    assign("eigraph.tit","", env=evbase)
 
  if(class(eigraph.psiu0)!="try-error")
    assign("eigraph.psiu",eigraph.psiu0, env=evbase)
  if(class(psiu0)!="try-error")
    assign("psiu",psiu0, env=evbase)
  
}
### Helper function to tomog and eigraph to draw the MLK lines
### It assumes that plot.new has been called previously
### INPUT: psiu is a vector of 5 rows with the results of MLK estimates
###        tomogPCt how many nested contours and which values to draw

 MLKcontours <- function(psiu, tomogPct,sms=1000)
 {
    r <- rows(psiu)
    p <-(r-3)/2
    if((p*2+3)!=r){
      message("tomog: problem with psiu"); return(NULL)} 
    Bb <- colMeans(as.matrix(psiu[1:p]))
    Bw <- colMeans(as.matrix(psiu[(p+1):(2*p)])) 
    bb <- Bb
    bw <- Bw
  ###  psym <- matrix(c(Bb,Bw,tomogClr[1],5,13,1,0),ncol=1)
    sb <- psiu[r-2] 
    sw <- psiu[r-1] 
    rho <- psiu[r] 
    sb2 <- sb*sb 
    sw2 <- sw*sw 
    sbw <- rho*sb*sw 
    vc <- matrix(c(sb2,sbw,sbw,sw2),nrow=2, ncol=2, byrow=TRUE)
    tt <- matrix(seqas(0,2*3.1415927,1000),ncol=1)
    tt <- cbind(sin(tt),cos(tt))
    ch <-  try(chol(vc,pivot=FALSE), silent=TRUE)
    if(class(ch)=="try-error") 
      ch <- eichol(vc,pivot=TRUE,mess="MLKcontours..eigraph")
    v <- tt %*% ch
    
    ###   /* decide on how far out to draw elipses */
    bnds <- matrix(c(0,1,0,1),ncol=2,byrow=TRUE)
    d <- rndbtn(bb,bw,sb,sw,rho,bnds,sms) 
    d1 <- d[,1]
    d2 <- d[,2]
    z1 <- (d1-bb)/sb
    z2 <- (d2-bw)/sw
    omr <- 1-rho*rho
    dist2 <- ((z1*z1+z2*z2-2*rho*(z1*z2))/omr)
    dist2 <- sortc(dist2,1)
    
    tomogVals <- tomogPct 
    for(n in  1:length(tomogPct)){
       tomogC <- sqrt(dist2[tomogPct[n]*sms])
       tomogVals[n] <- tomogC
       x <- bb+tomogC*v[,1]
       y <- bw+tomogC*v[,2]
           
       if(length(x) <2 || length(y) <2 || length(x) != length(y)) {
         message("tomog:input error")
         return(NULL)
       }
       lst <- plines(x,y,nolag=FALSE)
      
 
    
     }
}
    
### As in the Gauss code


plines <- function(x,y,color="black",nolag=TRUE){
  xst <- trimr(lag(x),1,0)
  yst <- trimr(lag(y),1,0)
  xend <- xst <- trimr(x,1,0)
  yend <- yst <- trimr(y,1,0)
  lst <- c(list(xst), list(yst),list(xend),list(yend))
  names(lst) <- c("xst","yst","xend","yend")
 matx <- cbind(lst[["xst"]],lst[["xend"]])
  maty <- cbind(lst[["yst"]],lst[["yend"]])
  lst1 <- c(list(matx),list(maty))
  names(lst1) <- c("x", "y")
   ###    lst1 <- pline(lst[["xst"]],lst[["yst"]],lst[["xend"]],lst[["yend"]])
    
       if( length(lst1[[1]]) != length(lst1[[2]])) {
         message("PLINES:input error")
         return(NULL)
      }
 if(nolag)
       lines(lst1[[1]][,1],lst1[[2]][,1],col=color)
  else
     lines(lst1[[1]],lst1[[2]],col=color)
}
###  As in the gauss equivalent. Draws straight lines starting at (xst, yst)
###  and ending at (xend, yend)
###  proc 0=pline(xst,yst,xend,yend);
###
###  INPUT: xst or start of x-range (scalar, vector or 1 column matrix)
###         yst or start of y-range (scalar, vector or 1 column matrix)
###         xend or stop of x-range (scalar, vector or 1 column matrix)
###         yend or stop of y-range (scalar, vector or 1 column matrix)
###         negslope boolean that assumes a  negative slope
###
### OUTPUT plot multiple lines 

pline <- function(xst,yst,xend,yend,negslope=TRUE,
                  xlim=c(0,1), ylim=c(0,1),new=TRUE,color="black",npts=100){
  if(negslope){
     tmp <- xend
     xend <- xst
     xst <- tmp
   }
  num <- as.vector(yend - yst)
  den <- as.vector(xend - xst)
  vert <- horz <- FALSE
  if(all(den<=1.e-5))vert <- TRUE
  if(all(num<=1.e-5)) horz <- TRUE
  slope <- num/(den+1.e-10)
  intercp <- as.vector(yst) - slope*as.vector(xst)
  max.x <- round(max(c(xend,xst)))
  if(max.x>1) max.x <- 1
  min.x <- round(min(c(xst,xend)))
  if(min.x <0) min.x <- 0
  max.y <- round(max(c(yend,yst)))
  if(max.y>1) max.y <- 1
  min.y <- round(min(c(yst,yend)))
  if(min.y < 0) min.y <- 0
  x <- seq(from=min.x,to=max.x, length.out=10)
  y <- seq(from=min.y,to=max.y, length.out=10)
  if(!new)
  plot(x, y,type="n")
  lns <- length(slope)
  lnc <- length(color)
  if(lnc > 1 && lnc < lns)
    color <- rep(color,round(lns/lns))
  
  for(n in 1:length(slope)){
    if(length(color) >1)
      abline(a=intercp[n], b=slope[n],col=color[n])
    else
     abline(a=intercp[n], b=slope[n],col=color)
  }
  
}
###  Draws straight lines verticals or horizontals starting at (xst, yst)
###  and ending at (xend, yend) either xst=xend (vertical) or yst=yend (horizontal) 
###  
###
###  INPUT: xst or start of x-range (scalar, vector or 1 column matrix)
###         yst or start of y-range (scalar, vector or 1 column matrix)
###         xend or stop of x-range (scalar, vector or 1 column matrix)
###         yend or stop of y-range (scalar, vector or 1 column matrix)
###        
###
### OUTPUT plot multiple lines verical or horizontal

 pline.bounds <- function(xst,yst,xend,yend,negslope=TRUE,
                         xlim=c(0,1), ylim=c(0,1),new=TRUE,
                         color=NULL,npts=100,
                         xylabels=c("",""),title=""){
   num <- as.vector(yend - yst)
   den <- as.vector(xend - xst)
   vert <- horz <- FALSE
   if(all(den<=1.e-5))vert <- TRUE
   if(all(num<=1.e-5)) horz <- TRUE
   if(!vert && !horz) stop("Use the pline method instead")
   if(length(color)<=0) color <- rainbow(length(xst))
   max.x <- round(max(c(xend,xst)))
   if(max.x>1) max.x <- 1
   min.x <- round(min(c(xst,xend)))
   if(min.x <0) min.x <- 0
   max.y <- round(max(c(yend,yst)))
   if(max.y>1) max.y <- 1
   min.y <- round(min(c(yst,yend)))
   if(min.y < 0) min.y <- 0
   lns <- length(xst)
   lnc <- length(color)
   if(lnc > 1 && lnc < lns)
     color <- rep(color,round(lns/lns))
   if(vert){  
   lstx <- lapply(xst,rep,npts)
   mat <- cbind(yst,yend)
   lsty <- lapply(as.list(1:length(yst)),function(n,mat){
        seq(from=mat[n,1],to=mat[n,2], by=(abs(mat[n,1]-mat[n,2]))/(npts-1))},mat)
   matx <- matrix(unlist(lstx), nrow=length(xst), byrow=TRUE)
   maty <- matrix(unlist(lsty),  nrow=length(xst), byrow=TRUE)
   dm <- dim(matx)
   phcev <- 0:18
   ln <- length(phcev)
   cl <- dm[1]
   phcev <- rep(phcev,floor(cl/ln)+1)
   matplot(t(matx),t(maty),type="l",xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           xlab=xylabels[1], ylab=xylabels[2],main=title,pch=phcev)
  matpoints(t(matx)[c(1,npts),],t(maty)[c(1,npts),],col="black",pch=19)
   return("vertical bounds")
 }
   if(horz){  
   lsty <- lapply(yst,rep,npts)
   mat <- cbind(xst,xend)
   lstx <- lapply(as.list(1:length(xst)),function(n,mat){
        seq(from=mat[n,1],to=mat[n,2], by=(abs(mat[n,1]-mat[n,2]))/(npts-1))},mat)
   matx <- matrix(unlist(lstx), nrow=length(yst), byrow=TRUE)
   maty <- matrix(unlist(lsty),  nrow=length(yst), byrow=TRUE)
   dm <- dim(matx)
   phcev <- 0:18
   ln <- length(phcev)
   cl <- dm[1]
   phcev <- rep(phcev,floor(cl/ln)+1)
   matplot(t(matx),t(maty),type="l",xlim=c(min.x,max.x),ylim=c(min.y,max.y),
           xlab=xylabels[1], ylab=xylabels[2],main=title,pch=phcev)
   matpoints(t(matx)[c(1,npts),],t(maty)[c(1,npts),],col="black",pch=19)
  
   return("bounds horizontal")
 }
 }
### Helper to tomog plots to draw multiple lines but instead
### of displaying them in a graph finds the sequence of points
### It is assumed that the slope is negative 
### It returns two matrices, which columns are those sequences
### Using matplot it will display the multiple lines.
###
### INPUT: b1 vector or scalar with the lower boundary for x-coordinate
###        b2 vector or scalar with the upper boundary for x-coordinate
###        b3 vector or scalar with the lower boundary for y-coordinate
###        b4 vector or scalar with the upper boundary for y-coordinate
###
### OUTPUT: List with two matrices matx and maty for the matplot graphics
###
drawlines <- function(b1,b2,b3,b4,npts=150){
  data <- data.frame(b1,b2,b3,b4)
  names(data) <-c("LoB", "UpB","LoW","UpW")

### y - LoW = (UpW- LoW)/(LoB-UpB) * (x-UpB)
  den <- data[1]-data[2]+1.e-10
  num <- data[4]-data[3]
  slope <- (num / den)[[1]] ###(UpW- LoW)/(LoB-UpB)
 ### hist(slope)
 
###  print(slope)
  clst <- as.list(1:length(b1))
  xvec <- lapply(clst, function(n,b1,b2){
    return( seq(from=b1[n],to=b2[n],by=(b2[n]-b1[n])/(npts-1)))},b1,b2)

  yvec <- lapply(clst, function(n,xvec,slope,b3){
    return((xvec[[n]]-b2[n])*slope[n]+b3[n])},xvec,slope,b3)
  ln <- length(xvec[[1]])
  maty <- matrix(unlist(yvec),nrow=ln)###each column a precint
  matx <- matrix(unlist(xvec),nrow=ln)
  lst <- c(list(matx),list(maty))
  names(lst)<-c("matx", "maty")
  return(lst)
}
### Draws confidence intervals for the precints
### INPUT b is a data.frame or matrix of 4 columns
###       with the coordinates for confidence intervals
###
### OUTPUT Draw confidence intervals in red and
###        the tomography lines in blue

 conf.intv <- function(b,tit,npts=500,dbuf)
  {
    evbase <- get("evbase", env=parent.frame())
    eigraph.bblo <- get("eigraph.bblo",env=evbase) 
    eigraph.bbhi <- get("eigraph.bbhi", env=evbase)
    eigraph.bwlo <- get("eigraph.bwlo", env=evbase)
    eigraph.bwhi <- get("eigraph.bwhi", env=evbase)
    eigraph.bb <- get("eigraph.bb", env=evbase)
    eigraph.bw <- get("eigraph.bw", env=evbase)
    b10 <- missrv(b[,1],0)
    b11 <- missrv(b[,1],1)
    b21 <- missrv(b[,2],1)
    b20 <- missrv(b[,2],0)
    b31 <- missrv(b[,3],1)
    b30 <- missrv(b[,3],0)
    b40 <- missrv(b[,4],0)
    b41 <- missrv(b[,4],0)
    lst1 <- drawlines(b10,b21,b31,b40,npts=500)
    matx1 <- lst1[[1]]
    maty1 <- lst1[[2]]
#### the tomography lines
    a <- eiread(dbuf,"bounds")
    a10 <- missrv(a[,1],0)
    a11 <- missrv(a[,1],1)
    a21 <- missrv(a[,2],1)
    a20 <- missrv(a[,2],0)
    a31 <- missrv(a[,3],1)
    a30 <- missrv(a[,3],0)
    a40 <- missrv(a[,4],0)
    a41 <- missrv(a[,4],1)
    lst2 <- drawlines(a10,b11,b41,a40,npts=500)
    matx2 <- lst2[[1]]
    maty2 <- lst2[[2]]
    lst3 <- drawlines(b20,a20,a31,b30,500)
    matx3 <- lst3[[1]]
    maty3 <- lst3[[2]]
###     par(new=TRUE)
    phcev <- 0:18
    ln <- length(phcev)
    dm1 <- dim(matx1)[2]
    dm2 <- dim(maty1)[2]
    cl <- ifelse(dm1>=dm2, dm1,dm2)
    phcev <- rep(phcev, floor(cl/ln)+1)
    matplot(matx1,maty1,type="l",col="red",lwd=2,pch=phcev,
            xlim=c(eigraph.bblo,eigraph.bbhi), ylim=c(eigraph.bwlo,eigraph.bwhi),
            xlab= eigraph.bb,ylab= eigraph.bw, main=tit)
     
    matlines(matx2,maty2,type="l",col="blue",lwd=1,pch=phcev)
    matlines(matx3,maty3,type="l",col="blue",lwd=1,pch=phcev)
         
    
  }
   
##/*
##** profileit(dbuf)
##**
##** Procedure to plot log-posterior function by each parameter, 
##** holding other parameter values at their MLE's.  Parameters
##** are plotted with range eiread(dbuf,"Ebounds") (CML bounds)
##** 
##** Input: dbuf = EI data buffer
##**        R    = 0 for log-posterior, 1 for just the R (cdfbvn) function
##**
###   eigraph.pro=c(list(c(-1,-1,-1,-1,-1,0,0))*number1,list(c(1,1,1,1,1,0,0))*number2)
##** 
##** Outputs: Graphs for each Parameter.
##**
##** Global:  _eigraph_pro, replaces eiread(dbuf,"ebounds") for plotting
##**          default={.}, use ebounds.
##** 
##** Error Checking: Generates message if infinities are encountered, 
##** pushes negative (positive) infinities down (up) to -60,000 (60,000)
##**
##*/

profileit <- function(dbuf,r,npts=100,eigraph.pro=NA,...){
 ### print(eigraph.pro)
  Erho <- eiread(dbuf,"Erho")
  Esigma <- eiread(dbuf,"Esigma")
  Ebeta <- eiread(dbuf,"Ebeta")
  EalphaB <-eiread(dbuf,"EalphaB") ## depends on eimodels NA
  EalphaW <- eiread(dbuf,"EalphaW") ## depends on eimodels NA
  Eeta <- eiread(dbuf,"Eeta")
  EnumTol <- eiread(dbuf,"Enumtol")
  phi <- eiread(dbuf,"phi")
  if(all(is.na(phi))){
    message("ei: Graph 'profile' is not set for non-parametric estimations")
    return(NA)
  }
  evbase <- try(get("evbase", env=parent.frame()), silent=TRUE)
  if(class(evbase)=="try-error") evbase <- eiset()
  if(!exists("eigraph.pro")||is.na(eigraph.pro))
    eigraph.pro <- get("eigraph.pro", env=evbase)
  if(is.list(eigraph.pro))
    eigraph.pro <- matrix(unlist(eigraph.pro),nrow=length(eigraph.pro[[1]]))
  print(eigraph.pro)
  
  Ez <- eiread(dbuf,"Ez")
  Zb <- eiread(dbuf,"Zb")
  Zw <- eiread(dbuf,"Zw")
  x <- eiread(dbuf,"x")
  ebounds <- Ebounds <- eiread(dbuf,"Ebounds",calculate=TRUE)
  
  data <- eiread(dbuf,"DataSet")
  phi <- eiread(dbuf,"phi")
  vc <- eiread(dbuf,"VCphi")
  pnames<- as.vector(eiread(dbuf,"ParNames"))
  v <- as.matrix(vc[col(vc)==row(vc)])
  rr <- rows(phi)
  
  if(any(is.na(eigraph.pro))){
    if(length(Ebounds)<= 1) ebounds <- get("cml.bounds", env=evbase)
    min <- ebounds[,1]
    max <- ebounds[,2]
  }else{  
    min <- eigraph.pro[,1]
    max <- eigraph.pro[,2]

  }
  numpgs <- ceiling((rr-2)/9) ###	@number of pages of graphs@
  cnt <- 0
  for( k in (1:numpgs)){ ###numpges to graphics each has  
    message("graphing page ", k)
    rs <- 9*(k-1)+1
    re <- ifelse( (k==numpgs),rr-2,9*k)
    cl <- rs*re -3
    par(mfrow=c(cl,3))

    for( i in rs:re){ ###each subplot in the page
       
      vect <- seqas(min[i],max[i],npts)
      size <- rows(vect)
      lik <- matrix(1,nrow=size,ncol=1)
      lik <- sapply(as.list(1:size), function(j)
                    {
                      if( i==1)
                        parms <- c(vect[j],phi[(i+1):rr])
                      else
                        parms <- c(phi[1:(i-1)],vect[j],phi[(i+1):rr])
                      
                      if (r==0)
                        res <- sumc(na.omit(eiloglik(parms,data,evbase)))
                      else{
                        ## function(bb,bw,sb,sw,rho){
                        lst <- eirepar(parms,as.vector(Zb),as.vector(Zw),as.vector(x),Ez,evbase)
                        res <- sumc(na.omit(lncdfbvnu(lst[[1]],lst[[2]],lst[[3]],lst[[4]],lst[[5]])))
                      }
                      return(res)
                    } )
      lik <- unlist(lik)
      ylabel <- "Log-posterior"
      xlabel <- pnames[i]
      plot(vect,lik,xlab=xlabel,ylab=ylabel,type="p", col="blue",lwd=1.3)
      
      pline(phi[i],-1e10,phi[i],maxc(lik),negslope=FALSE)

      cnt <- cnt +1
      if(cnt%%9<= 0)  user.prompt()  
    }###end for(i in rs:re) 
   
  }###ends for(k in1:numpgs)
  
}



chkPackages <- function(pkgName){
  res <- search()
  if(length(grep("package", pkgName))<= 0)
    str <- paste("package",pkgName,sep=":")
  else
    str <- pkgName 
  b <- str %in% res
  no <- grep(str, res)
  ret <- ifelse((b || length(no) >0), TRUE, FALSE)
  return(ret)
}
