
##/*
##**  This archive is part of the program EI
##**  (C) Copyright 1995-2001 Gary King
##**  All Rights Reserved.
##*/
##/*
##   v = eiread(dbuf,"name");
##**
##** Extracts, computes, or prints results from EI or EI2 output data buffers
##**
##** When used with an EI2 output data buffer (i.e., for 2xC tables), most
##** eiread items use the mean posterior estimate for X; the
##** multiply imputed values (x2) are used only for these options and their
##** many derivatives: betaBs, betaWs, CI80bw, CI95bw, coverage, aggs.
##**
##** INPUTS:
##** dbuf     a data buffer created by ei()
##** name     a string with the name of the element to read or compute.
##**          choose from this lists below (*=option prints output to screen
##**          if _Eprt>0)
##** OUTPUT
##** v = chosen item, or missing value if item is not available
##**
##** GLOBAL
##** _Eprt  = 0 don't print; 1 = do print selected items.
##** _EIMetaR = If dbuf is a meta-data buffer (output from ei2, eimodels_run, or
##**           eimodels_def), this global denotes which of the imputed data buffers 
##**           stored in dbuf should be accessed when running this procedure (default=1)
##**
##** STORED (can also be retrieved with vread):
##** _EalphaB  value of this global (prior on alphaB)
##** _EalphaW  value of this global (prior on alphaW)
##** _Ebeta    value of this global (priors on bb,bw)
##** _Ebounds  value of this global (bounds)
##** _Ecdfbvn  value of this global (method for CDF of bivariate normal calc)
##** _EdirTol  value of this global (tolerance of CML convergence)
##** _EdoML    value of this global (do maxlik)
##** _EdoML_phi value of this global (input phi's)
##** _EdoML_vcphi value of this global  (input vc of phi's)
##** _Eeta     value of this global (slope coeff on X in betab|betaw eqns)
##** _EIgraph_bvsmth value of this global (bivariate smoothing parameter)
##** _EisChk   value of this global (check importance sampling)
##** _EiLlikS  value of this global (log-likelihood at each simulation)
##** _EisFac   value of this global (variance factor in importance sampling)
##** _EisN     value of this global (extra sims factor for importance sampling)
##** _Eist     value of this global (multivar t or normal for imptce samplng)
##** _EmaxIter value of this global (maximum number of iterations for CML)
##** _EnonEval value of this global (nonpar point evaluations)
##** _EnonNumInt value of this global (nonpar numerical integration points)
##** _EnonPar  value of this global (nonparametric estimation)
##** _EnumTol  value of this global (numerical tolerance for homogeneous pcts)
##** _Erho     value of this global (prior on rho)
##** _Eselect  value of this global (observations to select)
##** _EselRnd  value of this global (delete randomly selected obs)
##** _Esigma   value of this global (priors on sb,sw)
##** _Esims    value of this global (simulations)
##** _Estval   value of this global (starting values or grid search)
##** _EI_vc    value of this global (parameters for variance computation)
##** date      a string containing the date and time execution completed,
##**           and the EI version number and date
##** t         p x 1: outcome var proportion (turnout) (this gives V if ei2)
##** x         p x 1: explanatory var proportion (black vap)
##** n         p x 1: number of individuals per observation (precinct)
##** GhActual  value of output global _GhActual from gvc()
##** retcode   CML return code or 3333 for grid search
##** parnames  character vector of names for phi (=_cml_parnames)
##** phi       MLE's from CML
##** vcphi     global vc matrix of coeff's phi from gvc(). If _EI_vc={-1 0},
##**           it returns the inverse of vc matrix.
##** loglik    value of log-likelihood at the maximum (unnormalized)
##** PhiSims   if _EisChk==1:  _Esims x rows(phi): sims of phi; 
##**           else            meanc(phi)~stdc(phi)
##** lnir      if _EisChk==1: ln(Importance Ratio)~(normal simulations of phi') 
##**           [_Esims*_Eisn x rows(phi)+1], or scalar zero otherwise.
##** lliksims   _Esims x 1: log-likelihood, or scalar mean if _EiLlikS=0
##** resamp    scalar: number of resampling tries
##** betaBs    p x _Esims: simulations of betaB
##**
##** CALCULATED OR CHANGED:
##** Zb        matrix of covariates for betaB or 1 for none
##** Zw        matrix of covariates for betaW or 1 for none
##** meanIR    scalar: ln of mean importance ratio
##** betaWs    p x _Esims: simulations of betaW
##** RNbetaBs  p x _Esims: randomly horizontally permuted simulations of betaB
##** RNbetaWs  p x _Esims: randomly horizontally permuted simulations of betaW
##** STbetaBs  p x _Esims: SORTED simulations of betaB 
##**           (e.g., 80% CI, lower bound is STbetaBs[int(0.1*_Esims)]);
##** STbetaWs  p x _Esims: SORTED simulations of betaW
##**           (e.g., 80% CI, upper bound is STbetaWs[int(0.9*_Esims)]);
##** CI95b     p x 2: lower~upper 95% confidence intervals for betaB
##** CI95w     p x 2: lower~upper 95% confidence intervals for betaW
##** CI95bw    p x 4: lowerB~upperB~lowerW~upperW 95% conf ints for betaB betaW
##** CI80b     p x 2: lower~upper 80% confidence intervals for betaB
##** CI80w     p x 2: lower~upper 80% confidence intervals for betaW
##** CI80bw    p x 4: lowerB~upperB~lowerW~upperW 80% conf ints for betaB betaW
##** CI50b     p x 2: lower~upper 50% confidence intervals for betaB
##** CI50w     p x 2: lower~upper 50% confidence intervals for betaW
##** checkR    rows(phi)x2 precision of R (for +/- _Edirtol) is ok (1) or not (0)
##** R         scalar: sum(ln(R)), where R=volume above the unit square
##** Ri        p x 1: ln(R), where R=volume above the unit square
##** dataset   Zb~Zw~x~t, used for input to eiloglik(); _EselRnd<1 is ignored
##** loglikS   value of log-likelihood at the maximum for each i (unnormalized)
##** Ebounds   params x 2: lower~upper constraints on maxlik searching routines
##** Nb        p x 1: denominator of x and t; x.*n (blacks of voting age)
##** Nw        p x 1: (1-x).*n = n-Nb (whites of voting age)
##** Nt        p x 1: n.*t (number of people who Turnout)
##** beta      p x 2 E(betaB)~E(betaW) for each precinct
##** betaB     p x 1 E(betaB) for each precinct
##** betaW     p x 1 E(betaW) for each precinct
##** sbetaB    p x 1 stand deviation of betaB
##** sbetaW    p x 1 stand deviation of betaW
##** CsbetaB   p x 1 CI-based stand deviation of betaB
##** CsbetaW   p x 1 CI-based stand deviation of betaW
##** GEbw      p x 3 betaB~betaW~Nsims based on sims where betaB>=betaW
##** GEbwa     2 x 1 aggregate B^b ~ B^w based on sims where betaB>=betaW
##** GEwb      p x 3 betaB~betaW~Nsims based on sims where betaW>=betaB
##** GEwba     2 x 1 aggregate B^w ~ B^b based on sims where betaW>=betaB
##** bounds    p x 4: bounds on betaB & betaW, lowerB~upperB~lowerW~upperW
##** bounds2   p x 4: same as bounds but for ei2 bounds on the lambdas
##** abounds   2 x 2: aggregate bounds rows:lower,upper; cols:betab,betaw
##** abounds2  2 x 2: aggregate bounds, ei2 rows:lower,upper;cols:lambdab,lambdaw
##** Pphi      2 x 5 row 1: phi, row 2: stand errors
##** psiu      reparameterized phi into untruncated scale (average for bb bw)
##** mpPsiu    psiu from Mean Posterior rather than MLEs
##** psi       reparameter'd phi into ultimate truncated scale (avg for bb bw)
##** aggs      _Esims x 2: sims of district-level weighted avg of betaBs~betaWs
##** Maggs     2 x 1: point est of 2 dist-level parameters: meanc(aggs);
##** VCaggs    2 x 2: var-cov matrix of 2 dist-level params vcx(aggBs~aggWs)
##** Paggs     2 x 2: row 1: ec inf coeff's; row 2: standard errors
##** Goodman   2 x 2: row 1: Goodman's Regression coeff's, row 2: stan errors
##** double    2 x 1: double Regression coefficients (takes ei2 dbuf as input)
##** Thomsen   2 x 1: Thomsen's Ecological logit Estimates (bb|bw)
##** Neighbor  2 x 1: Freedman et al.'s neighborhood model estimates (bb|bw)
##** Palmquist scalar: Palmquist's Inflation factor
##** Eaggbias  4 x 2:(b~se)|(b~se) for regs of est'd betaB,betaW on const. & X
##** nobs      scalar: number of observations
##** tsims     100 x _Esims+1:rows=seqase(0,1,100), cols=X~(sorted sims of T|X)
##** tsims0    p x _Esims+1: cols=T~(sorted sims of T|X_obs,Z_obs)
##** expvarci  100 x 4: X~20%CI~mean~80%CI of sims from p(T|X)
##** expvarci0 p x 4: T~20%CI~mean~80%CI of sims from p(T|X_obs,Z_obs)
##** expvarcis 100 x 4: X~20%CI~mean~80%CI of sims from p(T|X) LOESS Smoothed
##** etaC      2x1 etaB|etaW to fix eta's at (derived from _Eeta)
##** etaS      2x1 se(etaB|etaW) fixed se's for eta (derived from _Eeta)
##** _Ez       2x1: n of covariates, including implied constant term for Zb|Zw
##**                also sets this global
##** sum       prints a summary of selected results
##**
##** OPTIONAL: if desired, must be vput these into data buffer before eiread()
##** titl      string: descriptive info about run
##** truth     p x 2: true precinct betaB~betaW
##**
##** CALCULATED FROM OPTIONAL:
##** truthB    truth[.,1] truth for blacks
##** truthW    truth[.,2] truth for whites
##** NbT       px1 number of blacks who Turnout
##** NbN       px1 number of blacks who don't Turnout
##** NwT       px1 number of whites who Turnout
##** NwN       px1 number of whites who don't Turnout
##** aggtruth  2 x 1: true aggregate coeff for blacks|whites
##** psitruth  5 x 1: true values of psi on truncated scale
##** coverage  2 x 4: %w/in CI's: 50b~80b~50w~80w (1st row=means,2nd=wtd means)
##** aggbias   4 x 2:(b~se)|(b~se) for regressions of betaB,betaW on const. & X
##** truPtile  p x 2: true percentile at which true value falls, betaB~betaW
##**
##** ADDITIONAL OPTIONS IF DATA BUFFER WAS CREATED BY EI2:
##** x2        p x _Esims: simulations of x from prior stage analysis
##** x2rn      p x _Esims: x2 randomly horizontally permuted 
##** x         p x 1: redefined as mean posterior
##** Nb2       p x _Esims: denominator of x and t; x.*n (blacks of voting age)
##** Nw2       p x _Esims: (1-x).*n = n-Nb (whites of voting age)
##** _t        p x 1: the original variable T from the first stage
##** _x        p x 1: the original variable X from the first stage
##** _n        p x 1: the original variable N from the first stage
##**
##** ADDITIONAL OPTIONS IF DATA BUFFER WAD CREATED BY EIMODELS_AVG
##** t = t used for the estimation; if models use different t's, only t 
##**     from the first model is stored and used for eiread() calculation.
##** x = x used for the estimation; if models use different x's, only x 
##**     from the first model is stored and used for eiread() calculation.
##** n = n used for the estimation; if models use different n's, only n 
##**     from the first model is stored and used for eiread() calculation.     
##** betabs = betab draws from its BMA posterior (p x _Esims)
##** betaws = betaw draws from its BMA posterior (p x _Esims)
##** mfreq  = first column contains the model number and the second 
##**          column contains the frequency of sampling from each 
##**          model (# of model x 2)
##** postprob = first column contains the model number and the second 
##**            column contains the posterior model probabilities 
##**            (# of model x 2)
##** prprob = first column contains the model number and the second 
##**          column contains the prior model probability used for 
##**          calculation (# of model x 2).
##** margllik = first column contains the model number and the second 
##**            column contains the marginal log-likelihood for each 
##**            model (# of model x 2)
##**
##*/

eiread <- function(dbuf, str, formula=NA,calculate=FALSE,...){
  require(mvtnorm)
  str0 <- str
  str <- tolower(str)
  nm <- names(dbuf)
  nm <- sapply(nm, tolower)
  names(dbuf) <- nm
  if(vin(dbuf,str)&&!calculate)
    return(vread(dbuf,str))
 
###setting global variables
  evbase <- dbuf[["evbase"]]
  if(is.null(evbase)){
    evbase <- try(get("evbase", env=parent.frame()), silent=TRUE)
  if(class(evbase) == "try-error"|| length(evbase)<=0)
    evbase <- eiset(t=NULL,x=NULL,tvap=NULL,Zb=1,Zw=1,...)  ##environment with default parameters
  }
  drvdot <- match.call(expand.dots=TRUE)
  drv  <-  match.call(expand.dots=FALSE)
  if(exists("n"))
    tvap <- n
  else if(exists("tvap"))  
    n <- tvap
###extra parameters supplied with the function call:...
  entries <- expanddots(drvdot,drv,evbase)
 
### copy variables from evbase to local environment evei
  evei <- getEnvVar(evbase, environment())  ##environment
### R CMD check needs to know the explicit values of global variables
  if(exists("tvap")) tvap <- n <- get("tvap", env=evei)
  else if (exists("n")) tvap <- n <- get("n", env=evei)
    if(exists("Eptr"))  Eptr <- get("Eptr", env=evei)
    if(exists("Esims"))  Esims <- get("Esims", env=evei)
    if(exists("ei.bma.est"))  ei.bma.est <- get("ei.bma.est", env=evei)
    if(exists("aggs"))  aggs <- get("aggs", env=evei)
    if(exists("eagbias"))  eagbias <- get("eagbias", env=evei)
    if(exists("aggtruth"))  aggtruth <- get("aggtruth", env=evei)
    
### end of setting
  nmdbuf <- names(dbuf)
  nmdbuf <- sapply(nmdbuf, tolower)
  names(dbuf) <- nmdbuf
  if(length(entries)){
    for(x in entries){
      ln <- length(dbuf)
      nm <- names(dbuf)
      if(tolower(x) %in% nmdbuf)
        next
      dbuf[[ln +1]] <-  eval(as.symbol(x))
      names(dbuf) <- c(nm, x)
    }
  }  
      
  nmdbuf <- names(dbuf)
  nmdbuf <- sapply(nmdbuf,tolower)
  names(dbuf) <- nmdbuf
  #message("Variables in list are")
  #prettynm <- paste(names(dbuf), collapse="\t")
  #message(prettynm)
  
  if (vin(dbuf,"eimetar"))
    eimetar <- vread(dbuf,"eimetar")
  else
    eimetar <- NA
 
  if(vin(dbuf, "titl")){
    titl <- vread(dbuf, "titl")
    titl <- trim.blanks(titl)
    nc  <- nchar(as.character(floor(as.numeric(eimetar))))
 
    fmt <- formatC(as.numeric(eimetar),width=nc,digits=0, format="f")
 
    if (identical(titl, "*MDB* Meta-Data Buffer from 2nd Stage *MDB*")){ 
      dbuf <- vread(dbuf,paste("dbuf",fmt,sep=""))
    }else if(identical(titl, "*MDB* Meta-Data Buffer from eimodels_def() *MDB*")){
      if(!vin(dbuf,paste("mod.d",fmt, sep=""))){
        message("eiread: Model " + eimetar + " is not stored in this data buffer.") 
        return(NA)
      }
      dbuf <- vread(dbuf,paste("mod.d",fmt,sep=""))
      if(Eptr==3)
        message("Reading ", str, " from Model ", eimetar, "...") 
       
    }else if(identical(titl,"*MDB* Meta-Data Buffer from eimodels_run() *MDB*")){
      if(!vin(dbuf, paste("mod.r"+fmt,sep=""))){
        mess <- paste("eiread: Model", eimetar, "is not stored in this data buffer")
        message(mess)
        return(NA)
      }
      dbuf <- vread(dbuf, paste("mod.r",fmt, sep=""))
      if(Eptr == 3)
        message("Reading ", str, " from Model ", eimetar,"...")
       
    }
    
  }### if(vin(dbuf, "titl")
  strtmp <- str
  str <- tolower(str)
  str <- trim.blanks(str)
###  cv <- vnamecv(dbuf)
  cv <- names(dbuf)

  cv <- unlist(sapply(cv,tolower))
  names(dbuf) <- cv
  names(dbuf) <- sapply(names(dbuf),trim.blanks)
  strc <- str
###    /* for the output of eimodels_avg, prevent users from reading inappropriate parameters */ 
###    /* almost all of the global variables and some parameters */
  if(vin(dbuf, "titl")){
    
    titl <- vread(dbuf, "titl")
    if(identical(titl, "*DB* Data Buffer from eimodels_avg() *DB*")){
      vrs <- c("ealphab","ealphaw","ebeta","ebounds","ecdfbvn","edirtol","edoml", 
               "edoml.phi","edoml.vcphi","eeta","ei.vc","eigraph.bvsmth","eischk","eisfac",
               "eisn","eist","emaxiter","enoneval","enonnumint","enonpar","enumtol","erho",
               "eselrnd", "esigma","estval","evtol","ei2.m","eimetar","eimodels.save","zb", "zw",
               "t","x","ez"
               ###aqui unknown
               ,"x2","x2rn", "checkr", "dataset", "etac", "etas", "expvarci","expvarcis","lnir"
               ,"loglik", "logliks", "meanir", "mppsiu", "parnames","phi","phisims", "pphi"
               ,"psi", "psitruth", "psiu", "r",  "ri", "retcode", "tsims","vcphi") 
      vrs <- as.matrix(tolower(vrs))
     if(identical(str, "esims")){
     
       res <- ncol(as.matrix(eiread(dbuf,"betabs",formula,calculate=FALSE,...)))
       return(res)
     }
     if((str %inG% vrs) && !identical(str, "expvarci0")){
       mess <- paste("eiread: You cannot eiread", strtmp, "from the output buffer of eimodels.avg")
       message(mess)
       return(NA)
     }
    }
  }### if(vin(dbuf, "titl"))
  
###  /* for the output of eimodels_avg, prevent users from reading inappropriate parameters */ 
###  /* almost all of the global variables and some parameters */  
     
  if (vin(dbuf,"titl")){
    titl <- vread(dbuf, "titl")
    if(identical(titl,"*MDB* Meta-Data Buffer from eimodels_def() *MDB*")){
      vrs <- c("under.t", "or",  "under.x", "under.ez", "abounds", "abounds2",
               "aggbias", "beta", "betab", "betabs",
               "betaw", "betaws", "bounds", "checkr", "ci50b", "ci50w", "ci80b", "ci80bw", "ci95b", "ci95bw",
               "coverage", "csbetab", "csbetaw", "eaggbias", "etac", "etas", "expvarci", "expvarcis", "gebw",
               "gebwa", "gewb", "gewba", "goodman", "lnir", "loglik", "logliks", "maggs", "meanir", "mppsiu",
               "neighb", "paggs", "palmquist", "parnames", "phi", "phisims", "pphi", "psi", "psiu", "r", "ri",
               "resamp", "retcode", "rnbetabs", "rnbetaws", "sbetab", "sbetaw", "stbetabs", "stbetaws", "sum",
               "thomsen", "tsims", "vcaggs", "vcphi") 
      
      vrs <- as.matrix(tolower(vrs))
      if(str %inG% vrs){
        mess <- paste("eiread: You cannot eiread", strtmp, "from the output buffer of eimodels.def")
        message(mess)
        return(NA)
      }
    }
  }  ### if (vin(dbuf,"titl"))
  
### changes in stored globals
  res <- as.matrix(NA)
  if(vin(dbuf, str)){
    res <- vread(dbuf, str) ##returns element or NA
    res <- as.matrix(res)
  }
 
  if(identical(str, "eeta")){ ###	@ _Eeta @
    if(!vin(dbuf,"eeta")) res <- matrix(0,nrow=4,ncol=1)
    if(scalmiss(res))
      res <- matrix(0, nrow=4, ncol=1)
    else if(ncol(as.matrix(res)) == 2)
      res <- matrix(as.vector(res), ncol=1)
    else if(rows(res)==2)
      res <- rbind(res, matrix(0, nrow=2, ncol=1))
  }else if(identical(str, "zb")){ ###@ Zb supplemented with _Eeta @

    e <- eiread(dbuf, "eeta",formula,calculate=FALSE,...)
    if(any(e %in% c(1,3,4)))
      res <- vread(dbuf, "x")
         
  }else if(identical(str, "zw")){
    e <- eiread(dbuf, "eeta",formula,calculate=FALSE,...)
 
    if(any(e %in% c(5,2,3)))
      res <- vread(dbuf,"x")
         
  }else if(identical(str, "titl") || identical(str, "undertitle") && is.na(res)){
    if(!vin(dbuf,"titl")) res <- ""
    
### /***** computed results *****/
  }else if(identical(str, "x2")){
    if(!vin(dbuf, "x2"))
      message("eiread: 'x2' option is only available in data buffers created by ei2")
   
###    	@ horizontally randomly permuted x2  @
  }else if(identical(str, "x2rn")){
    
    if(vin(dbuf, "x2")){
    
      res <- as.matrix(vread(dbuf, "x2"))
      a <- rows(res)
      c <- cols(res)
      for(n in 1:a){
        if(exists("mock_runif"))
            res[n, ] <- res[n, order(mock_runif(c))]
        else
            res[n, ] <- res[n, order(runif(c))]
      }
    }else
      message("eiread: 'x2' option is only available in data buffers created by ei2")
      
  ###else if identical(str, "emaxiter"): already taken care of  
    
  }else if(identical(str,"eigraph.bvsmth") || identical(str,"bvsmth") ){
    
    if(vin(dbuf, "bvsmth"))
      res <- vread(dbuf, "bvsmth") 
  }else if(identical(str,"eimodels.save")){
    
    if(vin(dbuf, "eimsave"))
      res <- vread(dbuf, "eimsave")
  }else if(identical(str,"ei.bma.est")){
    
    if(vin(dbuf, "ei.bma.est")|| identical(dbuf, "eibmaest"))
      res <- vread(dbuf, "eibmaest")
  }else if(identical(str,"ei.bma.prior")|| identical(str,"eibmaprior")){
    
    if(vin(dbuf, "prprob")){
      priorp <- vread(dbuf,"prprob")
      res <- priorp[, 2]
    }
  
  }else if(identical(str,"doml.phi")||identical(str,"underedoml.phi")){
    if(vin(dbuf, "doml.phi"))
      res <- vread(dbuf, "doml.phi")
    
  }else if(identical(str,"enoneval")||identical(str,"enoneva")||
           identical(str,"underenoneval")||identical(str,"underenoneva")){
       
    if(vin(dbuf, "enoneva"))
      res <- vread(dbuf, "enoneva")
    else if( vin(dbuf,"underenoneva"))
      res <- vread(dbuf, "underenoneva")
    
  }else if( identical(str,"enonnumint") || identical(str,"enonnum")){
    if (vin(dbuf,"enonnum"))
      res <- vread(dbuf,"enonnum")
 

  }else if (identical(str, "meanir") && is.na(res)){
   
    if(vin(dbuf, "emeanIR") || vin(dbuf, "emeanir"))
      res <- vread(dbuf, "emeanir")
    else if(vin(dbuf, "lnir") &&
            length(na.omit(eiread(dbuf, "eischk", formula,calculate=FALSE,...)))
            && (eiread(dbuf,"eischk",formula,calculate=FALSE,...) == 0))
      res <- vread(dbuf, "lnir")
    else if(vin(dbuf, "lnir") && length(na.omit(eiread(dbuf, "eischk", formula,calculate=FALSE,...))) && eiread(dbuf,"eischk",formula,calculate=FALSE,...) == 1){
      a <- eiread(dbuf, "lnir",formula,calculate=FALSE,...)
      a <-  a[,1]
      max <- max(a)
      res <- max + log(colMeans(as.matrix(exp(a-max)))) ###  @ = ln(meanc(exp(lnir)) @
      
    }
         
  }else if (identical(str, "logliks")){ ##@ loglik for each i @

    a <- eiread(dbuf, "phi",formula,calculate=FALSE,...)
    b <- eiread(dbuf, "dataset",formula,calculate=FALSE,...)
    if(exists("evbase"))
      res <- eiloglik(a,b, evbase)
    else
      res <-  eiloglik(a,b,eiset())
  
  }else if(identical(tolower(str), "resamp")){
    res <- NA
    if(vin(dbuf, "resamp") )
      res <- vread(dbuf, "resamp")
    if(vin(dbuf, "eresamp"))
      res <- vread(dbuf, "eresamp")
       
  }else if(identical(tolower(str), "enonpar")){
    if(!vin(dbuf,"enonpar"))
      res <- 0
    else
      res <- vread(dbuf, "enonpar")
                
  }else if (identical(tolower(str),"etac")){ ###	@ 2x1 vect to fix coeff's at @
    res <- eiread(dbuf, "eeta",formula,calculate=FALSE,...)
    if(all(is.na(res))) return(res)  
    if(rows(res) == 1)
      res <- matrix(0, nrow=2, ncol=1)
    else if(rows(res) == 3 && res[1]==4)
      res <- as.matrix(c(0, res[2]))
    else if(rows(res) ==3 && res[1] == 5)
      res <- as.matrix(c(res[2], 0))
    else if(rows(res)==4)
      res <- as.matrix(res[1:2])
         
                       
  }else if(identical(tolower(str),"etas")){ ###	@ 2x1 vect to fix se's of eta at @
    res <- eiread(dbuf, "eeta",formula,calculate=FALSE,...)
    if(all(is.na(res))) return(res)  
    if(rows(res) == 1)
      res <- matrix(0, nrow=2, ncol=1)
    else if (rows(res) == 3 && res[1]==4)
      res <- as.matrix(c(0, res[3]))
    else if (rows(res) == 3 && res[1]==5)
      res <- as.matrix(c(res[3], 0))
    else if(rows(res) == 4)
      res <- as.matrix(res[3:4])
    
###  @ n of covariates, incl. implied constant for Zb|Zw @
  }else if(identical(tolower(str), "ez")){
  
    zb <- Zb <- eiread(dbuf, "zb",formula,calculate=FALSE,...)
    zw <- Zw <- eiread(dbuf,"zw", formula,calculate=FALSE,...)
    Ez <- as.matrix(c(cols(Zb) + 1 - as.numeric(Zb == 1), cols(Zw) + 1 - as.numeric(Zw == 1)))
    if(exists("evbase"))
      assign("Ez", Ez, env=evbase)
    res <- Ez
 ###@ _cml_bounds (repeated code from eicml.src @
  }else if(identical(tolower(str),"ebounds")){ ###  @ _cml_bounds (repeated code from eicml.src @
    dummy <- eiread(dbuf, "ez",formula,calculate=FALSE,...)
    b <- eiread(dbuf,"ebounds",formula,calculate=FALSE,...)
       
    if(scalzero(b)) ####@ don't use bounds  @
      res <- matrix(c(-20, 20), nrow=1, ncol=2)
    else if(cols(b) == 2)
      res <- b
    else if (scalone(b)){ ### @ automatic bounds calculation @
      e <- matrix(c(-10, 10), nrow=1, ncol=2)
      f <-  matrix(c(-20, 20), nrow=1, ncol=2) ## @{ -1e256 1e256 };@
      if(dummy[1] == 1)
        res <- e
      else
        res <- f %dot*% matrix(1, nrow=dummy[1], ncol=1)
      if(dummy[2]==1)
        res <- rbind(res, e)
      else
        res <- rbind(res, f%dot*% matrix(1,nrow=Ez[2], ncol=1))
      mat1 <-  matrix(c(-6, 3), nrow=1, ncol=2)
      mat2 <-  matrix(c(-2, 2), nrow=1, ncol=2)
      res <- rbind(res,mat1, mat1,mat2)
    }else
    message("eiread: problem with Ebounds")
  }else if(identical(tolower(str), "nobs")){ ##	@ number of observations  @
    res <- rows(vread(dbuf, "t"))
    if(Eprt > 0)
      message( "number of observations: ", res)
  }else if(identical(tolower(str), "tvap") ||identical(tolower(str), "n") ){ ##	@ total vap, old and new notation@
    strc <- str
    if(strc %inG% cv)
      res <- vread(dbuf, "tvap")
    else
      res <- vread(dbuf, "n")

  }else if(identical(tolower(str), "bvap") || identical(tolower(str), "nb")){ ###	@ black vap  @{
    x <- vread(dbuf, "x")
    n <- eiread(dbuf, "n",formula,calculate=FALSE,...)
    res <- x%dot*% n
  }else if(identical(tolower(str), "nb2")){
###	@ number of blacks turning out to vote: sims from ei2 @
    x <- eiread(dbuf, "x2",formula,calculate=FALSE,...)
    n <- eiread(dbuf, "n",formula,calculate=FALSE,...)
    res <- x%dot*% n
  }else if( identical(tolower(str),"wvap") || identical(tolower(str),"nw")) {###@ white vap  @
    x <- vread(dbuf, "x")
    n <- eiread(dbuf, "n",formula,calculate=FALSE,...)
    res <- (1-x)%dot*% n
  } else if(identical(tolower(str),"nw2")){	###@ white vap  @
    x <- eiread(dbuf,"x2",formula,calculate=FALSE,...)
    n <- eiread(dbuf,"n",formula,calculate=FALSE,...)
    res <- (1-x)%dot*%n

  } else if(identical(tolower(str),"nt")){		###	@ number of people who Turnout @
    t <- eiread(dbuf,"t",formula,calculate=FALSE,...)
    n <- eiread(dbuf,"n",formula,calculate=FALSE,...)
    res <- t%dot*%n

  } else if(identical(str,"dataset")){	###	@ dataset for input to eiloglik() @
    x <- vread(dbuf,"x")                  ### @ ignore _EselRnd @
    t <- vread(dbuf,"t")
    Zb <- eiread(dbuf,"zb",formula,calculate=FALSE,...)
    Zw <- eiread(dbuf,"zw",formula,calculate=FALSE,...)
    Ez <- eiread(dbuf,"ez",formula,calculate=FALSE,...)
    Eselect <- eiread(dbuf,"eselect",formula,calculate=FALSE,...)
    if(length(evbase))
      assign("Eselect",Eselect, env=evbase) 
    res <- packdta(x,Zb,Zw,t,evbase=evbase,Ez,Eselect) 
  }else if (identical(tolower(str),"betaws")){	###	@ betaW simulations @
    t <- vread(dbuf,"t") 
    if (vin(dbuf,"x2"))
      x <- vread(dbuf,"x2")
    else
      x <- vread(dbuf,"x")
    
    betaBs <- betabs <- eiread(dbuf,"betabs",formula,calculate=FALSE,...)
    res <- NA
    if(!scalmiss(betabs))
      res <- betab2w(t,x,betaBs) 
  } else if(identical(tolower(str),"beta")){ ###   @ E(betab)~E(betaw) by precinct @
    res <- cbind(eiread(dbuf,"betab",formula,calculate=FALSE,...),eiread(dbuf,"betaw", formula,calculate=FALSE,...))
  } else if(identical(tolower(str),"betab")){	###		@ E(betab) for each precinct @
    res <- colMeans(t(eiread(dbuf,"betabs",formula,calculate=FALSE,...)))
  } else if(identical(tolower(str),"betaw")){ ###			@ E(betaW) for each precinct @
    
    res <- colMeans(t(eiread(dbuf,"betaws",formula,calculate=FALSE,...)))
  }  else if (identical(tolower(str), "sbetab")){###		@ sd(betab) for each precinct @
    res <- sd(as.data.frame(t(eiread(dbuf,"betabs",formula,calculate=FALSE,...))))
  } else if( identical(tolower(str), "sbetaw")){ ###		@ sd(betaW) for each precinct @
    res <- sd(as.data.frame(t(eiread(dbuf,"betaws",formula,calculate=FALSE,...))))
  } else if(identical(tolower(str), "rnbetabs")){###		@ randomly permuted betabs sims @
    res <- eiread(dbuf,"betabs",formula,calculate=FALSE,...)
    a <- rows(res)
    c <- cols(res)
    for( i in 1:a){
      if(exists("mock_runif"))
          res[i,] <- res[i,order(mock_runif(c))]   
      else
          res[i,] <- res[i,order(runif(c))]   
    }
          
   
  } else if(identical(tolower(str), "rnbetaws")){ ###		@ randomly permuted betabs sims @
    res <- eiread(dbuf,"betaws",formula,calculate=FALSE,...)
    a <- rows(res)
    c <- cols(res)
    for( i in 1:a){
      if(exists("mock_runif"))
          res[i,] <- res[i,order(mock_runif(c))]
      else
          res[i,] <- res[i,order(runif(c))]
    }
  } else if( identical(tolower(str), "stbetabs")){ ###		@ sorted betaB simulations  @
    betabs <- betaBs <- eiread(dbuf,"betaBs",formula,calculate=FALSE,...) 
    res <- NA

    if(!scalmiss(betabs)){
      a <- rows(betabs)
      for(n in 1:a)
        if(all(!is.na(betabs[n,])))
          betabs[n,] <- sort(betabs[n,])
      betaBs <- betabs
###same as  betabs[n,] <- betaBs[n,] <- sortc(as.matrix(betabs[n,]))
      res <- betabs
    }
  } else if(identical(tolower(str), "stbetaws")){ ###		@ sorted betaW simulations  @
    betaWs <- betaws <- eiread(dbuf,"betaws",formula,calculate=FALSE,...) 
    res <- NA

    if(!scalmiss(betaws)){
      for(n in 1:rows(betaws))
        if(all(!is.na(betaws[n,])))
          betaws[n,] <- sort(betaws[n,])
      res <- betaWs <- betaws 
    }
    
  }else if(identical(tolower(str),"truptile")){###	@ percentile at which true value falls @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first")
      return(NA)
    }
          
    stbetabs <- eiread(dbuf,"stbetabs",formula,calculate=FALSE,...)
    stbetaws <- eiread(dbuf,"stbetaws",formula,calculate=FALSE,...)
    truth <- eiread(dbuf,"truth",formula,calculate=FALSE,...)
    b <- cols(stbetabs)
       
    res <- minindc(t(abs(stbetabs-truth[,1])))/b     
    res <- cbind(res, minindc(t(abs(stbetaws-truth[,2])))/b) 
    EnumTol <- 1.e-4
    if(length(evbase)){
      tol <- try(get("EnumTol", env=evbase),silent=TRUE)
      if(!inherits(tol,"try-error"))
        EnumTol <- as.vector(tol)
    }
       
    res[,1] <- recode(res[,1],stdc(t(stbetabs)) <= EnumTol,0.5)  ###@ homog prects @
    res[,2] <- recode(res[,2],stdc(t(stbetaws))<= EnumTol,0.5) 
          
  }else if(identical(tolower(str),"ci50b")){### @ 50% confidence intervals for betab @
    stbetabs <- eiread(dbuf,"stbetabs",formula,calculate=FALSE,...) 
    e <- cols(stbetabs) 
    res <- cbind(stbetabs[,floor(0.25*e)], stbetabs[,floor(0.75*e)])
  }else if(identical(tolower(str),"ci80b"))	{ ###@ 80% confidence intervals for betab @
    stbetabs <- eiread(dbuf,"stbetabs",formula,calculate=FALSE,...) 
    e <- cols(stbetabs) 
    res <- cbind(stbetabs[,floor(0.1*e)], stbetabs[,floor(0.9*e)]) 
    
  } else if(identical(tolower(str), "ci95b")){ ###@ 95% confidence intervals for betab @
    stbetabs <- eiread(dbuf,"stbetabs",formula,calculate=FALSE,...) 
    e <- cols(stbetabs) 
    res <- cbind(stbetabs[,floor(0.05*e)],stbetabs[,floor(0.95*e)]) 
  } else if(identical(tolower(str),"ci50w")){###	@ 50% confidence intervals for betaw @
    stbetaws <- eiread(dbuf,"stbetaws",formula,calculate=FALSE,...) 
    e <- cols(stbetaws) 
    res <- cbind(stbetaws[,floor(0.25*e)], stbetaws[,floor(0.75*e)]) 
  } else if(identical(tolower(str),"ci80w")){###		@ 80% confidence intervals for betaw @
    stbetaws <- eiread(dbuf,"stbetaws",formula,calculate=FALSE,...) 
    e <- cols(stbetaws) 
    res <- cbind(stbetaws[,floor(0.1*e)],stbetaws[,floor(0.9*e)]) 
  } else if(identical(tolower(str),"ci95w")){ ###	@ 95% confidence intervals for betaw @
    stbetaws <- eiread(dbuf,"stbetaws",formula,calculate=FALSE,...)
    e <- cols(stbetaws) 
    res <- cbind(stbetaws[,floor(0.05*e)],stbetaws[,floor(0.95*e)]) 
  }  else if(identical(tolower(str),"ci80bw")){ ###		@ 80% conf intervals for betab betaw  @
    res <- cbind(eiread(dbuf,"ci80b",formula,calculate=FALSE,...), eiread(dbuf,"ci80w",formula,calculate=FALSE,...)) 
  } else if(identical(tolower(str),"ci95bw")){###		@ 95% conf intervals for betab betaw  @
    res <- cbind(eiread(dbuf,"ci95b",formula,calculate=FALSE,...),eiread(dbuf,"ci95w",formula,calculate=FALSE,...)) 

  } else if(identical(tolower(str),"coverage")){ ###		@ CI coverage @
    truth <- eiread(dbuf,"truth",formula,calculate=FALSE,...)
    res <- NA
    if(!scalmiss(truth)){
      if(vin(dbuf,"x2")){
        nb <- eiread(dbuf,"nb2",formula,calculate=FALSE,...) 
        nw <- eiread(dbuf,"nw2",formula,calculate=FALSE,...) 
      }else{
        nb <- eiread(dbuf,"nb",formula,calculate=FALSE,...) 
        nw <- eiread(dbuf,"nw",formula,calculate=FALSE,...)  
      }
      a <- eiread(dbuf,"ci50b",formula,calculate=FALSE,...)   
      bb <- is.na(colSums(t(cbind(a,truth)))) 
      a <- subset(a,subset=!bb) 
      e <- subset(truth,subset=!bb) 
      f <- as.numeric((e[,1]>=a[,1]) & (e[,1] <=a[,2])) 
      res <- rbind(colMeans(as.matrix(f)),t(colMeans(as.matrix(meanwc(f,delif(nb,bb)))))) 
      a <- eiread(dbuf,"ci80b",formula,calculate=FALSE,...)   
      bb <- is.na(colSums(t(cbind(a,truth)))) 
      a <- delif(a,bb) 
      e <- delif(truth,bb) 
      f <- as.numeric((e[,1]>=a[,1]) & (e[,1]<=a[,2])) 
      res <- cbind(res, rbind(colMeans(as.matrix(f)), t(colMeans(as.matrix(meanwc(f,delif(nb,bb)))))))  
      a <- eiread(dbuf,"ci50w",formula,calculate=FALSE,...)  
      bb <- is.na(colSums(t(cbind(a,truth)))) 
      a <- delif(a,bb) 
      e <- delif(truth,bb) 
      f <- as.numeric((e[,2]>=a[,1]) & (e[,2] <=a[,2]))
      res <- cbind(res, rbind(colMeans(as.matrix(f)), t(colMeans(as.matrix(meanwc(f,delif(nw,bb))))))) 
      a <- eiread(dbuf,"ci80w",formula,calculate=FALSE,...)   
      bb <- is.na(colSums(t(cbind(a,truth)))) 
      a <- delif(a,bb) 
      e <- delif(truth,bb) 
      f <- as.numeric((e[,2]>=a[,1]) &(e[,2] <=a[,2])) 
      res <- cbind(res, rbind(colMeans(as.matrix(as.numeric(f))), t(colMeans(as.matrix(meanwc(as.numeric(f),delif(nw,bb))))))) 
      if(Eprt>0){
        message("CI coverage; % true values within each confidence interval") 
        vrs <- c("      ", "50%Black", "80%Black", "50%White", "80%White") 
        vrs <- as.matrix(vrs)
        print(t(vrs))
        message("    %: ")
        print(res[1,]) 
        message("Wtd %: ")
        print(res[2,]) 
      }
    }
  }  else if(identical(tolower(str),"checkr")){###    @ check R function precision @
      Edirtol <- eiread(dbuf,"Edirtol",formula,calculate=FALSE,...)  
      if(class(try(Ez,silent=T))!="try-error")
        res <- checkr(dbuf,as.vector(Edirtol),Ez)
      else
         res <- checkr(dbuf,as.vector(Edirtol),Ez=matrix(1,nrow=2, ncol=1))
    } else if(identical(tolower(str),"ri")) { ###                    @ ln(R) @    
      a <- eiread(dbuf,"phi",formula,calculate=FALSE,...)  
       if(class(try(Ez,silent=T))!="try-error")
         lst <- pluckdta(eiread(dbuf,"dataset",formula,calculate=FALSE,...),Ez)
      else{
        Ez <- matrix(1,nrow=2, ncol=1)
        lst <- pluckdta(eiread(dbuf,"dataset",formula,calculate=FALSE,...),Ez)
      }
      Zb <- zb <- lst$Zb
      Zw <- zw <- lst$Zw
      x <- lst$x
      t <- lst$y
      lst <- eirepar(a,zb,zw,x, Ez)
      bb <- lst[[1]]
      bw <- lst[[2]]
      sb <- lst[[3]]
      sw <- lst[[4]]
      rho <- lst[[5]]
   
      res <- na.omit(lncdfbvnu(bb,bw,sb,sw,rho)) 
    } else if(identical(tolower(str),"r")){ ###                    @ sum(ln(R)) @    
      res <- colSums(eiread(dbuf,"ri",formula,calculate=FALSE,...)  ) ###EV I have changed it from colSums(eiread(dbuf,"r",formula,calculate=FALSE,...)  ) 
    } else if(identical(tolower(str),"aggbias")){ ###		@ aggregation bias regressions @
      truth <- eiread(dbuf,"truth",formula,calculate=FALSE,...)  
    
      res <- NA
    if(!scalmiss(truth)){
      x <- vread(dbuf,"x")
      assign("Routput",0, env=evbase)
      assign("Rconst",1,env=evbase)
      Rconst <- 1
      Routput <- 0
      lstreg <- reg(formula,x,truth[,1],...)
      b <- lstreg$coefficients
      bb <- summary(lstreg)$coefficients[,2]
      res <- cbind(b,bb)
      lstreg <- reg(formula,x,truth[,2],...)
      b <- lstreg$coefficients
      bb <- summary(lstreg)$coefficients[,2]
      res <- rbind(res,cbind(b,bb))
      Eprt <- -1
      if(Eprt>0){
        vrs <- (c("TRUEDepV", "       ", "coeffs", "se's")) 
        
        print(vrs)
        const <- 
      ##  vrs <- as.matrix(c(const, x, const, x)) 
        a <- cbind(vrs,res) 
        vrs <- as.matrix(c("betaB", "   ", "betaW", "   ")) 
        a <- cbind(vrs,a)
        print(a)
        ###for gauss display of pretty formatting
        b <- matrix(1, nrow=4,ncol=1) 
        mask <- cbind(matrix(0,nrow=4,ncol=2), matrix(1,nrow=4,ncol=2)) 
	
      fmt <- matrix(c("-*.*s ", 8, 8,
                     "-*.*s ", 8, 8,
                     "*.*lf", 7, 4,
                     "*.*lf", 7, 4), nrow=4, ncol=3)
 ###       call printfm(a,mask,fmt);		  
      }
    }
    }else if(identical(tolower(str),"eaggbias")){###	@ estimated aggregation bias regressions @
      betab <- eiread(dbuf,"betab",formula,calculate=FALSE,...)   
      betaw <- eiread(dbuf,"betaw",formula,calculate=FALSE,...)  
      res <- NA
    if(!scalmiss(betaw)){
      x <- eiread(dbuf,"x",formula,calculate=FALSE,...)  
      assign("Routput",0, env=evbase)
      assign("Rconst",1, env=evbase)
      Rconst <- 1 
      Routput <- 0
      formula <- get("formula",env=environment())
      if(is.function(formula)) formula <- NA
      lstreg <- reg(formula,x,betab,...)
      b <- lstreg$coefficients
      bb <- summary(lstreg)$coefficients[,2]
      res <- cbind(b,bb)
      lstreg <- reg(formula,x,betaw,...)
      b <- lstreg$coefficients
      bb <- summary(lstreg)$coefficients[,2]
      res <- rbind(res, cbind(b,bb));
      Eprt <- -1
      if(Eprt>0){
        vrs <- c("ESTDepV", "       ", "coeffs", "se's")
        print(vrs)
        vrs <- as.matrix(c(const, x, const, x))
        a <- cbind(vrs,res)
        vrs <- as.matrix(c("betaB", "   ", "betaW", "   ")) 
        a <- cbind(vrs,a)
        print(a)
        ###gauss formatting 
        b <- matrix(1,nrow=4,ncol=1) 
        mask <- cbind(matrix(0,nrow=4,ncol=2), matrix(1,nrow=4,ncol=2)) 
        fmt <- matrix(c("-*.*s ", 8, 8,
                        "-*.*s ", 8, 8,
                        "*.*lf", 7, 4,
                        "*.*lf", 7, 4), nrow=4, ncol=3)
   ##     call printfm(a,mask,fmt);		   
      }
    }
    } else if(identical(tolower(str),"csbetab")){###		@ CI-based sd(betaB)  @
      stbetabs <- stbetaBs <- eiread(dbuf,"stbetabs",formula,calculate=FALSE,...)   
      a <- stbetaBs[,floor(cols(stbetabs)*0.3413)] ### @ 34th percentile @
      b <- stbetaBs[,floor(cols(stbetabs)*0.6827)]### @ 68th percentile @
      res <- (b-a)/2 
    } else if(identical(tolower(str),"csbetaw")){ ###		@ CI-based sd(betaW)  @
      stbetaws <- stbetaWs <- eiread(dbuf,"stbetaws",formula,calculate=FALSE,...)   
      a <- stbetaWs[,floor(cols(stbetaws)*0.3413)] ### @ 34th percentile @
      b <- stbetaWs[,floor(cols(stbetaws)*0.6827)] ### @ 68th percentile @
      res <- (b-a)/2 
    } else if(identical(tolower(str),"gebw")){ ###                 @ betaB~betaB for sims betab>=betaw @
      betaBs <- eiread(dbuf,"betabs",formula,calculate=FALSE,...)  
      betaWs <- eiread(dbuf,"betaws",formula,calculate=FALSE,...)  
      a <- betaBs < betaWs
   
      betaBs[as.logical(a)] <- NA
    ###  betaBs <- mkmissm(betaBs,a)
      betaWs[as.logical(a)] <- NA
   ###   betaWs <- mkmissm(betaWs,a);
      res <- cbind(meanwc(t(as.matrix(betaBs)),1),meanwc(t(as.matrix(betaWs)),1), colSums(as.matrix(1-t(a))))
    }  else if(identical(tolower(str),"gebwa")){ ###                @ B^b ~ B^w for sims betaB >= betaW @
      a <- eiread(dbuf,"gebw",formula,calculate=FALSE,...)  
      res <- cbind(meanwc(a[,1],a[,3]),meanwc(a[,2],a[,3]))
    } else if(identical(tolower(str),"gewb")){ ###  @ betaB~betaW for sims betaW>=betaB @
      betaBs <- eiread(dbuf,"betabs",formula,calculate=FALSE,...)  
      betaWs <- eiread(dbuf,"betaws",formula,calculate=FALSE,...)  
      a <- betaBs> betaWs
      betaBs[as.logical(a)] <- NA
      betaWs[as.logical(a)] <- NA
      res <- cbind(meanwc(t(as.matrix(betaBs)),1), meanwc(t(as.matrix(betaWs)),1), colSums(as.matrix(1-t(a))))
    } else if(identical(tolower(str),"gewba")){###                 @ B^b ~ B^w for sims betaW >= betaB @
      a <- eiread(dbuf,"gewb",formula,calculate=FALSE,...)  
      res <- cbind(meanwc(a[,1],a[,3]), meanwc(a[,2],a[,3]))
    } else if (identical(tolower(str),"bounds")){ ###		@ compute precinct bounds  @
      t <- vread(dbuf,"t")
      x <- vread(dbuf,"x")
      n <- eiread(dbuf,"n",formula,calculate=FALSE,...)  
      if(exists("evbase")) EnumTol <- get("EnumTol", env=evbase)
      else  EnumTol <- 1.e-4
        
      lst <- bounds1(t,x,n,EnumTol)
      res <- lst$bs
      a <- lst$aggs
      
    } else if(identical(tolower(str),"bounds2")){	###	@ compute precinct bounds for ei2 @
      res <- NA
      if(vin(dbuf,"undert")){
         v <- vread(dbuf,"t")
         x <- vread(dbuf,"underx")
         n <- eiread(dbuf,"undern",formula,calculate=FALSE,...)  
         t <- eiread(dbuf,"undert",formula,calculate=FALSE,...)  
         if(exists("evbase")) EnumTol <- get("EnumTol", env=evbase)
         else  EnumTol <- 1.e-4
         lst <- bounds2(v,t,x,n,EnumTol)
         res <- lst[[1]]
         a <- lst$aggs
       }
    } else if(identical(tolower(str),"abounds")) { ###		@ compute aggregate bounds  @
      t <- vread(dbuf,"t")
      x <- vread(dbuf,"x")
      n <- eiread(dbuf,"n",formula,calculate=FALSE,...)  
      if(exists("evbase")) EnumTol <- get("EnumTol", env=evbase)
      else  EnumTol <- 1.e-4
      lst <- bounds1(t,x,n,EnumTol)
      print(names(lst))
      a <- lst[[1]]
      res <- lst[[2]]
   
   if(Eprt>0){
  ##    vrs <- c("       ", betaBb, etaW, "Aggregate bounds")
  ##    print(vrs)
      message("Lower: ", res[,1])
      message("Upper: ", res[,2])
    }
    res <- t(res)
      
    } else if (identical(tolower(str),"abounds2")){###		@ compute aggregate bounds for ei2  @
      res <- NA
    if(vin(dbuf,"undert")){
      v <- vread(dbuf,"t")
      x <- vread(dbuf,"underx")
      n <- eiread(dbuf,"undern",formula,calculate=FALSE,...)  
      t <- eiread(dbuf,"undert",formula,calculate=FALSE,...)  
      if(exists("evbase")) EnumTol <- get("EnumTol", env=evbase)
      else  EnumTol <- 1.e-4
      lst <- bounds2(v,t,x,n,EnumTol)
      a <- lst[[1]]
      res <- lst[[2]]
      
    }
     
    if(Eprt>0){
   ##   vrs <- c("       ", lambdaB, lambdaW, "Aggregate bounds")
   ##   print(t(vrs))
      message("Lower: ", res[,1])
      message(" Upper: ", res[,2])
    }
    res <- t(res);
    } else if(identical(tolower(str), "pphi")){ ###			@ prints phi and se's  @
      b <- eiread(dbuf,"phi",formula,calculate=FALSE,...)  
      res <- NA
     if(!scalmiss(b)){
       e <- eiread(dbuf,"etas",formula,calculate=FALSE,...)  
       ei.vc <- eiread(dbuf,"ei.vc",formula,calculate=FALSE,...)  
       ind <- as.vector(eiread(dbuf, "ghactual",formula,calculate=FALSE,...)  )
       if(ind <= 0) stop("eiread: ghactual not found")
      if (ei.vc[ind,1]!=-1)
        a <- rbind(sqrt(extract.diag(vread(dbuf,"vcphi"))),e)
      else
        a <- matrix(NA, nrow=rows(b),ncol=1)
     
      res <- rbind(t(as.matrix(b)), t(a));
      if(Eprt>0){
        message("Maximum likelihood results in scale of estimation (and se's)")
         cat("? ")
        if( vin(dbuf,"parnames")){
          a <-t( vread(dbuf,"parnames"))
          print(a)
          
        }
	 
	 print(res)
      }
     }
    }else if(identical(tolower(str),"psiu")){ ###			@ untruncated psi  @
      if(!vin(dbuf,"phi"))
        return(NA)
      b <- eiread(dbuf,"phi",formula,calculate=FALSE,...)  
      if(all(is.na(b))) return(NA)
      res <- read.par(dbuf,b, par="psiu",mess="Untruncated psi's", evbase)
      
      
          
    } else if(identical(tolower(str), "mppsiu")){###		@ Mean Posterior untruncated psi  @
      b <- as.matrix(eiread(dbuf,"phisims",formula,calculate=FALSE,...)  )
      res <- NA
    if(!scalmiss(b)){
      if (cols(b)==2)##  @ i.e., if _EisChk @
        b <- b[,1]
      else
        b <- colMeans(b)
     res <- read.par(dbuf,b,"mppsiu", "Mean Posterior Untruncated psi's", evbase)
    } 
    }else if(identical(tolower(str), "psi")){###			@ ultimate truncated psi  @
      b <- as.matrix(eiread(dbuf,"phi",formula,calculate=FALSE,...)  )
      res <- NA
      if(!scalmiss(b)){
      ###     clearg _Erho;
      Erho <- eiread(dbuf,"Erho",formula,calculate=FALSE,...)  
      assign("Erho", Erho, env=evbase)
      if(Erho[1]==0){
        c <- eiread(dbuf,"parnames",formula,calculate=FALSE,...)  
        cl <- sapply(c,tolower)
        e <- unlist(sapply(cl, identical, "rho"))
        b <- subset(b, subset=!e)
      }else b <- as.vector(b)
      Zb <- as.vector(eiread(dbuf,"Zb",formula,calculate=FALSE,...)  )
      Zw <- as.vector(eiread(dbuf,"Zw",formula,calculate=FALSE,...)  )
      x <- vread(dbuf,"x")
      Ez <- eiread(dbuf,"Ez",formula,calculate=FALSE,...)  
      assign("Ez", Ez, env=evbase)
      res <- eirepart(b,Zb,Zw,as.vector(x),Ez)
    
      if(Eprt>0){
        vrs <- cbind('bb', 'bw', 'sb', 'sw', 'rho')
        message("Truncated psi's (ultimate scale)")
        print(vrs)
        print(as.vector(res))
      }
    }
    }else if(identical(tolower(str), "aggs")){###			@ sims x 2 of wtd mean of (betaB&W)  @
      if(vin(dbuf,"x2")){
        nb <- eiread(dbuf,"nb2",formula,calculate=FALSE,...)  
        nw <- eiread(dbuf,"nw2",formula,calculate=FALSE,...)  
      }else{
        nb <- eiread(dbuf,"nb",formula,calculate=FALSE,...)  
        nw <- eiread(dbuf,"nw",formula,calculate=FALSE,...)  
      }
    res <- cbind(meanWc(vread(dbuf,"betabs"),nb), 
                 meanWc(eiread(dbuf,"betaws",formula,calculate=FALSE,...),nw))
    }else if(identical(tolower(str),"maggs")){###		        @ 2 x 1: meanc(aggBs~aggWs) @
      a <- eiread(dbuf,"aggs",formula,calculate=FALSE,...)  
      res <- colMeans(a)
    } else if(identical(tolower(str),"vcaggs")){###		@ 2 x 2: vcx(aggBs~aggWs)  @
      a <- eiread(dbuf,"aggs",formula,calculate=FALSE,...)  
      res <- var(a)
    } else if(identical(tolower(str), "paggs")){##			@ 2 x 2: ests (se's)  @
      a <- eiread(dbuf,"Maggs",formula,calculate=FALSE,...)  
      b <- eiread(dbuf,"VCaggs",formula,calculate=FALSE,...)        
      b <- sqrt(extract.diag(b))
      res <- rbind(t(as.matrix(a)), t(as.matrix(b)))
      if(Eprt>0){
      message("Estimates of Aggregate Quantities of Interest")
      vrs <- cbind('betab', 'betaw')
      print(as.vector(vrs))
      print(res)
    }
    } else if(identical(tolower(str),"goodman")){###		@ 2 x 2 of Goodman's ests (se's)  @
  ###  clearg _Routput,_Rconst;
      assign("Routput", 0,env=evbase)
      assign("Rconst", 0,env=evbase)
      x <- vread(dbuf,"x")
      omx <- 1-x
      t <- vread(dbuf,"t")
   
      #lstreg <- reg(formula,cbind(x,omx),t,...)
      lstreg <- lm(t ~ x+omx-1)
      res <- lstreg$coefficients
      #a <- summary(lstreg)$sigma
      a <- summary(lstreg)$coefficients[,2]

      res <- rbind(as.vector(res), as.vector(a))
      if(Eprt>0){
        message("Goodman's Regression")
        vrs <- cbind('betab', 'betaw')
        print(as.vector(vrs))
        print(as.vector(res))
      }
    } else if(identical(tolower(str), "double")){###		@ 2 x 1 of Double regression ests @
### clearg _Routput,_Rconst;
      assign("Routput", 0,env=evbase)
      assign("Rconst", 0,env=evbase)
      res <- NA
    if(vin(dbuf,"undert")){
      t <- vread(dbuf,"undert")
      x <- vread(dbuf,"underx")
      v <- eiread(dbuf,"t",formula,calculate=FALSE,...)  
      omx <- cbind(x,(1-x))
      lstreg <- reg(formula,omx,t,...) 
      a <- lstreg$coefficients
      lstreg <- reg(formula,omx,v%dot*%t,...)
      b <- lstreg$coefficients
      res <- b%dot/%a
      if(Eprt>0){
        message("Double Regression")
        vrs <- rbind('lambdaB', 'lambdaW')
        print(as.vector(vrs))
        print(as.vector(res))
      }
    }
     
    } else if(identical(tolower(str), "neighbor")){ ###             @ 2 x 1 of Neighborhood estimates @
      t <- eiread(dbuf,"t",formula,calculate=FALSE,...)  
      b <- eiread(dbuf,"Nb",formula,calculate=FALSE,...)  
      a <- eiread(dbuf,"Nw",formula,calculate=FALSE,...)  
      res <- rbind(meanwc(t,b),meanwc(t,a))
      if(Eprt>0){
        message("Freedman et al.'s Neighborhood Model Estimates");
        vrs <- rbind('betab', 'betaw')
        print(as.vector(vrs))
        print(as.vector(res))
    }
    } else if(identical(tolower(str),"thomsen")){###		@ 2 x 1 of Thomsen's estimates @
      x <- vread(dbuf,"x")
      t <- vread(dbuf,"t")
      a <- (x==0)| (x==1) | (t ==0) | (t==1)
      invX <- qnorm(subset(x, subset=x!=a))
      invT <- qnorm(subset(t, subset=t!=a))
      meanX <- -colMeans(as.matrix(invX))
      meanT <- -colMeans(as.matrix(invT))
      rho <- cor(cbind(invX,invT))
      rho <- rho[1,2]
      p00 <- cdfbvn(meanX,meanT,rho) ##bivariate normal pmvnorm in probs.R 
      p10 <- pnorm(meanT)-p00
      p01 <- pnorm(meanX)-p00
      p11 <- 1-p00-p10-p01
      bb <- p11/(p11+p10)
      bw <- p01/(p01+p00)
      res <- rbind(bb,bw)
      if(Eprt>0){
        message("Thomsen's Ecological Logit Approach Estimates");
        vrs <- cbind('betab','betaw');
        print(as.vector(vrs))
        print(as.vector(res))
      }
    } else if(identical(tolower(str),"palmquist")){###		@ Palmquist's inflation factor  @
      x <- vread(dbuf,"x")
      n <- eiread(dbuf,"n",formula,calculate=FALSE,...)  
      b <- meanwc(x,n)
      a <- (meanwc(x^2,n)-b^2)/(b*(1-b))
      res <- (1/a)-1
      if (Eprt>0)
        message("Palmquist's Inflation Factor: ", res)
    }else if(identical(tolower(str),"truthb")){###		@ true betab  @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first")
      return(NA)
    }
    truth <- eiread(dbuf,"truth",formula,calculate=FALSE,...)  
    res <- truth[,1]
  }  else if(identical(tolower(str),"truthw")){###		@ true betaw @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first");
    
      return(NA)
    }
    truth <- eiread(dbuf,"truth",formula,calculate=FALSE,...)  
    res <- truth[,2]
  }  else if(identical(tolower(str),"nbv") || identical(tolower(str),"nbt")){###	@ number of blacks who Turnout @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first");
      return(NA);
    }
    b <- eiread(dbuf,"truthb",formula,calculate=FALSE,...)  
    nb <- eiread(dbuf,"nb",formula,calculate=FALSE,...)  
    res <- nb%dot*%b
  } else if(identical(tolower(str),"nbn")){###			@ number of blacks who don't vote @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first");
      return(NA);
    }
    b <- eiread(dbuf,"truthb",formula,calculate=FALSE,...)  
    nb <- eiread(dbuf,"nb",formula,calculate=FALSE,...)  
    res <- nb%dot*%(1-b)
  } else if(identical(tolower(str),"nwv") || identical(tolower(str),"nwt")){###	@ number of whites who Turnout @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first");
      return(NA);
    }
    b <- eiread(dbuf,"truthw",formula,calculate=FALSE,...)  
    nw <- eiread(dbuf,"nw",formula,calculate=FALSE,...)  
    res <- nw%dot*%b
  } else if(identical(tolower(str),"nwn")){###			@ number of blacks who don't vote @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first");
      return(NA)
    }
     
    b <- eiread(dbuf,"truthw",formula,calculate=FALSE,...)  
    nw <- eiread(dbuf,"nw",formula,calculate=FALSE,...)  
    res <- nw%dot*%(1-b)
  } else if(identical(tolower(str),"aggtruth")){###		@ aggregate truths  @
    if(!("truth" %inG% cv)){
      message("eiread: truth needs to be stored first")
      return(NA)
    }
   
    n <- eiread(dbuf,"n",formula,calculate=FALSE,...)  
    x <- vread(dbuf,"x")
    nb <- x%dot*%n
    nw <- (1-x)%dot*%n
    res <- vread(dbuf,"truth")
    res <- rbind(meanwc(res[,1],nb), meanwc(res[,2],nw))
    if(Eprt>0){
      message("Aggregate Truth")
      vrs <- cbind('betab', 'betaw')
      print(as.vector(vrs))
      print(as.vector(res))
    }
  } else if(identical(tolower(str),"psitruth")){###		@ true psi's, truncated scale  @
    a <- as.matrix(eiread(dbuf,"truth",formula,calculate=FALSE,...)  )
    res <- NA
    if(!scalmiss(a)){
      b <- cor(na.omit(a));
      res <- rbind(colMeans(na.omit(a)), sd(as.data.frame(na.omit(a))),b[2,1])
      if (Eprt>0){
        message("TRUE truncated psi's (ultimate scale)")
        vrs <- cbind('bb', 'bw', 'sb', 'sw', 'rho')
       
        print(as.vector(vrs))
        print(as.vector(res))
      }
    }
    }else if(identical(tolower(str),"tsims")){###			@ sims from p(T|X=seqas(0,1,100))  @
      c <- eiread(dbuf,"Eeta",formula,calculate=FALSE,...)  
    if(scalzero(c) && 
      (!(scalone(eiread(dbuf,"Zb",formula,calculate=FALSE,...)))
       || !(scalone(eiread(dbuf,"Zw",formula,calculate=FALSE,...) )))){
      message("eiread: tsims only works without covariates.")
      return(NA);
    }
      a <- 100
      x <- seqase(0,1,a)
      b <- eiread(dbuf,"phi",formula,calculate=FALSE,...) 
      res <- NA
      if(!scalmiss(b)){
    
      if (c==1 || c==4)
        b[rows(b)-1] <- b[2]
      else if( c==2 || c==5)
        b[rows(b)] <- b[3]
      else if( c==3 && rows(c)==1){
        b[rows(b)-1] <- b[2]
        b[rows(b)] <- b[4]
      }
      Ez <- eiread(dbuf,"Ez",formula,calculate=FALSE,...) 
      assign("Ez", Ez, env=evbase)
      zb <- eiread(dbuf,"Zb",formula,calculate=FALSE,...) 
      zw <- eiread(dbuf,"Zw",formula,calculate=FALSE,...) 
      lst <- eirepar(b,0,0,0,Ez)
      bb0 <- lst$Bb
      bw0 <- lst$Bw
      sb <- lst$sb
      sw <- lst$sw
      rho <- lst$rho
      res <- matrix(0,nrow=a,ncol=Esims)
      bnds <- matrix(c(0,1,0,1),nrow=2, byrow=T)
      e <- colMeans(as.matrix(vread(dbuf,"x")))
      for(i in 1:a){
        bb <- bb0+b[rows(b)-1]*(x[i+0]-e);
        bw <- bw0+b[rows(b)]*(x[i+0]-e);
        bs <- rndbtn(as.vector(bb),as.vector(bw),as.vector(sb),as.vector(sw),as.vector(rho),bnds,Esims)
        t <- bs[,1]%dot*%x[i+0]+ bs[,2]%dot*%(1-x[i+0])
        res[i+0,] <- t(as.matrix(sortc(t,1)))
      }
      res <- cbind(x, res)
    }

    } else if( identical(tolower(str), "tsims0")){ ###          @ sims from p(T|X_obs,Z_obs) @
      b <- eiread(dbuf,"phi",formula,calculate=FALSE,...) 
      if(scalmiss(b)) return(res <- NA)
      
      Ez <- eiread(dbuf,"Ez",formula,calculate=FALSE,...) 
      t <- eiread(dbuf,"t",formula,calculate=FALSE,...) 
      x <- eiread(dbuf,"x",formula,calculate=FALSE,...) 
      zb <- eiread(dbuf,"Zb",formula,calculate=FALSE,...) 
      zw <- eiread(dbuf,"Zw",formula,calculate=FALSE,...) 
      a <- ifelse(is.matrix(x), rows(x), length(x))
      
      lst <- eirepar(b,zb,zw,as.vector(x),Ez)
      bb <- lst$Bb
      bw <- lst$Bw
      sb <- lst$sb
      sw <- lst$sw
      rho <- lst$rho
      c <- eiread(dbuf,"Esims",formula,calculate=FALSE,...) 
      res <- matrix(0, nrow=a,ncol=c);
      bnds <- matrix(c(0, 1, 0, 1),nrow=2, ncol=2, byrow=T)
      for (i in 1:a){
        bs <- rndbtn(bb[i],bw[i],sb,sw,rho,bnds,c)
        p <- bs[,1]%dot*%x[i+0]+bs[,2]%dot*%(1-x[i+0])
        res[i+0,] <- t(as.matrix(sortc(p,1)))
      }
      res <- cbind(t,res)
      
    }else if(identical(tolower(str),"expvarci")){ ###		@ x~20%CI~mean~80%ci @
      if( rows(eiread(dbuf,"Eeta",formula,calculate=FALSE,...) )==4 &&
      (!scalone(eiread(dbuf,"Zb",formula,calculate=FALSE,...) )
       || !scalone(eiread(dbuf,"Zw",formula,calculate=FALSE,...) ))){
      message("eiread: expvarci only works without covariates.");
      return(NA)
    }
      b <- eiread(dbuf,"tsims",formula,calculate=FALSE,...) 
      res <- NA
    if(!scalmiss(b)){
      x <- as.matrix(b[,1])
      b <- as.matrix(b[,2:(Esims+1)])
      res <- cbind(x,b[,floor(0.2*Esims)], colMeans(t(b)), b[,floor(0.8*Esims)])
    }
    } else if(identical(tolower(str), "expvarci0")){###		@ t~20%CI~mean~80%ci @
      b <- eiread(dbuf,"tsims0",formula,calculate=FALSE,...) 
      c <- eiread(dbuf,"Esims",formula,calculate=FALSE,...) 
      res <- NA
      if (!scalmiss(b)){
        t <- as.matrix(b[,1])
        b <- as.matrix(b[,2:(c+1)])
        res <- cbind(t, b[,floor(0.2*c)], colMeans(t(b)), b[,floor(0.8*c)])
      }
    
  }else if(identical(tolower(str), "expvarcis")){###		@ x~20%CI~mean~80%ci LOESS smoothed @
    b <- eiread(dbuf,"expvarci",formula,calculate=FALSE,...) 
    res <- NA
    if(!scalmiss(b)){
      e <- output
      output <- 0
      loess.WgtType <- rep(2, length(b[,1]))
      loess.span <- 0.6667
      loess.degree <- 1
      colnames(b) <- paste("b", 1:ncol(b),sep="")
 
      y.loess <- loess(b2~b1, as.data.frame(b), weights= loess.WgtType ,span=loess.span,degree=loess.degree)
### predicted values for b2:Yhat
      yhat <- y.loess$fitted
### independent variable coordinates(vs b1) 
      c <- y.loess$x
    
      res <- cbind(c,yhat)
      y.loess <- loess(formula=b3 ~ b1, as.data.frame(b),weights= loess.WgtType,span= loess.span,degree=loess.degree)
      yhat <- y.loess$fitted
      c <- y.loess$x
  
      res <- cbind(res,yhat);
      y.loess <- loess(b4 ~ b1, as.data.frame(b), weights= loess.WgtType,span= loess.span,degree=loess.degree)
      yhat <- y.loess$fitted
      c <- y.loess$x
  
      res <- cbind(res,yhat)
      
      output <- e
    }

  } else if(identical(tolower(str),"sum")){ ###			@ prints all printable items @
    if(Eprt<1)
      Eprt <- 1
    if("titl    " %inG% cv){
      titl <- vread(dbuf,"titl")
      print(paste("** ",titl," **", sep=""))
      if(identical(titl,"*DB* Data Buffer from eimodels_avg() *DB*")){
        postp <- vread(dbuf,"postprob")
        priorp <- vread(dbuf,"prprob")
        mrgllk <- vread(dbuf,"margllik")
       ## ?; 
        message("The number of model averaged:", rows(vread(dbuf,"postprob")))
       ## ?;
        message("_EI_bma_est: ", ei.bma.est)
       ## ?;
        message("Model  Posterior Prior  Marginal")
        message("Number   Prob     Prob   LogLik ", cbind(postp, priorp[,2], mrgllk[,2]))
       ## ?;
        abounds <- eiread(dbuf,"abounds",formula,calculate=FALSE,...)  
        assign("abounds",abounds,env=evbase)
        ##?;
        paggs <- eiread(dbuf,"paggs",formula,calculate=FALSE,...)  
         assign("paggs",aggs,env=evbase)
        res <- ""
        return(res)
      }
    }
  
    if (scalone(eiread(dbuf,"Enonpar",formula,calculate=FALSE,...))){
      message("Nonparametric Estimation");
      message("EnonNumInt:   ", eiread(dbuf,"EnonNum",formula,calculate=FALSE,...))
      message("EnonEval:     ", eiread(dbuf,"EnonEva",formula,calculate=FALSE,...))
      message("N:             ", rows(vread(dbuf,"x")))
      message("Esims:        ", vread(dbuf,"Esims"))
     ## ?;
    }else{
    
      message("CML return: ", vread(dbuf,"retcode"),"     ")
      message("N:          ", rows(vread(dbuf,"x")), "      ")
      message("Esims:     ", vread(dbuf,"Esims"))
      if ("ebeta  " %inG% cv)
        message("Ebeta      ", vread(dbuf,"Ebeta"),"     ")
     
      message("Esigma:    ", vread(dbuf,"Esigma"),"       ")
      message("Erho:      ", as.vector((vread(dbuf,"Erho"))))
      message("Eisn:      ", vread(dbuf,"Eisn"),"     ")
      message("resamp:     ",eiread(dbuf,"resamp",formula,calculate=FALSE,...))
      message("GhActual:  ",vread(dbuf,"ghactual"))
      message("Estval:    ",vread(dbuf,"Estval"))
      if("eeta   " %inG% cv)
        message("Eeta:      ", as.vector(eiread(dbuf,"Eeta",formula,calculate=FALSE,...) ))
      
      message("log-likelihood:         ", vread(dbuf,"loglik"))
      message("ln(mean(Imptce Ratio)): ", eiread(dbuf,"meanIR",formula,calculate=FALSE,...) )
     ## ?;
      pphi <- eiread(dbuf,"pphi",formula,calculate=FALSE,...)  
      assign("pphi", pphi, env=evbase)
     ## ?;
      psiu <- eiread(dbuf,"psiu",formula,calculate=FALSE,...)  
      assign("pphi", psiu, env=evbase)
      ##?;
      psi <- eiread(dbuf,"psi",formula,calculate=FALSE,...)  
      assign("pphi", psi, env=evbase)
      ##?;
    
    }
      if ("truth   "%inG% cv)
        {
          psitruth <- eiread(dbuf,"psitruth",formula,calculate=FALSE,...)  
          assign("psitruth", psitruth, env=evbase)
        ##?;
          aggbias <- eiread(dbuf,"aggbias",formula,calculate=FALSE,...)  
          assign("aggbias", aggbias, env=evbase)
          ##?;
          eiread(dbuf,"eaggbias",formula,calculate=FALSE,...)  
           assign("eaggbias", eagbias, env=evbase)
          ##?;
          coverage <- eiread(dbuf,"coverage",formula,calculate=FALSE,...)  
           assign("coverage", coverage, env=evbase)
###?; 
          eiread(dbuf,"aggtruth",formula,calculate=FALSE,...)  
           assign("aggtruth", aggtruth, env=evbase)
###?;
        }
    abounds <- eiread(dbuf,"abounds",formula,calculate=FALSE,...)  
     assign("abounds", abounds, env=evbase)
    ###?;
    paggs <- eiread(dbuf,"paggs",formula,calculate=FALSE,...)  
     assign("paggs", paggs, env=evbase)
    res <- ""
    
  }else{
    if (Eprt>0)
      message("eiread: no such name, ",str)
    
    res <- NA
  
  }
  
  
  return(res)
}
  ###helper to eiread to avoid repetition
  read.par <- function(dbuf,b, par, mess, evbase=NULL,...){ 
  
    Erho <- eiread(dbuf,"Erho",formula=NA,calculate=FALSE,...)  
    if(length(evbase))
      assign("Erho", Erho, env=evbase)  
    if (Erho[1]==0){
      c <- eiread(dbuf,"parnames",formula=NA,calculate=FALSE,...) 
      cl <- sapply(c,tolower)
      ss <- unlist(sapply(cl,identical,"rho"))
        
      if(length(ss))
        b <- as.vector(subset(b, subset=!ss))
    }else b <- as.vector(b)
   
    Zb <- as.vector(vread(dbuf,"Zb"))
    Zw <- as.vector(vread(dbuf,"Zw"))
    x <- as.vector(vread(dbuf,"x"))
   
    Ez <- eiread(dbuf,"Ez",formula=NA,calculate=FALSE,...) 
    nm <- names(dbuf)
    dbuf <- c(dbuf,list(Ez))
    names(dbuf) <- c(nm, "ez")
    if(length(evbase))
      assign("Ez", Ez, env=evbase)
         
    lst <- eirepar(b,Zb,Zw,x,Ez)
    
    bb <- as.matrix(lst$Bb)
    bw <- as.matrix(lst$Bw)
    sb <- as.matrix(lst$sb)
    sw <- as.matrix(lst$sw)
    rho <- as.matrix(lst$rho)
    res <- rbind(bb,bw,sb,sw,rho)
    Eprt <- get("Eprt", env=evbase)
    if(Eprt>0){
      message(mess)
      a <- rbind(colMeans(bb),colMeans(bw),sb,sw,rho)
      vrs <- rbind(bb, bw, sb, sw, rho)
   ###   print(vrs)
    print(as.vector(a))
    }
  }

##
## checkr(dbuf,eps)
##
## Procedure to calculate value of R function  +/- eps
## for each parameter, holding other parameter values at their MLEs and
## report whether cdfbvn is sufficiently precise for each parameter.
## 
## Inputs: dbuf = EI data buffer
##         eps  = tolerance check (probably use _EdirTol)
## 
## Output: rchk, rows(phi)x2 matrix with rows correspnding to phi, cols
## corresponding to slightly less~more (by eps) than the MLEs, and
## each element indicating that the CDFBVN function is sufficiently
## precise (when 1) and insufficiently precise (when 0)
##
checkr <- function(dbuf,eps,Ez,...){
###  local phi,R,rr,zb,zw,x,y,k,kk,loparms,hiparms,loR,hiR,rchk;

  phi <- eiread(dbuf,"phi",formula=NA,calculate=FALSE,...) 
  rr <- rows(phi)
  lst = pluckdta(eiread(dbuf,"dataset",formula=NA,calculate=FALSE,...) ,Ez=Ez)
  Zb <- zb <- lst$Zb
  Zw <- zw <- lst$Zw
  x <- lst$x
  y <- lst$y
  loR <- matrix(1,nrow=rr-2,ncol=1)
  hiR <- matrix(1, nrow=rr-2,ncol=1)
  loparms <- phi-eps
  hiparms <- phi+eps
  lst <- eirepar(phi,zb,zw,x,Ez)
  bb <- lst[[1]]
  bw <- lst[[2]]
  sb <- lst[[3]]
  sw <- lst[[4]]
  rho <- lst[[5]]
  
  R <- colSums(as.matrix(na.omit(lncdfbvnu(bb,bw,sb,sw,rho))))

  for(kk in 1:(rr-2)){
    k <- kk;
    lst <- eirepar(loparms,zb,zw,x,Ez)
    bb <- lst[[1]]
    bw <- lst[[2]]
    sb <- lst[[3]]
    sw <- lst[[4]]
    rho <- lst[[5]]
    loR[k] <- colSums(as.matrix(na.omit(lncdfbvnu(bb,bw,sb,sw,rho))))
    lst <- eirepar(hiparms,zb,zw,x,Ez)
    bb <- lst[[1]]
    bw <- lst[[2]]
    sb <- lst[[3]]
    sw <- lst[[4]]
    rho <- lst[[5]]
    hiR[k] <- colSums(as.matrix(na.omit(lncdfbvnu(bb,bw,sb,sw,rho))))
  }

  rchk <- (cbind(loR,hiR) !=R)

return(rchk)
}

##/*
##** {b,se,vc,yhat,e,sig} = reg(x,y);
##**
##**  input:  x = NxK explanatory variable matrix (no constant term by default)
##**          y = Nx1 dependent variable
##**
##**  output: b = regression coefficients
##**          se = sqrt(diag(vc)), standard errors
##**          vc = variance-covariance matrix (see _Rrobust, below)
##**          yhat = X'b, predicted values
##**          e = y-yhat = residuals
##**          sig = sqrt(e'e/(n-k))
##**
##** globals:
##**          _Rweight = sd(y), the weight vector for Weighted Least Squares
##**          _RXnames = vector of x variable names
##**          _RYname  = vector of y variable names
##**          _Rfmt    = digits to the right of decimal point for output
##**          _Rselect = vector of 1's to select and 0's to delete
##**          _Rrobust = defines VC.  if -1 (default) use usual method;
##**                     if 0 use White's heteroskedasticity-consistent VC matrx;
##**                     if an integer >0, use MA(N)-time-series and
##**                        heteroskedasticity-consistent VC matrix
##**          _Rtheta = if _Rrobust>0 doesn't produce a non-neg.definate VC matrix
##**                    then set this to some value >0 and usually <1.  Possibly
##**                    also increase _Rrobust.  (Default _Rtheta=0)
##**          _Routput = 0 if none; 1 (default) if print output
##**          _Rconst =  1 if include constant term (default), 0 don't include
##**
##** example:     let _RYname = income;
##**              let _RXnames = educat race;
##**              call reg(educat~race,income);
##**
##*/
reg <- function(formula=NA,x,y=NULL,...)
{
  y <- if(length(y)) as.matrix(y)
  x <- as.matrix(x)
  if(length(y) && rows(x) != rows(y))
    stop("reg: rows in depvar and covariates do not agree")
  drvdot <- match.call(expand.dots=TRUE)
  drv    <-  match.call(expand.dots=FALSE)
  ev <- environment()
  entries <- expanddots(drvdot,drv,ev)
  
  vec <- c("subset","weights","na.action", "offset", "contrasts")
  if(length(entries) && !all(entries %in% vec))
    stop("reg:Invalid argument to lm")
 
  wght <- if("weights" %in% entries) get("weights",env=ev)
 
  subr <- if("subset" %in% entries) get("subset",env=ev)
   
  if("na.action" %in% entries)
    namiss <- get("na.action",env=ev)
  else
    namiss <- options()[["na.action"]]
    
  
  dat <- cbind(y,x)
  nc  <- cols(x)  
  nr  <- rows(x)
  nmx <- colnames(x)
  nmy <- if(length(y)) colnames(y)
  if(dim(y)[2] > 1) stop("reg: regression depvar only one column")
  if(!any(is.na(as.character(formula))) && length(c(nmx,nmy)) <=0)
    stop("reg:Formula requires named data for covariates x and depvar y")
  if(!any(is.na(as.character(formula))) && length(c(nmx,nmy)) >0 )
    {
      dat <- as.data.frame(dat)
      return(lm.reg <- lm(formula, data=dat,subset=subr,weights=wght,na.action=namiss))
      
    }
  if (length(y)){
      
      nms <- c("y", paste("x",1:nc,sep=""))
      xs  <- paste("x", 1:nc,sep="",collapse=" + ")
      ff <- paste("y ~",xs)
    
      colnames(dat) <- c(nmy,nms)
      dat <- as.data.frame(dat)
      lm.reg <- lm(formula=as.formula(ff), data=dat,subset=subr, weights=wght,na.action=namiss)
   ###   lm0 <<- lm.reg 
      return(lm.reg)
    }
    stop("reg:Provide formula and names data")
}
###DESCRIPTION S3 method that calls eiread
####           if only one argument, then calls "paggs"
###            if two arguments then calls eiread(object,str)
summary.ei <- function(object,...){
  
  drvdot <- match.call(expand.dots=TRUE)
  if(length(drvdot)>2) str <- drvdot[[3]]
  else str <- "paggs"
  return(eiread(object, str))
}
