library(mvtnorm)


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
eiread <- function(dbuf, str, compute =FALSE){

  if(!compute)
   return(vread(dbuf,str))
  evbase  <- get("evbase", env=parent.frame())
  getEnvVar(evbase, environment())###, vecvar=c("eimetar"))
  eimetar <- get("eimetar", env=evbase)

  if(vin(dbuf, "titl")){
    titl <- vread(dbuf, "titl")
    nc  <- nchar(as.character(floor(eimetar))) 
    fmt <- formatC(eimetar,width=nc,digits=0, format="f")
    if (identical(titl, "*MDB* Meta-Data Buffer from 2nd Stage *MDB*")){ 
      dbuf <- vread(dbuf,paste("dbuf",fmt,sep=""))
    }else if(identical(titl, "*MDB* Meta-Data Buffer from eimodels_def() *MDB*")){
      if(!vin(dbuf,paste("mod.d",fmt, sep=""))){
        message("eiread: Model " + eimetar + " is not stored in this data buffer.");
        return(NA)
      }
      dbuf <- vread(dbuf,paste("mod.d",fmt,sep=""))
      if(Eptr==3)
        message("Reading ", str, " from Model ", eimetar, "...");
       
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
###  cv <- vnamecv(dbuf)
  cv <- names(dbuf)
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
               ,"psi", "psitruth", "psiu", "r",  "ri", "retcode", "tsims","vcphi");
      vrs <- as.matrix(tolower(vrs))
     if(identical(str, "esims")){
       res <- ncol(as.matrix(eiread(dbuf,"betabs")))
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
               "thomsen", "tsims", "vcaggs", "vcphi");
      
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
    res <- vread(dbuf, str) ##returns lement or NA
    res <- as.matrix(res)
  }
  if(identical(str, "eeta")){
   
    if(scalmiss(res))
      res <- matrix(0, nrow=4, ncol=1)
    else if(ncol(as.matrix(res)) == 2)
      res <- matrix(as.vector(res), ncol=1)
    else if(rows(res)==2)
      res <- rbind(res, matrix(0, nrow=2, ncol=1))
  }else if(identical(str, "zb")){
    e <- eiread(dbuf, "eeta")
    if(e %in% c(1,3,4))
      res <- vread(dbuf, "x")
         
  }else if(identical(str, "zw")){
    e <- eiread(dbuf, "Eeta")
 
    if(e %in% c(5,2,3))
      res <- vread(dbuf,"x")
         
  }else if(identical(str, "titl") || identical(str, "undertitle") && is.na(res)){
    res <- ""
  
     ###      /***** stored globals *****/
  }else if(strc %inG% cv){
    res <- vread(dbuf,str)
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
      for(n in 1:a)
        res[n, ] <- res[n, order(runif(c))]
    }else
      message("eiread: 'x2' option is only available in data buffers created by ei2")
      
  ###else if identical(str, "emaxiter"): already taken care of  
    
  }else if(identical(str,"eigraph.bvsmth")){
    
    if(vin(dbuf, "bvsmth"))
      res <- vread(dbuf, "bvsmth") 
  }else if(identical(str,"eimodels.save")){
    
    if(vin(dbuf, "eimsave"))
      res <- vread(dbuf, "eimsave")
  }else if(identical(str,"ei.bma.est")){
    
    if(vin(dbuf, "eibmaest"))
      res <- vread(dbuf, "eibmaest")
  }else if(identical(str,"ei.bma.prior")){
    
    if(vin(dbuf, "prprob")){
      priorp <- vread(dbuf,"prprob")
      res <- priorp[, 2]
    }
  
  }else if(identical(str,"edoml.phi")){
    if(vin(dbuf, "doml.phi"))
      res <- vread(dbuf, "doml.phi")
    
  }else if(identical(str,"enoneval")){
       
    if(vin(dbuf, "enoneva"))
      res <- vread(dbuf, "enoneva")
    
  }else if( identical(str,"enonnumint")){
    if (vin(dbuf,"enonnum"))
      res <- vread(dbuf,"enonnum")
 
#### left checking here#######################################    
  }else if (identical(str, "meanir") && is.na(res)){
   
    if(vin(dbuf, "EmeanIR"))
      res <- vread(dbuf, "emeanir")
    else if(vin(dbuf, "lnir") && eiread(dbuf,"EisChk") == 0)
      res <- vread(dbuf, "lnir")
    else if(vin(dbuf, "lnir") && eiread(dbuf,"EisChk") == 1){
      a <- eiread(dbuf, "lnir")
      a <-  a[,1]
      max <- max(a)
      res <- max + log(meanc(exp(a-max)))
      
    }
         
  }else if (identical(str, "logliks")){

    a <- eiread(dbuf, "phi")
    b <- eiread(dbuf, "dataset")
    res <- eiloglik(a,b)
  }else if(identical(tolower(str), "resamp")){
         res <- NA
         if(vin(dbuf, "resamp") || vin(dbuf, "Eresamp"))
           res <- vread(dbuf, "resamp")
       
       }else if(identical(tolower(str, "enonpar"))){
         if(!vin(dbuf,"enonpar"))
           res <- 0
         else
           res <- vread(dbuf, "enonpar")
                
       }else if (identical(tolower(str),"etac")){
         res <- eiread(dbuf, "Eeta")
         if(rows(res) == 1)
           res <- matrix(0, nrow=2, ncol=1)
         else if(rows(res) == 3 && res[1]==4)
           res <- as.matric(c(0, res[2]))
         else if(rows(res) ==3 && res[1] == 5)
           res <- as.matrix(c(res[2], 0))
         else if(rows(res)==4)
           res <- as.matrix(res[1:2])
         
                       
       }else if(identical(tolower(str),"etas")){
         res <- eiread(dbuf, "Eeta")
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
         zb <- eiread(dbuf, "zb")
         zw <- eiread(dbuf,"zw")
         assign("Ez", 0, env=evbase)
         Ez <- as.matrix(c(cols(Zb) + 1 - (Zb == 1), cols(Zw) + 1 - (Zw == 1)))
         assign("Ez", Ez, env=evbase)
         res <- Ez

       }else if(identical(tolower(str),"ebounds")){
         dummy <- eiread(dbuf, "Ez")
         b <- eiread(dbuf,"Ebounds")
         if(scalzero(b))
           res <- matrix(c(-20, 20), nrow=1, ncol=2)
         else if(cols(b) == 2)
           res <- b
         else if (scalone(b)){
           e <- matrix(c(-10, 10), nrow=1, ncol=2)
           f <-  matrix(c(-20, 20), nrow=1, ncol=2) ## @{ -1e256 1e256 };@
           if(Ez[1] == 1)
             res <- e
           else
             res <-matrix(f, nrow=Ez[1], ncol=1)
           if(Ez[2]==1)
             res <- as.matrix(c(res, e))
           else
             res <- rbind(res,matrix(f, nrow=Ez[2], ncol=1))
           mat1 <-  matrix(c(-6, 3), nrow=1, ncol=2)
           mat2 <-  matrix(c(-2, 2), nrow=1, ncol=2)
           res <- rbind(res,mat1, mat1,mat2)
         }else
         message("eiread: problem with Ebounds")
       }else if(identical(tolower(str), "nobs")){
         res <- rows(vread(dbuf, "t"))
         if(Eprt > 0)
           message( "number of observations: ", res)
       }else if(identical(tolower(str), "tvap")){
         strc <- str
         if(strc %inG% cv)
           res <- vread(dbuf, "tvap")
         else
           res <- vread(dbuf, "n")
         
       }else if(identical(tolower(str), "n")){
         strc <- str
         if(strc %inG% cv)
           res <- vread(dbuf, "n")
         else
           res <- vread(dbuf, "tvap")
       }else if(identical(tolower(str), "bvap") || identical(tolower(str), "nb")){ ###	@ black vap  @{
         x <- vread(dbuf, "x")
         n <- eiread(dbuf, "n")
         res <- x * n
       }else if(identical(tolower(str), "nb2")){
###	@ number of blacks turning out to vote: sims from ei2 @
         x <- eiread(dbuf, "x2")
         n <- eiread(dbuf, "n")
         res <- x * n
       }else if( identical(tolower(str),"wvap") || identical(tolower(str),"nw")) {###@ white vap  @
         x <- vread(dbuf, "x")
         n <- eiread(dbuf, "n")
         res <- (1-x) * n
       } else if(identical(tolower(str),"nw2")){	###@ white vap  @
         x <- eiread(dbuf,"x2");
         n <- eiread(dbuf,"n");
         res <- (1-x)*n;

       } else if(identical(tolower(str),"nt")){		###	@ number of people who Turnout @
           t <- eiread(dbuf,"t");
           n <- eiread(dbuf,"n");
           res <- t*n;

       } else if(tolower(str,"dataset")){	###	@ dataset for input to eiloglik() @
         x <- vread(dbuf,"x");                  ### @ ignore _EselRnd @
         t <- vread(dbuf,"t");
         Zb <- eiread(dbuf,"Zb");
         Zw <- eiread(dbuf,"Zw");
         eiread(dbuf,"Ez");
         Eselect <- eiread(dbuf,"Eselect");
         assign("Eselect",Eselect, env=evbase) 
         res <- packdta(x,Zb,Zw,t);
        }else if (identical(tolower(str),"betaws")){	###	@ betaW simulations @
          t <- vread(dbuf,"t");
          if (vin(dbuf,"x2"))
            x <- vread(dbuf,"x2")
          else
            x <- vread(dbuf,"x")
    
          betaBs <- eiread(dbuf,"betabs")
          res <- NA
          if(!scalmiss(betabs))
            res <- betab2w(t,x,betaBs);
        } else if(identical(tolower(str),"beta")){ ###   @ E(betab)~E(betaw) by precinct @
          res <- cbind(eiread(dbuf,"betab"),eiread(dbuf,"betaw"))
        } else if(identical(tolower(str),"betab")){	###		@ E(betab) for each precinct @
          res <- colMeans(t(eiread(dbuf,"betabs")))
        } else if(identical(tolower(str),"betaw")){ ###			@ E(betaW) for each precinct @
          res <- colMeans(t(eiread(dbuf,"betaws")))
        }  else if (identical(tolower(str), "sbetab")){###		@ sd(betab) for each precinct @
          res <- sd(as.data.frame(t(eiread(dbuf,"betabs"))))
        } else if( identical(tolower(str), "sbetaw")){ ###		@ sd(betaW) for each precinct @
          res <- sd(as.data.frame(t(eiread(dbuf,"betaws"))))
        } else if(identical(tolower(str), "rnbetabs")){###		@ randomly permuted betabs sims @
          res <- eiread(dbuf,"betabs")
          a <- rows(res)
          c <- cols(res);
          for( i in 1:a)
            res[i,] <- res[i,order(runif(c))]   
   
        } else if(identical(tolower(str), "rnbetaws")){ ###		@ randomly permuted betabs sims @
          res <- eiread(dbuf,"betaws");
          a <- rows(res);
          c <- cols(res);
          for( i in 1:a)
            res[i,] <- res[i,order(runif(c))]
        } else if( identical(tolower(str), "stbetabs")){ ###		@ sorted betaB simulations  @
          betabs <- betaBs <- eiread(dbuf,"betaBs");
          res <- NA
          if(!scalmiss(betabs))
            res <- betabs <- betaBs <- sortbyRow(betabs)
          
        } else if(identical(tolower(str), "stbetaws")){ ###		@ sorted betaW simulations  @
          betaWs <- betaws <- eiread(dbuf,"betaWs");
          res <- NA
          if(!scalmiss(betaws))
             res <- betaws <- betaWs <- sortbyRow(betaws)
        
        }else if(identical(tolower(str,"truptile"))){###	@ percentile at which true value falls @
          if(!vin(dbuf,"truth")){
           message("eiread: truth needs to be stored first")
           return(NA)
         }
          
          stbetabs <- eiread(dbuf,"stbetabs");
          stbetaws <- eiread(dbuf,"stbetaws");
          truth <- eiread(dbuf,"truth");
          b <- cols(stbetabs);
          
          res <- minindc(abs(t(stbetabs-truth[,1])))/b;    
          res <- cbind(res, minindc(abs(t(stbetaws-truth[,2])))/b);
          res[,1] <- recode(res[,1],stdc(t(stbetabs)) <= EnumTol,0.5); ###@ homog prects @
          res[,2] <- recode(res[,2],stdc(t(stbetaws))<= EnumTol,0.5);
          
        }else if(identical(tolower(str),"ci50b")){### @ 50% confidence intervals for betab @
          stbetabs <- eiread(dbuf,"stbetabs");
          e <- cols(stbetabs);
          res <- st
        }else if(identical(tolower(str),"ci80b"))	{ ###@ 80% confidence intervals for betab @
          stbetabs <- eiread(dbuf,"stbetabs");
          e <- cols(stbetabs);
          res <- cbind(stbetabs[,floor(0.1*e)], stbetabs[,floor(0.9*e)]);

        } else if(identical(tolower(str), "ci95b")){ ###@ 95% confidence intervals for betab @
          stbetabs <- eiread(dbuf,"stbetabs");
          e <- cols(stbetabs);
          res <- cbind(stbetabs[,floor(0.05*e)],stbetabs[,floor(0.95*e)]);
        } else if(identical(tolower(str),"ci50w")){###	@ 50% confidence intervals for betaw @
          stbetaws <- eiread(dbuf,"stbetaws");
          e <- cols(stbetaws);
          res <- cbind(stbetaws[,floor(0.25*e)], stbetaws[,floor(0.75*e)]);
        } else if(identical(tolower(str),"ci80w")){###		@ 80% confidence intervals for betaw @
          stbetaws <- eiread(dbuf,"stbetaws");
          e <- cols(stbetaws);
          res <- cbind(stbetaws[,floor(0.1*e)],stbetaws[,floor(0.9*e)]);
        } else if(identical(tolower(str),"ci95w")){ ###	@ 95% confidence intervals for betaw @
          stbetaws <- eiread(dbuf,"stbetaws");
          e <- cols(stbetaws);
          res <- cbind(stbetaws[,floor(0.05*e)],stbetaws[,floor(0.95*e)]);
        }  else if(identical(tolower(str),"ci80bw")){ ###		@ 80% conf intervals for betab betaw  @
          res <- cbind(eiread(dbuf,"ci80b"), eiread(dbuf,"ci80w"));
        } else if(identical(tolower(str),"ci95bw")){###		@ 95% conf intervals for betab betaw  @
          res <- cbind(eiread(dbuf,"ci95b"),eiread(dbuf,"ci95w"));

        } else if(identical(tolower(str),"coverage")){ ###		@ CI coverage @
          truth <- eiread(dbuf,"truth");
          res <- NA
          if(!scalmiss(truth)){
            if(vin(dbuf,"x2")){
              nb <- eiread(dbuf,"nb2");
              nw <- eiread(dbuf,"nw2");
            }else{
              nb <- eiread(dbuf,"nb");
              nw <- eiread(dbuf,"nw");
            }
            a <- eiread(dbuf,"ci50b");
            bb <- is.na(colSums(t(cbind(a,truth))));
            a <- subset(a,subset=!bb);
            e <- subset(truth,subset=!bb);
            f <- (e[,1]>=a[,1]) & (e[,1] <=a[,2]);
            res <- rbind(colMeans(f),t(colMeans(meanwc(f,delif(nb,bb)))));
            a <- eiread(dbuf,"ci80b");
            bb <- is.na(colSums(t(cbind(a,truth))));
            a <- delif(a,bb);
            e <- delif(truth,bb);
            f <- (e[,1]>=a[,1]) & (e[,1]<=a[,2]);
            res <- cbind(res, rbind(colMeans(f), t(colMeans(meanwc(f,delif(nb,bb)))))); 
            a <- eiread(dbuf,"ci50w");
            bb <- is.na(colSums(t(cbind(a,truth))));
            a <- delif(a,bb);
            e <- delif(truth,bb);
            f <- (e[,2]>=a[,1]) & (e[,2] <=a[,2])
            res <- cbind(res, rbind(colMeans(f), t(colMeans(meanwc(f,delif(nw,bb))))));
            a <- eiread(dbuf,"ci80w");
            bb <- is.na(colSums(t(cbind(a,truth))));
            a <- delif(a,bb);
            e <- delif(truth,bb);
            f <- (e[,2]>=a[,1]) &(e[,2] <=a[,2]);
            res <- cbind(res, rbind(meanc(f), t(colMeans(meanwc(f,delif(nw,bb))))));
            if(Eprt>0){
              message("CI coverage; % true values within each confidence interval");
              vrs <- c("      ", "50%Black", "80%Black", "50%White", "80%White");
              vrs <- as.matrix(vrs)
              print(t(vrs))
              message("    %: ")
              print(res[1,]);
              message("Wtd %: ")
              print(res[2,]);
            }
          }
    }  else if(identical(tolower(str),"checkr")){###    @ check R function precision @
      Edirtol <- eiread(dbuf,"Edirtol");
      res <- checkr(dbuf,Edirtol);
    } else if(identical(tolower(str),"ri")) { ###                    @ ln(R) @    
      a <- eiread(dbuf,"phi");
      lst <- pluckdta(eiread(dbuf,"dataset"));
      Zb <- lst$Zb
      Zw <- lst$Zw
      x <- lst$x
      t <- lst$t 
      res <- na.omit(lncdfbvnu(eirepar(a,zb,zw,x)));
    } else if(identical(tolower(str),"r")){ ###                    @ sum(ln(R)) @    
      res <- colSums(eiread(dbuf,"r"));
    } else if(identical(tolower(str),"aggbias")){ ###		@ aggregation bias regressions @
      truth <- eiread(dbuf,"truth");
      res <- NA
    if(!scalmiss(truth)){
    
      x <- vread(dbuf,"x");
      assign("Routput",0, env=evbase)
      assign("Rconst",0, env=evbase)
      
      Rconst <- 1;
    ###  {b,bb,t,t,t,t}=reg(x,truth[.,1]);
      res <- cbind(b,bb)
    ###  {b,bb,t,t,t,t}=reg(x,truth[.,2]);
      res <- rbind(res,cbind(b,bb));
      if(Eprt>0){
        vrs <- as.matrix(c("TRUEDepV", "       ", "coeffs", "se's"));
        
        print(t(vrs));
        vrs <- as.matrix(c(const, x, const, x));
        a <- cbind(vrs,res);
        vrs <- as.matrix(c("betaB", "   ", "betaW", "   "));
        a <- cbind(vrs,a)
        b <- matrix(1, nrow=4,ncol=1);
        mask <- cbind(matrix(0,nrow=4,ncol=2), matrix(1,nrow=4,ncol=2));
	
      fmt <- matrix(c("-*.*s ", 8, 8,
                     "-*.*s ", 8, 8,
                     "*.*lf", 7, 4,
                     "*.*lf", 7, 4), nrow=4, ncol=3)
 ###       call printfm(a,mask,fmt);		  
      }
    }
    }else if(identical(tolower(str),"eaggbias")){###	@ estimated aggregation bias regressions @
      beta <- eiread(dbuf,"betab");
      betaw <- eiread(dbuf,"betaw");
      res <- NA
    if(!scalmiss(betaw)){
      x <- eiread(dbuf,"x");
      assign("Routput",0, env=evbase)
      assign("Rconst",0, env=evbase)
      Rconst <- 1;
      
###      {b,bb,t,t,t,t}=reg(x,betab);
      res <- cbind(b,bb);
###      {b,bb,t,t,t,t}=reg(x,betaw);
      res <- rbind(res, cbind(b,bb));
      if(Eprt>0){
        vrs <- as.matrix(c("ESTDepV", "       ", "coeffs", "se's"))
        print(t(vrs))
        vrs <- as.matrix(c(const, x, const, x))
        a <- cbind(vrs,res)
        vrs <- as.matrix(c("betaB", "   ", "betaW", "   "));
        a <- cbind(vrs,a);
        b <- matrix(1,nrow=4,ncol=1);
        mask <- cbind(matrix(0,nrow=4,ncol=2), matrix(1,nrow=4,ncol=2));
        fmt <- matrix(c("-*.*s ", 8, 8,
                        "-*.*s ", 8, 8,
                        "*.*lf", 7, 4,
                        "*.*lf", 7, 4), nrow=4, ncol=3)
   ##     call printfm(a,mask,fmt);		   
      }
    }
    } else if(identical(tolower(str),"csbetab")){###		@ CI-based sd(betaB)  @
      stbetabs <- eiread(dbuf,"stbetabs");
      a <- stbetaBs[,floor(cols(stbetabs)*0.3413)] ### @ 34th percentile @
      b <- stbetaBs[,floor(cols(stbetabs)*0.6827)]### @ 68th percentile @
      res <- (b-a)/2;
    } else if(identical(tolower(str),"csbetaw")){ ###		@ CI-based sd(betaW)  @
      stbetaws <- eiread(dbuf,"stbetaws");
      a <- stbetaWs[,floor(cols(stbetaws)*0.3413)] ### @ 34th percentile @
      b <- stbetaWs[,floor(cols(stbetaws)*0.6827)] ### @ 68th percentile @
      res <- (b-a)/2;
    } else if(identical(tolower(str),"gebw")){ ###                 @ betaB~betaB for sims betab>=betaw @
      betaBs <- eiread(dbuf,"betabs");
      betaWs <- eiread(dbuf,"betaws");
      a <- betaBs< betaWs;
      betaBs <- mkmissm(betaBs,a);
      betaWs <- mkmissm(betaWs,a);
      res <- cbind(meanwc(t(betaBs),1),meanwc(t(betaWs),1), colSums(1-t(a)));
    }  else if(identical(tolower(str),"gebwa")){ ###                @ B^b ~ B^w for sims betaB >= betaW @
      a <- eiread(dbuf,"gebw");
      res <- cbind(meanwc(a[,1],a[,3])~meanwc(a[,2],a[,3]));
    } else if(identical(tolower(str),"gewb")){ ###  @ betaB~betaW for sims betaW>=betaB @
      betaBs <- eiread(dbuf,"betabs");
      betaWs <- eiread(dbuf,"betaws");
      a <- betaBs> betaWs;
      betaBs <- mkmissm(betaBs,a);
      betaWs <- mkmissm(betaWs,a);
      res <- cbind(meanwc(t(betaBs),1), meanwc(t(betaWs),1), colSums(1-t(a)));
    } else if(identical(tolower(str),"gewba")){###                 @ B^b ~ B^w for sims betaW >= betaB @
      a <- eiread(dbuf,"gewb");
      res <- cbind(meanwc(a[,1],a[,3]), meanwc(a[,2],a[,3]));
    } else if (identical(tolower(str),"bounds")){ ###		@ compute precinct bounds  @
      t <- vread(dbuf,"t");
      x <- vread(dbuf,"x");
      n <- eiread(dbuf,"n");
      lst <- bounds1(t,x,n);
      res <- lst$bs
      a <- lst$aggs
      
    } else if(identical(tolower(str),"bounds2")){	###	@ compute precinct bounds for ei2 @
      res <- NA
      if(vin(dbuf,"undert")){
         v <- vread(dbuf,"t");
         x <- vread(dbuf,"underx");
         n <- eiread(dbuf,"undern");
         t <- eiread(dbuf,"undert");
         lst <- bounds2(v,t,x,n);
         res <- lst[[1]]
         a <- lst$aggs
       }
    } else if(identical(tolower(str),"abounds")) { ###		@ compute aggregate bounds  @
      t <- vread(dbuf,"t");
      x <- vread(dbuf,"x");
      n <- eiread(dbuf,"n");
      lst <- bounds1(t,x,n);
      a <- lst$aggs
      res <- lst$bs
   
    if(Eprt>0){
      vrs <- as.matrix(c("       ", betaBb, etaW, "Aggregate bounds"));
      print(t(vrs));
      message("Lower: ", res[,1])
      message("Upper: ", res[,2]);
    }
    res <- t(res);
      
    } else if (identical(tolower(str),"abounds2")){###		@ compute aggregate bounds for ei2  @
      res <- NA
    if(vin(dbuf,"undert")){
      v <- vread(dbuf,"t");
      x <- vread(dbuf,"underx");
      n <- eiread(dbuf,"undern");
      t <- eiread(dbuf,"undert");
      lst <- bounds2(v,t,x,n);
      a <- lst[[1]]
      res <- lst$aggs
      
    }
     
    if(Eprt>0){
      vrs <- as.matrix(c("       ", lambdaB, lambdaW, "Aggregate bounds"))
      print(t(vrs))
      message("Lower: ", res[,1])
      message(" Upper: ", res[,2]);
    }
    res <- t(res);
    } else if(identical(tolower(str), "pphi")){ ###			@ prints phi and se's  @
      b <- eiread(dbuf,"phi");
      res <- NA
     if(!scalmiss(b)){
       e <- eiread(dbuf,"etas");
      if (ei.vc[eiread(dbuf,"ghactual"),1]!=-1)
        a <- rbind(sqrt(extract.diag(vread(dbuf,"vcphi"))),e)
      else
        a <- matrix(NA, nrow=rows(b),ncol=1);
     
      res <- rbind(t(b), t(a));
      if(Eprt>0){
        message("Maximum likelihood results in scale of estimation (and se's)")
         cat("? ")
        if( vin(dbuf,"parnames")){
          a <- vread(dbuf,"parnames")
          print(t(a))
        }
	 
	 print(res);
      }
     }
    }else if(identical(tolower(str),"psiu")){ ###			@ untruncated psi  @
      
      if(!vin(dbuf,"phi"))
        return(NA);
      b <- eiread(dbuf,"phi");
   ### clearg _Erho;
      Erho <- eiread(dbuf,"Erho");
      assign("Erho", Erho, env=evbase)  
      if (Erho[1]==0){
        c <- eiread(dbuf,"parnames");
        ss <- unlist(sapply(c,identical,"Rho"))
        b <- subset(bx, subset=!ss)
      }
      Zb <- vread(dbuf,"Zb");
      Zw <- vread(dbuf,"Zw");
      x <- vread(dbuf,"x");
      eiread(dbuf,"Ez");
      lst <- eirepar(b,Zb,Zw,x);
      bb <- lst$Bb
      bw <- lst$Bw
      sb <- lst$sb
      sw <- lst$sw
      rho <- lst$rho
      res <- rbind(bb,bw,sb,sw,rho);
    if(Eprt>0){
      message("Untruncated psi's");
      a <- rbind(colMeans(bb),colMeans(bw),colMeans(sb),colMeans(sw),colMeans(rho));
      vrs <- cbind(bb, bw, sb, sw, rho);
      print(t(vrs));
      print(t(a));
    }
      
    } else if(identical(tolower(str), "mppsiu")){###		@ Mean Posterior untruncated psi  @
      b <- eiread(dbuf,"phisims");
      res <- NA
    if(!scalmiss(b)){
      if (cols(b)==2)##  @ i.e., if _EisChk @
        b <- b[,1]
      else
        b <- colMeans(b);
     
     ### clearg _Erho;
      Erho <- eiread(dbuf,"Erho");
      assign("Erho", Erho, env=evbase)
      if(Erho[1]==0){
        c <- eiread(dbuf,"parnames");
        e <- unlist(sapply(c, identical, "Rho"))
        b <- subset(b, subset=!e)
      }
                                   
      zb <- eiread(dbuf,"Zb");
      Zw <- eiread(dbuf,"Zw");
      x <- vread(dbuf,"x");
      eiread(dbuf,"Ez");
      lst <- eirepar(b,Zb,Zw,x);   
      bb <- lst$Bb
      bw <- lst$Bw
      sb <- lst$sb
      sw <- lst$sw
      rho <- lst$rho
      res <- rbind(bb,bw,sb,sw,rho);
      if(Eprt>0){
        message("Mean Posterior Untruncated psi's")
        a <- rbind(colMeans(bb),colMeans(bw),sb,sw,rho);
        vrs <- rbind(bb, bw, sb, sw, rho);
        print(t(vrs));
        print(t(a))
      }
    }
    }else if(identical(tolower(str), "psi")){###			@ ultimate truncated psi  @
      b <- eiread(dbuf,"phi");
      res <- NA
      if(!scalmiss(b)){
      ###     clearg _Erho;
      Erho <- eiread(dbuf,"_Erho");
      assign("Erho", Erho, env=evbase)
      if(Erho[1]==0){
        c <- eiread(dbuf,"parnames");
        e <- unlist(sapply(c, identical, "Rho"))
        b <- subset(b, subset=!e)
      }
      Zb <- eiread(dbuf,"Zb");
      Zw <- eiread(dbuf,"Zw");
      x <- vread(dbuf,"x");
      eiread(dbuf,"Ez");
      res <- eirepar(b,Zb,Zw,x);
      bb <- res$Bb
      bw <- res$Bw
      sb <- res$sb
      sw <- res$sw
      rho <- res$rho
      if(Eprt>0){
        vrs <- rbind(bb, bw, sb, sw, rho);
        message("Truncated psi's (ultimate scale)")
        print(t(vrs));
        print(t(res))
      }
    }
    }else if(identical(tolower(str), "aggs")){###			@ sims x 2 of wtd mean of (betaB&W)  @
      if(vin(dbuf,"x2")){
        nb <- eiread(dbuf,"nb2");
        nw <- eiread(dbuf,"nw2");
      }else{
        nb <- eiread(dbuf,"nb");
        nw <- eiread(dbuf,"nw");
      }
    res <- cbind(meanwc(vread(dbuf,"betaBs"),nb), 
                 meanwc(eiread(dbuf,"betaWs"),nw))
    }else if(identical(tolower(str),"maggs")){###		        @ 2 x 1: meanc(aggBs~aggWs) @
      a <- eiread(dbuf,"aggs");
      res <- colMeans(a);
    } else if(identical(tolower(str),"vcaggs")){###		@ 2 x 2: vcx(aggBs~aggWs)  @
      a <- eiread(dbuf,"aggs");
      res <- var(a);
    } else if(identical(tolower(str), "paggs")){##			@ 2 x 2: ests (se's)  @
      a <- eiread(dbuf,"Maggs");
      b <- eiread(dbuf,"VCaggs");
      b <- sqrt(extract.diag(b));
      res <- rbind(t(a), t(b));
      if(Eprt>0){
      message("Estimates of Aggregate Quantities of Interest")
      vrs <- rbind(betab, betaw)
      print(t(vrs))
      print(res)
    }
    } else if(identical(tolower(str),"goodman")){###		@ 2 x 2 of Goodman's ests (se's)  @
  ###  clearg _Routput,_Rconst;
      assign("Routput", 0,env=evbase)
      assign("Rconst", 0,env=evbase)
      x <- vread(dbuf,"x");
      omx <- 1-x;
      t <- vread(dbuf,"t");
   ### {res,a,tt,tt,tt,tt}=reg(x~omx,t);
      res <- rbind(t(res), t(a));
      if(Eprt>0){
        message("Goodman's Regression")
        vrs <- rbind(betab, betaw);
        print(t(vrs));;
        print(res)
      }
    } else if(identical(tolower(str), "double")){###		@ 2 x 1 of Double regression ests @
### clearg _Routput,_Rconst;
      assign("Routput", 0,env=evbase)
      assign("Rconst", 0,env=evbase)
      res <- NA
    if(vin(dbuf,"undert")){
      t <- vread(dbuf,"undert");
      x <- vread(dbuf,"underx");
      v <- eiread(dbuf,"t");
      omx <- cbind(x,(1-x));
 ###     {a,tt,tt,tt,tt,tt}=reg(omx,t);
  ###    {b,tt,tt,tt,tt,tt}=reg(omx,v.*t);
      res <- b/a;
      if(Eprt>0){
        message("Double Regression")
        vrs <- rbind(lambdaB, lambdaW)
        print(t(vrs));
        print(t(res));
      }
    }
     
    } else if(identical(tolower(str), "neighbor")){ ###             @ 2 x 1 of Neighborhood estimates @
      t <- eiread(dbuf,"t");
      b <- eiread(dbuf,"Nb");
      a <- eiread(dbuf,"Nw");
      res <- rbind(meanwc(t,b),meanwc(t,a));
      if(Eprt>0){
        message("Freedman et al.'s Neighborhood Model Estimates");
        vrs <- rbind(betab, betaw);
        print(t(vrs));
        print(t(res));
    }
    } else if(identical(tolower(str),"thomsen")){###		@ 2 x 1 of Thomsen's estimates @
      x <- vread(dbuf,"x");
      t <- vread(dbuf,"t");
      a=(x==0)| (x==1) | (t ==0) | (t==1);
      invX <- qnorm(subset(x, subset=a));
      invT <- qnorm(subset(t, subset=a));
      meanX <- -colMeans(invX);
      meanT <- -colMeans(invT);
      rho <- cor(cbind(invX,invT));
      rho <- rho[1,2];
      p00 <- cdfbvn(meanX,meanT,rho);
      
    p10 <- pnorm(meanT)-p00;
    p01 <- pnorm(meanX)-p00;
    p11 <- 1-p00-p10-p01;
    bb <- p11/(p11+p10);
    bw <- p01/(p01+p00);
    res <- rbind(bb,bw)
    if(Eprt>0){
      message("Thomsen's Ecological Logit Approach Estimates");
      vrs <- cbind(betab,betaw);
      print(t(vrs));
      print(t(res));
    }
    } else if(identical(tolower(str),"palmquist")){###		@ Palmquist's inflation factor  @
      x <- vread(dbuf,"x");
      n <- eiread(dbuf,"n");
      b <- meanwc(x,n);
      a <- (meanwc(x^2,n)-b^2)/(b*(1-b));
      res <- (1/a)-1;
      if (Eprt>0)
        message("Palmquist's Inflation Factor: ", res)
    }else if(identical(tolower(str),"truthb")){###		@ true betab  @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first")
      return(NA)
    }
    truth <- eiread(dbuf,"truth");
    res <- truth[,1];
  }  else if(identical(tolower(str),"truthw")){###		@ true betaw @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first");
    
      return(NA);
    }
    truth <- eiread(dbuf,"truth");
    res <- truth[,2];
  }  else if(identical(tolower(str),"nbv") | identical(tolower(str),"nbt")){###	@ number of blacks who Turnout @
    b <- eiread(dbuf,"truthb");
    nb <- eiread(dbuf,"nb");
    res <- nb*b;
  } else if(identical(tolower(str),"nbn")){###			@ number of blacks who don't vote @
    b <- eiread(dbuf,"truthb");
    nb <- eiread(dbuf,"nb");
    res <- nb*(1-b);
  } else if(identical(tolower(str),"nwv") | identical(tolower(str),"nwt")){###	@ number of whites who Turnout @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first");
      return(NA);
    }
    b <- eiread(dbuf,"truthw");
    nw <- eiread(dbuf,"nw");
    res <- nw*b;
  } else if(identical(tolower(str),"nwn")){###			@ number of blacks who don't vote @
    if(!vin(dbuf,"truth")){
      message("eiread: truth needs to be stored first");
      return(NA)
    }
     
    b <- eiread(dbuf,"truthw");
    nw <- eiread(dbuf,"nw");
    res <- nw*(1-b);
  } else if(identical(tolower(str),"aggtruth")){###		@ aggregate truths  @
    if(!("truth   " %inG% cv)){
      message("eiread: truth needs to be stored first")
      return(NA);
    }
   
    n <- eiread(dbuf,"n");
    x <- vread(dbuf,"x");
    nb <- x*n;
    nw <- (1-x)*n;
    res <- vread(dbuf,"truth");
    res <- rbind(meanwc(res[,1],nb), meanwc(res[,2],nw));
    if(Eprt>0){
      message("Aggregate Truth");
      vrs <- cbind(betab, betaw)
      print(t(vrs));
      print(t(res));
    }
  } else if(identical(tolower(str),"psitruth")){###		@ true psi's, truncated scale  @
    a <- eiread(dbuf,"truth");
    res <- NA
    if(!scalmiss(a)){
      b <- cor(na.omit(a));
      res <- rbind(meanc(packr(a)), sd(as.data.frame(na.omit(a))),b[2,1])
      if (Eprt>0){
        message("TRUE truncated psi's (ultimate scale)")
        vrs=cbind(bb, bw, sb, sw, rho);
        print(t(vrs));
        print(t(res));
      }
    }
    }else if(identical(tolower(str),"tsims")){###			@ sims from p(T|X=seqas(0,1,100))  @
      c <- eiread(dbuf,"Eeta");
    if(scalzero(c) && 
      (!(scalone(eiread(dbuf,"Zb"))) || !(scalone(eiread(dbuf,"Zw"))))){
      message("eiread: tsims only works without covariates.")
      return(NA);
    }
      a <- 100;
      x <- seqase(0,1,a);
      b <- eiread(dbuf,"phi");
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
      eiread(dbuf,"Ez");
      zb <- eiread(dbuf,"zb");
      zw <- eiread(dbuf,"zw");
      lst <- eirepar(b,0,0,0);
      bb0 <- lst$Bb
      bw0 <- lst$Bw
      sb <- lst$sb
      sw <- lst$sw
      rho <- lst$rho
      res <- matrix(0,nrow=a,ncol=Esims);
      bnds <- matrix(c(0,1,0,1),nrow=2, byrow=T)
      e <- meanc(vread(dbuf,"x"));
      for(i in 1:a){
        bb <- bb0+b[rows(b)-1]*(x[i+0]-e);
        bw <- bw0+b[rows(b)]*(x[i+0]-e);
        bs <- rndbtn(bb,bw,sb,sw,rho,bnds,Esims);
        t <- bs[,1]*x[i+0]+bs[,2]*(1-x[i+0]);
        res[i+0,] <- t(sortc(t,1));
      }
      res <- cbind(x, res)
    }

    } else if( identical(tolower(str), "tsims0")){ ###          @ sims from p(T|X_obs,Z_obs) @
      b <- eiread(dbuf,"phi");
      if(scalmiss(b))
        res <- NA
      
      eiread(dbuf,"_Ez");
      t <- eiread(dbuf,"t");
      x <- eiread(dbuf,"x");
      zb <- eiread(dbuf,"zb");
      zw <- eiread(dbuf,"zw");
      a <- ifelse(is.matrix(x), rows(x), length(x));
      
      lst <- eirepar(b,zb,zw,x);
      bb <- lst$Bb
      bw <- lst$Bw
      sb <- lst$sb
      sw <- lst$sw
      rho <- lst$rho
      c <- eiread(dbuf,"Esims");
      res <- matrix(0, nrow=a,ncol=c);
      bnds <- matrix(c(0, 1, 0, 1),nrow=2, ncol=2, byrow=T)
      for (i in 1:a){
        bs <- rndbtn(bb[i],bw[i],sb,sw,rho,bnds,c);
        p <- bs[,1]*x[i+0]+bs[,2]*(1-x[i+0]);
        res[i+0,] <- t(sortc(p,1));
      }
      res <- cbind(t,res);
      
    }else if(identical(tolower(str),"expvarci")){ ###		@ x~20%CI~mean~80%ci @
      if( rows(eiread(dbuf,"Eeta"))==4 &&
      (!scalone(eiread(dbuf,"Zb")) || !scalone(eiread(dbuf,"Zw")))){
      message("eiread: expvarci only works without covariates.");
      return(NA)
    }
      b <- eiread(dbuf,"tsims");
      res <- NA
    if(!scalmiss(b)){
      x <- b[,1];
      b <- b[,2:(Esims+1)];
      res <- cbind(x,b[,floor(0.2*Esims)], colMeans(t(b)), b[,floor(0.8*Esims)]);
    }
    } else if(identical(tolower(str), "expvarci0")){###		@ t~20%CI~mean~80%ci @
      b <- eiread(dbuf,"tsims0");
      c <- eiread(dbuf,"Esims");
      res <- NA
      if (!scalmiss(b)){
        t <- b[,1];
        b <- b[,2:c+1];
        res <- cbind(t, b[,floor(0.2*c)], colMeans(t(b)), b[,floor(0.8*c)]);
      }
    
  }else if(identical(tolower(str), "expvarcis")){###		@ x~20%CI~mean~80%ci LOESS smoothed @
    b <- eiread(dbuf,"expvarci");
    res <- NA
    if(!scalmiss(b)){
      e <- output;
      output <- 0;
      loess.WgtType <- 2;
      loess.span <- 0.45;
      y.loess <- loess(b[,2] ~ b[,1], b, weights=loess.WgtType, span= loess.span,degree=1)
      tt <- y.loess$fitted
      c <- y.loess$y
    ##  a <- y.loess$x
      res <- cbind(c,tt);
      y.loess <- loess(b[,3] ~ b[,1], b, weights=loess.WgtType, span= loess.span,degree=1)
      tt <- y.loess$fitted
      c <- y.loess$y
   ###   a <- y.loess$x
      res <- cbind(res,tt);
      y.loess <- loess(b[,4] ~ b[,1], b, weights=loess.WgtType, span= loess.span,degree=1)
      tt <- y.loess$fitted
      c <- y.loess$y
    ###  a <- y.loess$x
      res <- cbind(res,tt);
      
      output <- e;
    }

  } else if(identical(tolower(str),"sum")){ ###			@ prints all printable items @
    if(Eprt<1)
      Eprt <- 1
   
    if("titl    " %inG% cv){
      titl <- vread(dbuf,"titl");
      print(paste("** ",titl," **", sep=""));
      if(identical(titl,"*DB* Data Buffer from eimodels_avg() *DB*")){
        postp <- vread(dbuf,"postprob");
        priorp <- vread(dbuf,"prprob");
        mrgllk <- vread(dbuf,"margllik");
       ## ?; 
        message("The number of model averaged:", rows(vread(dbuf,"postprob")));
       ## ?;
        message("_EI_bma_est: ", ei.bma.est);
       ## ?;
        message("Model  Posterior Prior  Marginal");
        message("Number   Prob     Prob   LogLik ", cbind(postp, priorp[,2], mrgllk[,2]));
       ## ?;
        eiread(dbuf,"abounds");
        ##?;
        eiread(dbuf,"paggs");
        res <- "";
        return(res);
      }
    }
  
    if (scalone(eiread(dbuf,"Enonpar"))){
      message("Nonparametric Estimation");
      message("EnonNumInt:   ", eiread(dbuf,"EnonNum"));
      message("EnonEval:     ", eiread(dbuf,"_EnonEva"));
      message("N:             ", rows(vread(dbuf,"x")));
      message("Esims:        ", vread(dbuf,"_Esims"));
     ## ?;
    }else{
    
      message("CML return: ", vread(dbuf,"retcode"),"     ")
      message("N:          ", rows(vread(dbuf,"x")), "      ")
      message("Esims:     ", vread(dbuf,"Esims"))
      if ("ebeta  " %inG% cv)
        message("Ebeta      ", vread(dbuf,"Ebeta"),"     ")
     
      message("Esigma:    ", vread(dbuf,"Esigma"),"       ")
      message("Erho:      ", t(vread(dbuf,"Erho")))
      message("Eisn:      ", vread(dbuf,"Eisn"),"     ")
      message("resamp:     ",eiread(dbuf,"resamp"))
      message("GhActual:  ",vread(dbuf,"ghactual"));
      message("Estval:    ",vread(dbuf,"_Estval"));
      if("eeta   " %inG% cv)
        message("Eeta:      ", eiread(dbuf,"_Eeta"))
      
      message("log-likelihood:         ", vread(dbuf,"loglik"));
      message("ln(mean(Imptce Ratio)): ", eiread(dbuf,"meanIR"));
     ## ?;
      eiread(dbuf,"pphi");
     ## ?;
      eiread(dbuf,"psiu");
      ##?;
      eiread(dbuf,"psi");
      ##?;
    
    }
      if ("truth   "%inG% cv)
        {
          eiread(dbuf,"psitruth");
        ##?;
          eiread(dbuf,"aggbias");
          ##?;
          eiread(dbuf,"eaggbias");
          ##?;
          eiread(dbuf,"coverage");
###?; 
          eiread(dbuf,"aggtruth");
###?;
        }
    eiread(dbuf,"abounds")
    ###?;
    eiread(dbuf,"paggs")
    res <- ""
    
  }else{
    if (Eprt>0)
      message("eiread: no such name, ",str);
    
    res <- NA;
  
  }
  
  
  return(res);
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
checkr <- function(dbuf,eps){
###  local phi,R,rr,zb,zw,x,y,k,kk,loparms,hiparms,loR,hiR,rchk;

  phi <- eiread(dbuf,"phi");
  rr <- rows(phi);
  lst = pluckdta(eiread(dbuf,"dataset"));
  Zb <- lst$Zb
  Zw <- lst$Zw
  x <- lst$x
  y <- lst$y
  loR <- matrix(1,nrow=rr-2,ncol=1);
  hiR <- matrix(1, nrow=rr-2,ncol=1);
  loparms <- phi-eps;
  hiparms <- phi+eps;
  R <- colSums(na.omit(lncdfbvnu(eirepar(phi,zb,zw,x))));
  for(kk in 1:(rr-2)){
    k <- kk;
    loR[k] <- colSums(na.omit(lncdfbvnu(eirepar(loparms,zb,zw,x))));
    hiR[k] <- colSums(na.omit(lncdfbvnu(eirepar(hiparms,zb,zw,x))));
  }

  rchk=(cbind(loR,hiR) !=R);

return(rchk);
}
