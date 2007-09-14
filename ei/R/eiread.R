eiread <- function(dbuf, str, compute =FALSE){

  if(!compute)
   return(vread(dbuf,str))
  evbase  <- get("evbase", env=parent.frame())
  getEnvVar(evbase, environment())###, vecvar=c("eimetar"))
  eimetar <- get("eimetar", env=evbase)
 ### the computation part needs work and is not implemented yet   
   if(vin(dbuf, "titl")){
      titl <- vread(dbuf, "titl")
      nc  <- nchar(as.character(floor(eimetar))) 
      fmt <- formatC(eimetar,width=nc,digits=0, format="f")
     if (identical(titl, "*MDB* Meta-Data Buffer from 2nd Stage *MDB*")){ 
       dbuf <- vread(dbuf,paste("dbuf",fmt,sep=""))
     }else if(identical(titl, "*MDB* Meta-Data Buffer from eimodels_def() *MDB*")){
       if(vin(dbuf,paste("mod.d",fmt, sep=""))){
         dbuf <- vread(dbuf,paste("mod.d",fmt,sep=""))
         if(Eptr==3)
           message("Reading ", str, " from Model ", eimetar, "...");
       }else{
         message("eiread: Model " + eimetar + " is not stored in this data buffer.");
         res <- NA
         return(res)
       }
              
      }else if(identical(titl,"*MDB* Meta-Data Buffer from eimodels_run() *MDB*")){
        if(vin(dbuf, paste("mod.r"+fmt,sep=""))){
          dbuf <- vread(dbuf, paste("mod.r",fmt, sep=""))
          if(Eptr == 3)
            message("Reading ", str, " from Model ", eimetar,"...")
        }else{
          mess <- paste("eiread: Model", eimetar, "is not stored in this data buffer")
          message(mess)
          res <- NA
          return(res)
        }
      }
    }
  strtmp <- str
  str <- tolower(str)
###  cv <- vnamecv(dbuf)
  cv <- names(dbuf)
  strc <- str
  if(vin(dbuf, "titl")){
    
    titl <- vread(dbuf, "titl")
    if(identical(titl, "*DB* Data Buffer from eimodels_avg() *DB*")){
      vrs <- c("ealphab","ealphaw","ebeta","ebounds","ecdfbvn","edirtol","edoml", 
               "edoml.phi","edoml.vcphi","eeta","ei.vc","eigraph.bvsmth","eischk","eisfac",
               "eisn","eist","emaxiter","enoneval","enonnumint","enonpar","enumtol","erho",
               "eselrnd", "esigma","estval","evtol","ei2.m","eimetar","eimodels.save","zb", "zw",
               "t","x","ez"
###  ,"x2","x2rn", "checkr", "dataset", "etac", "etas", "expvarci","expvarcis","lnir"
### ,"loglik", "logliks", "meanir", "mppsiu", "parnames","phi","phisims", "pphi"
### ,"psi", "psitruth", "psiu", "r",  "ri", "retcode", "tsims","vcphi"
     );
     vrs <- as.matrix(tolower(vrs))
     if(identical(str, "esims")){
       res <- ncol(eiread(dbuf,"betabs"))
       return(res)
     }
     if((str %inG% vrs) && !identical(str, "expvarci0")){
       mess <- paste("eiread: You cannot eiread", strtmp, "from the output buffer of eimodels.avg")
       message(mess)
       return(NA)
     }
    }
  }
   
     
   if (vin(dbuf,"titl")){
    if(identical(titl,"*MDB* Meta-Data Buffer from eimodels_def() *MDB*")){
      vrs <- c("under.t", "or",  "under.x", "under.ez", "abounds", "abounds2", "aggbias", "beta", "betab", "betabs",
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
  }
 ### changes in stored globals
       if(identical(tolower(str), "eeta")){

         if(vin(dbuf,"Eeta"))
           res <- vread(dbuf, "eeta")
         else
           
           res <- matrix(0, nrow=4, ncol=1)
         
         if(scalmiss(res))
           res <- matrix(0, nrow=4, ncol=1)
         else if(ncol(res) == 2)
           res <- as.matrix(as.vector(res))
         else if(nrow(res)==2)
           res <- rbind(res, matrix0, nrow=2, ncol=1)
       }else if(identical(tolower(str), "zb")){
         e <- eiread(dbuf, "eeta")
         res <- vread(dbuf, "Zb")
         if(e %in% 3:4 || e == 1)
           res <- vread(dbuf, "x")
         
       }else if(identical(tolower(str), "zw")){
         e <- eiread(dbuf, "eeta")
         res <- vread(dbuf, "ZW")
         if(e %in% 2:3 || e ==5)
           res <- vread(dbuf,"x")
         
       }else if(identical(tolower(str), "titl") || identical(tolower(str), "title")){
          res <- ""
         if(vin(dbuf,"titl"))
           res <- vread(dbuf, "titl")
     ###      /***** stored globals *****/
        }else if(strc %inG% cv){
          res <- vread(dbuf,str)
          ### /***** computed results *****/
       }else if(identical(tolower(str), "x2")){
         if(vin(dbuf, "x2"))
           res <- vread(dbuf,"x2")
         else{
           message("eiread: 'x2' option is only available in data buffers created by ei2")
           res <- NA
         }
     ###    	@ horizontally randomly permuted x2  @
       }else if(identical(tolower(str), "x2rn")){
         if(vin(dbuf, "x2")){
           res <- vread(dbuf, "x2")
           a <- nrow(res)
           c <- ncol(res)
           for(n in 1:a)
             res[n, ] <- res[n, order(runif(c))]
         }else{
           message("eiread: 'x2' option is only available in data buffers created by ei2")
           res <- NA
         }
       }else if(identical(tolower(str), "emaxiter")){
         res <- NA
         if(vin(dbuf, "emaxiter"))
           res <- vread(dbuf, "emaxiter")
         
       }else if(identical(tolower(str),"eigraph.bvsmth") ||identical(tolower(str),"bvsmth")){
         res <- NA
         if(vin(dbuf, "bvsmth"))
           res <- vread(dbuf, "bvsmth") 
       }else if(identical(tolower(str),"eimodels.save")){
         res <- NA
         if(vin(dbuf, "eimsave"))
           res <- vread(dbuf, "eimsave")
       }else if(identical(tolower(str),"ei.bma.est")){
         res <- NA
         if(vin(dbuf, "eibmaest"))
           res <- vread(dbuf, "eibmaest")
       }else if(identical(tolower(str),"ei.bma.prior")){
         res <- NA
         if(vin(dbuf, "prprob")){
           priorp <- vread(dbuf,"prprob")
           res <- priorp[, 2]
         }


       }else if(identical(tolower(str),"edoml.phi") || identical(tolower(str), "doml.phi")){
         res <- NA
         if(vin(dbuf, "doml.phi"))
           res <- vread(dbuf, "doml.phi")

       }else if(identical(tolower(str),"enoneval") || identical(tolower(str), "enoneva")){
         res <- NA
         if(vin(dbuf, "enoneva"))
           res <- vread(dbuf, "enoneva")
       }else if (identical(tolower(str), "meanir")){
         res <- NA
       if(vin(dbuf, "meanir"))
         res <- vread(dbuf, "meanir")
       else if(vin(dbuf, "EmeanIR"))
         res <- vread(dbuf, "emeanir")
       else if(vin(dbuf, "lnir") && eiread(dbuf,"EisChk") == 0)
         res <- vread(dbuf, "lnir")
       else if(vin(dbuf, "lnir") && eiread(dbuf,"EisChk") == 1){
         a <- eiread(dbuf, "lnir")
         a <-  a[,1]
         max <- max(a)
         res <- max + log(meanc(exp(a-max)))
         
       }
         
       }else if (identical(tolower(str), "logliks")){

         a <- eiread(dbuf, "phi")
         b <- eiread(dbuf, "dataset")
      ###   res <- eiloglik(a,b)
       }else if((identical(tolower(str), "resamp")){
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
         if(nrow(res) == 1)
           res <- matrix(0, nrow=2, ncol=1)
         else if(nrow(res) == 3 && res[1]==4)
           res <- as.matric(c(0, res[2]))
         else if(nrow(res) ==3 && res[1] == 5)
           res <- as.matrix(c(res[2], 0))
         else if(nrow(res)==4)
           res <- as.matrix(res[1:2])
         
                       
       }else if(identical(tolower(str,"etas"))){
         res <- eiread(dbuf, "Eeta")
         if(nrow(res) == 1)
           res <- matrix(0, nrow=2, ncol=1)
         else if (nrow(res) == 3 && res[1]==4)
           res <- as.matrix(c(0, res[3]))
         else if (nrow(res) == 3 && res[1]==5)
           res <- as.matrix(c(res[3], 0))
         else if(nrow(res) == 4)
           res <- as.matrix(res[3:4])
       ###  @ n of covariates, incl. implied constant for Zb|Zw @
       }else if(identical(tolower(str), "ez"){
         zb <- eiread(dbuf, "zb")
         zw <- eiread(dbuf,"zw")
         assign("Ez", 0, env=evbase)
         Ez <- as.matrix(c(ncol(Zb) + 1 - (Zb == 1), ncol(Zw) + 1 - (Zw == 1)))
         assign("Ez", Ez, env=evbase)
         res <- Ez

       }else if(identical(tolower(str),"ebounds")){
         dummy <- eiread(dbuf, "Ez")
         b <- eiread(dbuf,"Ebounds")
         if(scalzero(b))
           res <- matrix(c(-20, 20), nrow=1, ncol=2)
         else if(ncol(b) == 2)
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
         message("eiread: problem with _Ebounds")
       }else if(identical(tolower(str), "nobs")){
         res <- nrow(vread(dbuf, "t"))
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
           res <- vread((dbuf, "tvap")
       }else if(identical(tolower(str), "bvap") || identical(tolower(str), "nb")) ###	@ black vap  @{
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
       } else if(identical(tolower(str),"nw2"){	###@ white vap  @
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
          a <- nrow(res)
          c <- ncol(res);
          for( i in 1:a)
            res[i,] <- res[i,order(runif(c))]   
   
        } else if(identical(tolower(str), "rnbetaws")){ ###		@ randomly permuted betabs sims @
          res <- eiread(dbuf,"betaws");
          a <- nrow(res);
          c <- ncol(res);
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
          b <- ncol(stbetabs);
          
          res <- minindc(t(abs(stbetabs-truth[,1])))/b;    
          res <- cbind(res, (minindc(abs(stbetaws-truth[.,2])')/b);
    res[.,1]=recode(res[.,1],stdc(stbetabs').<=_EnumTol,0.5); @ homog prects @
    res[.,2]=recode(res[.,2],stdc(stbetaws').<=_EnumTol,0.5);

      

  
    
    
         
                
         
     
}     
