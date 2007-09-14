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
         if(nrow(res) == 1)
           res <- matrix(0, nrow=2, ncol=1)
         else if(nrow(res) == 3 && res[1]==4)
           res <- as.matric(c(0, res[2]))
         else if(nrow(res) ==3 && res[1] == 5)
           res <- as.matrix(c(res[2], 0))
         else if(nrow(res)==4)
           res <- as.matrix(res[1:2])
         
                       
       }else if(identical(tolower(str),"etas")){
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
       }else if(identical(tolower(str), "ez")){
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
          
          res <- minindc(abs(t(stbetabs-truth[,1])))/b;    
          res <- cbind(res, minindc(abs(t(stbetaws-truth[,2])))/b);
          res[,1] <- recode(res[,1],stdc(t(stbetabs)) <= EnumTol,0.5); ###@ homog prects @
          res[,2] <- recode(res[,2],stdc(t(stbetaws))<= EnumTol,0.5);
          
        }else if(identical(tolower(str),"ci50b")){### @ 50% confidence intervals for betab @
          stbetabs <- eiread(dbuf,"stbetabs");
          e <- ncol(stbetabs);
          res <- st
        }else if(identical(tolower(str),"ci80b"))	{ ###@ 80% confidence intervals for betab @
          stbetabs <- eiread(dbuf,"stbetabs");
          e <- ncol(stbetabs);
          res <- cbind(stbetabs[,floor(0.1*e)], stbetabs[,floor(0.9*e)]);

        } else if(identical(tolower(str), "ci95b")){ ###@ 95% confidence intervals for betab @
          stbetabs <- eiread(dbuf,"stbetabs");
          e <- ncol(stbetabs);
          res <- cbind(stbetabs[,floor(0.05*e)],stbetabs[,floor(0.95*e)]);
        } else if(identical(tolower(str),"ci50w")){###	@ 50% confidence intervals for betaw @
          stbetaws <- eiread(dbuf,"stbetaws");
          e <- ncol(stbetaws);
          res <- cbind(stbetaws[,floor(0.25*e)], stbetaws[,floor(0.75*e)]);
        } else if(identical(tolower(str),"ci80w")){###		@ 80% confidence intervals for betaw @
          stbetaws <- eiread(dbuf,"stbetaws");
          e <- ncol(stbetaws);
          res <- cbind(stbetaws[,floor(0.1*e)],stbetaws[,floor(0.9*e)]);
        } else if(identical(tolower(str),"ci95w")){ ###	@ 95% confidence intervals for betaw @
          stbetaws <- eiread(dbuf,"stbetaws");
          e <- ncol(stbetaws);
          res <- cbind(stbetaws[,floor(0.05*e)],stbetaws[,floor(0.95*e)]);
        }  else if(identical(tolower(str),"ci80bw")){ ###		@ 80% conf intervals for betab betaw  @
          res <- cbind(eiread(dbuf,"ci80b"), eiread(dbuf,"ci80w"));
        } else if(identical(tolower(str),"ci95bw")){###		@ 95% conf intervals for betab betaw  @
          res <- cbind(eiread(dbuf,"ci95b"),eiread(dbuf,"ci95w"));
####################################
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
      {b,bb,t,t,t,t}=reg(x,truth[.,1]);
      res <- cbind(b,bb)
      {b,bb,t,t,t,t}=reg(x,truth[.,2]);
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
	let fmt[4,3]=
      fmt <- matrix(c("-*.*s ", 8, 8,
                     "-*.*s ", 8, 8,
                     "*.*lf", 7, 4,
                     "*.*lf", 7, 4), nrow=4, ncol=3)
        call printfm(a,mask,fmt);		  
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
      
      {b,bb,t,t,t,t}=reg(x,betab);
      res <- cbind(b,bb);
      {b,bb,t,t,t,t}=reg(x,betaw);
      res <- rbind(res, cbind((b,bb)));
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
        call printfm(a,mask,fmt);		   
      }
    }
    } else if(identical(tolower(str),"csbetab")){###		@ CI-based sd(betaB)  @
      stbetabs <- eiread(dbuf,"stbetabs");
      a <- stbetaBs[,floor(ncol(stbetabs)*0.3413)] ### @ 34th percentile @
      b <- stbetaBs[,floor(ncol(stbetabs)*0.6827)]### @ 68th percentile @
      res <- (b-a)/2;
    } else if(identical(tolower(str),"csbetaw")){ ###		@ CI-based sd(betaW)  @
      stbetaws <- eiread(dbuf,"stbetaws");
      a <- stbetaWs[,floor(ncol(stbetaws)*0.3413)] ### @ 34th percentile @
      b <- stbetaWs[,floor(ncol(stbetaws)*0.6827)] ### @ 68th percentile @
      res <- (b-a)/2;
    } else if(identical(tolower(str),"gebw")){ ###                 @ betaB~betaB for sims betab>=betaw @
      betaBs <- eiread(dbuf,"betabs");
      betaWs <- eiread(dbuf,"betaws");
      a <- betaBs< betaWs;
      betaBs <- mkmissm(betaBs,a);
      betaWs <- mkmissm(betaWs,a);
      res <- cbind(meanwc(t(betaBs),1),meanwc(t(betaWs),1), colSums(1-t(a)));
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
checkr <- function(dbuf,eps){
###  local phi,R,rr,zb,zw,x,y,k,kk,loparms,hiparms,loR,hiR,rchk;

  phi <- eiread(dbuf,"phi");
  rr <- nrow(phi);
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
