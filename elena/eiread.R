library(mvtnorm)

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
        a <- matrix(NA, nrow=nrow(b),ncol=1);
     
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
      if (ncol(b)==2)##  @ i.e., if _EisChk @
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
      res <- eirepart(b,Zb,Zw,x);
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
        b[nrow(b)-1] <- b[2]
      else if( c==2 || c==5)
        b[nrow(b)] <- b[3]
      else if( c==3 && nrow(c)==1){
        b[nrow(b)-1] <- b[2]
        b[nrow(b)] <- b[4]
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
        bb <- bb0+b[nrow(b)-1]*(x[i+0]-e);
        bw <- bw0+b[nrow(b)]*(x[i+0]-e);
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
      a <- ifelse(is.matrix(x), nrow(x), length(x));
      
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
      if( nrow(eiread(dbuf,"Eeta"))==4 &&
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
      a <- y.loess$x
      res <- cbind(a,c);
      y.loess <- loess(b[,3] ~ b[,1], b, weights=loess.WgtType, span= loess.span,degree=1)
      tt <- y.loess$fitted
      c <- y.loess$y
      a <- y.loess$x
      res <- cbind(res,c);
      y.loess <- loess(b[,4] ~ b[,1], b, weights=loess.WgtType, span= loess.span,degree=1)
      tt <- y.loess$fitted
      c <- y.loess$y
      a <- y.loess$x
      res <- cbind(res,c);
      
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
        message("The number of model averaged:", nrow(vread(dbuf,"postprob")));
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
      message("N:             ", nrow(vread(dbuf,"x")));
      message("Esims:        ", vread(dbuf,"_Esims"));
     ## ?;
    }else{
    
      message("CML return: ", vread(dbuf,"retcode"),"     ")
      message("N:          ", nrow(vread(dbuf,"x")), "      ")
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
