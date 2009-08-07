### DESCRIPTION takes the match.call of a function and
###             use the variables in expand.dot to assign them
###             to to the global environment evbase.  Use in ei()
###
### INPUT: drvdot=macth.call(expand.dot=T)  and
###        drv=match.call(expand.dot=FALSE)
###        evbase environment
###
### OUTPUT modify evbase with new variables from expanding the dots
###
### AUTHOR Elena Villalon (evillalon@iq.harvard.edu)

expanddots <- function(drvdot, drv, evbase=.GlobalEnv){
 
  args <- names(drv)
  args <- args[-1]
  
  if(args[length(args)]=="...") ### counts for "..."
    args <- args[-length(args)]
  argsdots <- names(drvdot)
  argsdots <- argsdots[-1]

  if(length(argsdots) == length(args))
    return(NULL)  ##nothing in expand.dots
  ln <- length(args)
 
  expn  <- setdiff(argsdots, args)
 
  argsdots <- nm <- names(drvdot)
   
  ix <- sapply(expn, grep, argsdots)
  ix <- unlist(ix)
  
  for( n in 1:length(argsdots)){
    if(!n %in% ix)
      next
    if(argsdots[[n]]=="")
      stop("ei: entries in ... should be name")
      
#    if(length(as.list(drvdot[[n]])) > 1){
#      vec <- as.character(drvdot[[n]])
#      sign <- 1
#      vsgn <- trim.blanks(vec[1])
#      if(length(vec) ==2 && identical(vsgn,"-")) sign <- -1
#     
#      lst <- NULL
# 
#      for(m in 2:length(vec)){
#        print(as.list(drvdot[[n]])[[m]])
#     
#        lst <- c(lst,sign*(as.list(drvdot[[n]])[[m]]))
#      }
#      val <- as.matrix(lst)
#    }else
    #Trying this out! Only symbols were getting assigned!
    frame <- parent.frame()
    nframe <- 1
    #Kludge to crawl back up the frame stack to evaluate arguments.
    while(!identical(frame, .GlobalEnv)){
        val <- try(eval(drvdot[[n]], envir=frame), TRUE)
        if(class(val) != "try-error") break

        nframe <- nframe + 1
        frame <- parent.frame(n=nframe)
    }
    if(class(val) == "try-error")
       stop(val) 
        
    if(!identical(evbase,.GlobalEnv))
      assign(argsdots[n], val, env=evbase)
        
  }
    if(!identical(evbase,.GlobalEnv))
      assign("evbase", evbase,env=parent.frame())  ##not necessary 
 ### return(evbase)
  
  return(argsdots[ix])
 

}
###DESCRIPTION It assigns enviromental variables of env=evbase
###            to env=evto with the same names. If vecvar is nonnull
###            then it only assigns those variables in vecvar
###
### INPUT two environments and, optionally, vector with names of some variables
###
### OUTPUT it assigns variables from evfrom to evto
### AUTHOR Elena Villalon (evillalon@iq.harvard.edu)
###
getEnvVar <- function(evfrom, evto, vecvar=NULL){

  param <- ls(env=evfrom)
### select some variables from param, those in vecvar  
  if(length(vecvar)){
     vecvar <- sapply(vecvar,function(x) { paste("^",x,"$", sep="")})
     ix <- sapply(vecvar, grep, param, ignore.case=T, extended=T)
     ix <- unlist(ix)
     if(length(ix)) param <- param[ix]
   }
  evfrom <- as.environment(evfrom)
  evto <- as.environment(evto)
  
  envlapply <- lapply
  ass <- envlapply(param, function(att,evfrom,evto) {
  
    #val <- try(get(att, env=evfrom))
    val <- get(att, env=evfrom)
  
    #if(class(val) == "try-error")
    #  message("Error getting param")
    assign(att, val,env=evto)}, evfrom,evto)
  
return(evto)
}
###DESCRIPTION Changes the values of globals or
###            assigned them to env=evbase
### INPUT evbase an environmnet with global parameters
###       lst the vecor with values to be assigned,
###       nm the names of the variables for the assigment,
###          if null then those names are part of lst
### OUTPUT the new evbase (environmnet) with the additions

setGlobals <- function(evbase, lst, nm=NULL){

if(!length(nm)){
  if (is.vector(lst) || is.list(lst))
    nm <- names(lst)
  else if( is.matrix(lst)){
    nm <- ifelse(nrow(lst) > 1, dimnames(lst)[[1]],  dimnames(lst)[[2]])
    lst <- as.vector(lst)
    names(lst) <- nm
  }
}
if(!length(nm)){
  warning("ei: setGlobals no param provided for assigment")
  return(evbase)
}
ix <- 1:length(nm)

res <- sapply(ix,function(n,lst,nm, evbase){
  ch  <- nm[[n]]
  val <- lst[[n]]
  assign(ch,val, env=evbase)
})
return(evbase)
}

### A convenient function used in ei() to update some globals
###
add.to.Eres <- function(Eres=list(), round=1, evbase=NULL,...)
{
  if(!length(evbase))
    evbase <- eiset()
  
  evei <- getEnvVar(evbase, environment())
        
 ### inputs
    if(round == 1){
      ### vput(Eres,t,"t") will work as well BUT R CMD check complains 
      if(exists("t"))
        Eres <- vput(Eres,get("t",env=evei),"t") 
      if(exists("x"))
        Eres <- vput(Eres,get("x", env=evei),"x")
      if(exists("tvap"))
        Eres <- vput(Eres,get("tvap", env=evei),"n")
      if(exists("Zb"))
        Eres <- vput(Eres,get("Zb", env=evei),"Zb")
      if(exists("Zw"))
        Eres <- vput(Eres,get("Zw", env=evei),"Zw")
  
## essential globals
      if(exists("EalphaB"))
        Eres <- vput(Eres,get("EalphaB",env=evei),"EalphaB")
     if(exists("EalphaW")) 
       Eres <- vput(Eres,get("EalphaW", env=evei),"EalphaW")
      if(exists("Ebeta"))
        Eres <- vput(Eres,get("Ebeta",env=evei),"Ebeta")
      if(exists("Ebounds"))
        Eres <- vput(Eres,get("Ebounds", env=evei),"Ebounds")
      if(exists("Ecdfbvn"))
        Eres <- vput(Eres,get("Ecdfbvn",env=evei),"Ecdfbvn")
      if(exists("EdirTol"))
        Eres <- vput(Eres,get("EdirTol", env=evei),"EdirTol")
      if(exists("EcdfTol"))
        Eres <- vput(Eres,get("EcdfTol",env=evei),"EcdfTol")
      if(exists("EvTol"))
        Eres <- vput(Eres,get("EvTol", env=evei),"EvTol")
      if(exists("EdoML"))
        Eres <- vput(Eres,get("EdoML", env=evei),"EdoML")
      if(exists("EdoML.phi"))
        Eres <- vput(Eres,get("EdoML.phi", env=evei),"doml.phi")
      if(exists("EdoML.vcphi"))
        Eres <- vput(Eres,get("EdoML.vcphi", env=evei),"doml.vc")
      if(exists("EdoSim"))
        Eres <- vput(Eres,get("EdoSim", env=evei),"EdoSim")
      if(exists("Eeta"))
        Eres <- vput(Eres,get("Eeta", env=evei),"Eeta")
      if(exists("Eigraph.bvsmth")){
        Eres <- vput(Eres,get("Eigraph.bvsmth",env=evei),"eigraph.bvsmth")
        Eres <- vput(Eres,get("Eigraph.bvsmth", env=evei),  "bvsmth")
      }
      if(exists("EisChk"))
        Eres <- vput(Eres,get("EisChk", env=evei),"EisChk")
      if(exists("EiLliks"))
        Eres <- vput(Eres,get("EiLliks", env=evei),"EiLliks")
      if(exists("EisFac"))
        Eres <- vput(Eres,get("EisFac", env=evei),"EisFac")
      if(exists("Eisn"))
        Eres <- vput(Eres,get("Eisn", env=evei),"Eisn")
      if(exists("Eist"))
        Eres <- vput(Eres,get("Eist", env=evei),"Eist")
      Eres <- vput(Eres,get("EmaxIter", env=evei),"EmaxIter")
      Eres <- vput(Eres,get("EnonEval",  env=evei),"EnonEva") 
      Eres <- vput(Eres,get("EnonNumInt", env=evei),"EnonNum")
      Eres <- vput(Eres,get("EnonPar", env=evei),"EnonPar")
      Eres <- vput(Eres,get("EnumTol", env=evei),"EnumTol")
      Eres <- vput(Eres,get("Erho", env=evei),"Erho")
      Eres <- vput(Eres,get("Eselect",env=evei),"Eselect")
      Eres <- vput(Eres,get("EselRnd",env=evei),"EselRnd")
      Eres <- vput(Eres,get("Esigma",  env=evei),"Esigma")
      Eres <- vput(Eres,get("Esims",  env=evei),"Esims")
      Eres <- vput(Eres,get("Estval", env=evei),"Estval")
      Eres <- vput(Eres,get("ei.vc",  env=evei),"ei.vc")
      assign("Eres",Eres, env=evbase)
      return(Eres)
    }
    if (round == 2){
       Eres <- vput(Eres,get("betaBs", env=evei),"betabs");
       Eres <- vput(Eres,NA,"retcode")
       Eres <- vput(Eres,NA,"phi")
       Eres <- vput(Eres,NA,"loglik")
       Eres <- vput(Eres,NA,"ghactual");
       Eres <- vput(Eres,NA,"vcphi")
       if(!"Esims" %in% names(Eres)) 
         Eres <- vput(Eres,get("Esims", env=evei),"Esims")
       if(!"ei.vc" %in% names(Eres)) 
         Eres <- vput(Eres,get("ei.vc",  env=evei),"ei.vc")
       
       assign("Eres",Eres, env=evbase)
      return(Eres)
     }

     if(round == 3){
       Eres <- vput(Eres,NA,"retcode")
       Eres <- vput(Eres,get("EdoML.phi",  env=evei),"phi")
       Eres <- vput(Eres,NA,"loglik")
       Eres <- vput(Eres,NA,"ghactual")
       Eres <- vput(Eres,get("EdoML.vcphi",  env=evei),"vcphi")
       Eres <- vput(Eres,get("ei.vc", env=evei),"ei.vc")
       assign("Eres",Eres, env=evbase)
       return(Eres)
     }
  ####added by Ev
     if(round==4){
       
       Eres <- vput(Eres, get("phihat", env=evbase),"phihat");
       Eres <- vput(Eres,get("MLpsi",env=evbase),"phi");
       Eres <- vput(Eres,get("MLvc", env=evbase),"vcphi");
       return(Eres);
     }
       
    return(Eres)
        
  }
  
settTruthBnds <- function(betaB, betaW, bnd, truthB, truthW, bnds, evei)
{
  evbase <- get("evbase", env=parent.frame())
  truthB <- truthW <- bnds <- as.matrix(0)
  vec <- c(list(truthB), list(truthW), list(bnds))
  names(vec) <- c("truthB", "truthW", "bnds")
  evbase <- setGlobals(evbase, lst= vec)
  truthB <- betaB;
  assign("truthB", truthB, env=evei)
  truthW <- betaW;
  assign("truthW", truthW, env=evei)
  bnds <- bnd;
  assign("bnds", bnds, env=evei)
  return(evei)

}
