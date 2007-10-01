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

expanddots <- function(drvdot, drv, evbase){
 
  args <- names(drv)
  args <- args[-1]
  argsdots <- names(drvdot)
  argsdots <- argsdots[-1]
 
  if(length(argsdots) == length(args))
    return(evbase)  ##nothing in expand.dots
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
      
    if(length(as.list(drvdot[[n]])) > 1){
      vec <- as.character(drvdot[[n]])
  
      lst <- NULL
 
      for(m in 2:length(vec))
        lst <- c(lst,as.list(drvdot[[n]])[[m]])
  
      val <- as.matrix(lst)
    }else
    val <- drvdot[[n]]
    
    assign(argsdots[n], val, env=evbase)
        
  }
  
  return(evbase)

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
  
  ass <- lapply(param, function(att,evfrom,evto) {
  
    val <- try(get(att, env=evfrom))
  
    if(class(val) == "try-error")
      message("Error getting param")
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
add.to.Eres <- function(Eres=list(), round=1, evbase=NULL)
{
  if(!length(evbase))
    evbase <- eiset()
  
  evei <- getEnvVar(evbase, environment())
        
 ### inputs
    if(round == 1){
      Eres <- vput(Eres,t,"t");
      Eres <- vput(Eres,x,"x");
      Eres <- vput(Eres,tvap,"n");
      Eres <- vput(Eres,Zb,"Zb");
      Eres <- vput(Eres,Zw,"Zw");
  
## essential globals 
      Eres <- vput(Eres,EalphaB,"EalphaB");
      Eres <- vput(Eres,EalphaW,"EalphaW");
      Eres <- vput(Eres,Ebeta,"Ebeta");
      Eres <- vput(Eres,Ebounds,"Ebounds");
      Eres <- vput(Eres,Ecdfbvn,"Ecdfbvn");
      Eres <- vput(Eres,EdirTol,"EdirTol");
      Eres <- vput(Eres,EcdfTol,"EcdfTol");
      Eres <- vput(Eres,EvTol,"EvTol");
      Eres <- vput(Eres,EdoML,"EdoML");
      Eres <- vput(Eres,EdoML.phi,"doml.phi");   
      Eres <- vput(Eres,EdoML.vcphi,"doml.vc");   
      Eres <- vput(Eres,EdoSim,"EdoSim");
      Eres <- vput(Eres,Eeta,"Eeta");
      Eres <- vput(Eres,Eigraph.bvsmth,   "eigraph.bvsmth"); 
      Eres <- vput(Eres,EisChk,"EisChk");
      Eres <- vput(Eres,EiLliks,"EiLliks");
      Eres <- vput(Eres,EisFac,"EisFac");
      Eres <- vput(Eres,Eisn,"Eisn");
      Eres <- vput(Eres,Eist,"Eist");
      Eres <- vput(Eres,EmaxIter,"EmaxIter");
      Eres <- vput(Eres,EnonEval,"EnonEva"); 
      Eres <- vput(Eres,EnonNumInt,"EnonNum");
      Eres <- vput(Eres,EnonPar,"EnonPar");
      Eres <- vput(Eres,EnumTol,"EnumTol");
      Eres <- vput(Eres,Erho,"Erho");
      Eres <- vput(Eres,Eselect,"Eselect");
      Eres <- vput(Eres,EselRnd,"EselRnd");
      Eres <- vput(Eres,Esigma,"Esigma");
      Eres <- vput(Eres,Esims,"Esims");
      Eres <- vput(Eres,Estval,"Estval");
      Eres <- vput(Eres,ei.vc,"ei.vc");
      assign("Eres",Eres, env=evbase)
      return(Eres)
    }
    if (round == 2){
    ###   Eres <- vput(Eres,betaBs,"betabs");
       Eres <- vput(Eres,NA,"retcode");
       Eres <- vput(Eres,NA,"phi");
       Eres <- vput(Eres,NA,"loglik");
       Eres <- vput(Eres,NA,"ghactual");
       Eres <- vput(Eres,NA,"vcphi");  
       Eres <- vput(Eres,Esims,"Esims");
       Eres <- vput(Eres,ei.vc,"ei.vc");
       assign("Eres",Eres, env=evbase)
      return(Eres)
     }

     if(round == 3){
       Eres <- vput(Eres,NA,"retcode");
       Eres <- vput(Eres,EdoML.phi,"phi");
       Eres <- vput(Eres,NA,"loglik");
       Eres <- vput(Eres,NA,"ghactual");
       Eres <- vput(Eres,EdoML.vcphi,"vcphi");
       Eres <- vput(Eres,ei.vc,"ei.vc");
       assign("Eres",Eres, env=evbase)
       return(Eres)
     }
     
    return(Eres)
        
  }
  
settTruthBnds <- function(betaB, betaW, bnd, truthB, truthW, bnds, evei)
{
  evbase <- get("evbase", env=parent.frame())
  truthB <- truthW <- bnds <- as.matrix(0)
  evbase <- setGlobals(evbase, vec = c(truthB=truthB, truthW=truthW, bnds=bnds))
  truthB <- betaB;
  assign("truthB", truthB, env=evei)
  truthW <- betaW;
  assign("truthW", truthW, env=evei)
  bnds <- bnd;
  assign("bnds", bnds, env=evei)
  return(evei)

}
