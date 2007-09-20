### DESCRIPTION corresponds to Gary's eidata.src
###             for pulling apart and putting together the data
###             Translation of Gary King Gauss code
###
### Elena Villalon (evillalon@iq.harvard.edu)

packdta <- function(x,Zb,Zw,t, evbase=NULL){
  if(!length(evbase))
    evbase <- get("evbase", env=parent.frame())
  Ez <- get("Ez", env=evbase)
  Eselect <- get("Eselect", env=evbase)
  dataset <- cbind(x,t)
 
  if (Ez[2]>1)
    dataset <- cbind(Zw,dataset)
  

  if (Ez[1]>1)
    dataset <- cbind(Zb,dataset)
 
  if( nrow(Eselect)!=1)
    dataset <- subset(dataset,subset=Eselect); ###selif is subset
  
  return(dataset)
}

pluckdta <- function(dta, evbase=NULL){
  if(!length(evbase))
    evbase <- get("evbase", env=parent.frame())
  Ez <- get("Ez", env=evbase)
  if(Ez[1]>1)
    Zb <- dta[,1:(Ez[1]-1)]
  else
    Zb <- 1;
  
  if(Ez[2]>1)
    Zw <- dta[,Ez[1]:(sum(Ez)-2)]
  else
    Zw <- 1;
  

  x <- dta[,ncol(dta)-1];
  y <- dta[,ncol(dta)];
  lst <- c(list(Zb=Zb), list(Zw=Zw), list(x=x), list(y=y))
  return(lst);
}


