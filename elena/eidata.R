##/* 
##**  This archive is part of the program EI
##**  (C) Copyright 1995-2001 Gary King
##**  All Rights Reserved.
##**
##**  For pulling apart and putting together the data
##**
##*/
##/*
##  {Zb,Zw,x,t} = pluckdta(dataset);
##**
##*/
### DESCRIPTION corresponds to Gary King's eidata.src
###             for pulling apart and putting together the data
###             Translation of Gary King Gauss code
###
### Elena Villalon (evillalon@iq.harvard.edu)

packdta <- function(x,Zb,Zw,t, evbase=NULL, Ez=matrix(1,nrow=2, ncol=1), Eselect=as.matrix(1)){
  if(!length(evbase))
    evbase <- try(get("evbase", env=parent.frame()),silent=TRUE)
  if(class(evbase)!="try-error"){
    Ez <- get("Ez", env=evbase)
    Eselect <- get("Eselect", env=evbase)
  }
  x <- as.matrix(x)
  t <- as.matrix(t)

  dataset <- cbind(x,t)
  
  if (length(Ez) > 1 && Ez[2]>1)
    
    dataset <- cbind(Zw,dataset)
  

  if (Ez[1]>1)
    dataset <- cbind(Zb,dataset)
 
  if( rows(Eselect)!=1)
    dataset <- subset(dataset,subset=Eselect); ###selif is subset
 
  return(dataset)
}

pluckdta <- function(dta, evbase=NULL,Ez=matrix(1,nrow=2, ncol=1)){
  if(!length(evbase)){
    evb <- try(get("evbase", env=parent.frame()))
    if(class(evb) != "try-error") evbase <- evb
    Ez <- get("Ez", env=evbase)
  }
  if(Ez[1]>1)
    Zb <- dta[,1:(Ez[1]-1)]
  else
    Zb <- 1;
  
  if(Ez[2]>1)
    Zw <- dta[,Ez[1]:(sum(Ez)-2)]
  else
    Zw <- 1;
  

  x <- dta[,cols(dta)-1];
  y <- dta[,cols(dta)];
  lst <- c(list(Zb=Zb), list(Zw=Zw), list(x=x), list(y=y))
  return(lst);
}


