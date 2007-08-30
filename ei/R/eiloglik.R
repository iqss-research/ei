###
##  l = homoindx(x);
##  INPUT: is a matrix of one column or array.
##
##  OUTPUT: a list with 3 elements
##          c  = vector of index numbers for heterogeneous precincts
##          c0 = vector of index numbers for homogenous white (x=0) precincts
##          c1 = vector of index numbers for homogenous black (x=1) precincts
##          Translates code in Gauss written by Gary King
## Ferdinand Alhimadi & Elena Villalon (evillalon@iq.harvard.edu)
## Date August 17th, 2007

homoindx<-function(x){
        EnumTol <- try(grep("EnumTol", env=grep("evbase", env=parent.frame())))
        if(class(EnumTol) == "try-error")
          EnumTol <- 0.0001
        x <- as.matrix(x)
        indx <- seq(1,nrow(x), 1)
        
        res<-list()
       
        c0<- x< EnumTol
        c1<- x>(1- EnumTol)
        c<-(1-c0-c1)
        cl <- as.logical(c)
        res$c0<-indx[x <EnumTol]
        res$c1<-indx[x >(1-EnumTol) ]
        res$c<-indx[c==TRUE]
        return (res)
}



                
