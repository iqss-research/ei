###
##  l = homoindx(x);
##  l is a list with 3 elements
##  c  = vector of index numbers for heterogeneous precincts
##  c0 = vector of index numbers for homogenous white (x=0) precincts
##  c1 = vector of index numbers for homogenous black (x=1) precincts
###

homoindx<-function(x){
        .EnumTol<-0.001
        res<-list()
        indx<-c(1:length(x))                    # depend how x is going to be passed, as matrix or array c()???
        c0<-(x<(.EnumTol))
        c1<-(x>(1-(.EnumTol)))
        c<-1-c0-c1
        res$c0<-indx[c0==TRUE]
        res$c1<-indx[c1==TRUE]
        res$c<-indx[c==TRUE]
        return (res)
}
