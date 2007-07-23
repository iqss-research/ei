###
##  tst = vin(dbuf,"var");
##
##  dbuf = R object
##  var  = name of variable
##
##  tst = 1 if var is in dbuf and
##        0 if not

vin<-function(dbuf,str){
        return(str %in% names(dbuf))
}

###
##  t = scalzero(y);
##
##  y = matrix
##  t = 1 if y is a scalar zero and 0 otw
##

scalzero<-function(y){
        if (is.na(y))
          return(0)
        if (nrow(y)==1 && ncol(y)==1 && y==0)
          return(1)
        else
          return(0)
}

###
##    t = scalone(y);
##
##    y = matrix
##    t = 1 if y is a scalar one and 0 otw
##

scalone<-function(y){
        if (is.na(y))
          return(0)
        if (nrow(y)==1 && ncol(y)==1 && y==1)
          return (1)
        else
          return (0)
}


###
##  y = meanWc(x,wt);
##
##  x = NxM matrix to be avg'd
##  wt = scalar 1 or Nx1 or NxM weight used in averaging
##
##  works with missing values; packs rowwise.

meanwc<-function(x,wt){
        if(all(is.na(wt)) || wt==1)
          wt<-rep(1,nrow(x))
        wwt<-wt
        wwt[is.na(wwt)]<-0
        xtmp<-x
        xtmp[is.na(xtmp)]<-0
        if(is.matrix(x))
          res<-(apply((xtmp*wwt),2,sum))/(apply((1-ismissm(x+wt))*wwt,2,sum))
        else
          res<-sum(xtmp*wwt)/(sum(1-ismissm(x+wt))*wwt)
        return(res)
}


###
##  y = ismissm(x)
##
##      x = an n x m matrix
##
##      y = an n x m matrix of 0's (indicating not missing) and 1's (missing)
##

ismissm<-function(x){
        x[!(is.na(x))]<-0
        x[is.na(x)]<-1
        return(x)
}


###
##  create vector of PTS evenly spaced points between STRT and ENDD,
##  including the end points.

seqase<-function(strt,endd,pts){

        t<-(endd-strt)/(pts-1)
        res<-seq(strt,endd,t)
        return (res)
}

