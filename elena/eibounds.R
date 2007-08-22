##    {Bs,aggs}=bounds1(t,x,n);
##
## computes bounds on parameters given aggregate data
##
## INPUTS: nx1 vectors, unit of analysis is precinct
## see output of sumvar()
##
## OUTPUTS: bounds on precinct-level parameters
## Bs   = cols: lower-black ~ upper-black ~ lower-white ~ upper-white
## aggs =  bounds on district aggregates
##         cols: lower ~ upper
##         rows: beta-b, beta-w
###
### USES: recode
###
### Translation of the Gauss code by Gary King
### AUTHOR: Ferdinand Alhimadi & Elena Villalon
###         evillalon@iq.harvard.edu
###
#include ei.ext;
bounds1<-function(t,x,n){
        ## local LbetaB,UbetaB,LbetaW,UbetaW,aggs,omx,Nb,Nw,c,c0,c1,p,tx,tomx,z,o,m;
        omx<-1-x;
        Nb<-x*n;                       #nr of black people
        Nw<-omx*n;                     #nr of white people
        ##  {c,c0,c1} = homoindx(x);       # proc in eiloglik.src
        cs<-homoindx(x)
        c<-cs$c
        c0<-cs$c0
        c1<-cs$c1
        p<-length(x);                    # nr of precinct could by nrow(data) if the args are dataset
        
        LbetaB<-matrix(0, nrow=p,ncol=1)
        UbetaB<-matrix(0,nrow=p,ncol=1)
        LbetaW<- matrix(0, nrow=p,ncol=1)
        UbetaW<-matrix(0, nrow=p,ncol=1)
        z<-matrix(0, nrow=p,ncol=1)
        o<-matrix(1, nrow=p,ncol=1)
        m<-o*NA;
        
        cna <- na.omit(c)
        if (length(c) >1 || length(cna)){
                tx<-t[c]/x[c]
                tomx<-t[c]/omx[c]
                LbetaB[c]<-maxr(z[c],tx-(omx[c]/x[c]))
                UbetaB[c]<-minr(tx,o[c])
                LbetaW[c]<-maxr(z[c],tomx-(x[c]/(1-x[c])))
                UbetaW[c]<-minr(tomx,o[c])
        }
        c0na <- na.omit(c0)
        if (length(c0) > 1 || length(c0na)){## homogeneously white 
                LbetaB[c0]<-m[c0]
                UbetaB[c0]<-m[c0]
                LbetaW[c0]<-t[c0]
                UbetaW[c0]<-t[c0]
        }
        c1na <- na.omit(c1)
        if (lengt(c1) > 1 || length(c1na)){## homogeneously black @
                LbetaB[c1]<-t[c1]
                UbetaB[c1]<-t[c1]
                LbetaW[c1]<-m[c1]
                UbetaW[c1]<-m[c1]
        }
        
### fix rounding errors due to machine precision */
### basically change any negative value to 0 and any value >1 to 1
        LbetaB=recode(LbetaB,cbind((LbetaB<0),(LbetaB>1)),c(0,1))
        UbetaB=recode(UbetaB,cbind((UbetaB<0),(UbetaB>1)),c(0,1))
        LbetaW=recode(LbetaW,cbind((LbetaW<0),(LbetaW>1)),c(0,1))
        UbetaW=recode(UbetaW,cbind((UbetaW<0),(UbetaW>1)),c(0,1))
        
        res<-list()
        res$aggs<-rbind(cbind(weighted.mean(LbetaB,Nb),weighted.mean(UbetaB,Nb)),
                        cbind(weighted.mean(LbetaW,Nw),weighted.mean(UbetaW,Nw)))
        res$bs<-cbind(LbetaB,UbetaB,LbetaW,UbetaW)
        
        return(res)
}


### support proc
## row maximum
###
maxr<-function(a,b){
        res<-apply(cbind(a,b),1,max)
        return (res)
}

### support proc
## row minimum
###
minr<-function(a,b){
        res<-apply(cbind(a,b),1,min)
        return (res)
}

