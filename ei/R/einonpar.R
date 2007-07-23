### dist = perpdist(t,x,px,py);
##
## INPUTS:  t,x   = (px1) turnout, race
##          px,py = (Mx1) coordinates of points to compute distance
##
## OUTPUT: dist = (PxM) perpendicular distance from each point in px,py
##                 to the tomography line betaW=(t./(1-x))-(x./(1-x))*betaB
##/

perpdist<-function(t,x,px,py){
        ##  local A,B,C,D,uB,uW,lB,lW,acarg,tst1,tst2,tt,bnds;
        smtmp<-bounds1(t,x,1)
        bnds<-smtmp$bnds
        tt<-bounds1$tt
        lB<-bnds[,1]
        uB<-bnds[,2]
        lW<-bnds[,3]
        uW<-bnds[,4]
        A<-sqrt((uW-lW)^2+(uB-lB)^2)
        A[A<=0.0001]<-0.0001
        B<-sqrt((uW-t(py))^2+(lB-t(px))^2)
        B[b<=0.00001]<-0.00001
        C<-sqrt((lW-t(py))^2+(uB-t(px))^2)
        acarg<-(a^2+b^2-c^2)/(2*a*b)
        tst1<-(acarg>1)
        tst2<-(acarg<(-1))
        x<-acarg
        e<-(tst1 | tst2)
        v<-tst1-tst2
        x[(e==TRUE)]<-v[e==TRUE]
        ## acarg<-substute(acarg,(tst1.or tst2),0*acarg+tst1-tst2)
        D<-B*sin(acos(x))
        return (D)
}

###
##  area = unitarea(t,x);
##
##  t,x = (px1) turnout, race
##
##  area = (px1) area within unit square of truncated nonparametric
##         normal-sheet kernel for each tomography line
##
unitarea<-function(t,x){
        area<-c()
        x<-recode(x,cbind((x<Enumtol),(x>(1-.Enumtol))),c(.Enumtol,(1-.Enumtol)))
        smtmp<-bounds1(t,x,1);
        bnds<-smtmp$bnds
        tt<-smtmp$tt
        lB<-bnds[,1]
        uB<-bnds[,2]
        lW<-bnds[,3]
        uW<-bnds[,4]
        x1<-1-x
        xx1<-x/x1
        x1x<-x1/x
        nobs<-length(x)                                         # its just nr of obseration, again is depend how the input will be
        area<-zeros(nobs,1)
        var<-(.eigraph.bvsmth)^2
        
        for(i in 1:nobs){
                c<-cbind(maxr(1-uB[i],lW[i]),maxr(lB[i],1-uW[i]))
                ##   c=substute(c,c.<0.00001,c*0+0.00001);               later you can implement the function if necessary .
                c[c<0.00001]<-0.00001
                d<-perpdist(t[i+0],x[i+0],,c(1,0),c(0,1))                 # ??? dont forget to implement
                a<-sqrt(c^2-d^2)
                k<-a*cos(asin(a/c))
                g<-sqrt(a^2-k^2)
                
                ## coordinates of tomography line extended so that the end points
                ## are perpendicular to the 1,0 0,1 coordinates of the unit square:
                ## (uB+k[.,1])~(lW-g[.,1])~(lB-k[.,2])~(uW+g[.,2]);  */
                ## points to evaluate on the tomography line */
                
                betabS<-seqas(lB[i]-k[,2],uB[i]+k[,1],.EnonNumInt);
                betawS<-(t[i]/x1[i])-(x[i]/x1[i])*betabS;
                z=zeros(.EnonNumInt,1);
                o=ones(.EnonNumInt,1);
                
                ## lengths of perpendicular lines within the unit square */
                minbetaw<-maxr(z,betawS-x1x[i]*betabS)
                maxbetaw<-minr(o,betawS+x1x[i]*(1-betabS))
                minbetab<-maxr(z,betabS-xx1[i+0]*betawS)
                maxbetab<-minr(o,betabS+xx1[i+0]*(1-betawS))
                a<-sqrt((maxbetaw-betawS)^2+(maxbetab-betabS)^2)
                b<-sqrt((betabS-minbetab)^2+(betawS-minbetaw)^2)
                S<-((betabS>=minbetab)& (betabS<=maxbetab))
                
                ## area within square for each perpendicular line */
                area[i]<-meanc(cdfnorm(maxr(a,b),0,var)+((S-0.5)*2)*cdfnorm(minr(a,b)-S,0,var))
        }
        return(area)
}



###
## pz = nonbiv(t,x,px,py);
##
## INPUTS: t, x, = (px1) turnout, race
##         px,py = (MxQ) coordinates to evaluate nonparametric biv density at
##
## OUTPUT: pz = (MxQ) height of the nonparametric bivariate density at px,py,
##               which is f(px,py), using a sheet-normal kernel.
##
## If memory is available, use vec(px) and vec(py) and reshape pz upon output
## to make this proc run faster.
##

nonbiv<-function(t,x,px,py){
        ## local b,c,d,pz,z,r,col,i;
        col<-ncol(px)
        pz<-zeros(nrow(px),col)
        c <-unitarea(t,x)         # scale factor to divide by
        r<-length(x)
        for (i in 1:col){
                d <- perpdist(t,x,px[,i],py[,i])  # perpendicular distance to line #
                
                ## sheet-normal kernel */
                z <-d / (.eigraph.bvsmth)
                pz[,i] <- ((1 / sqrt(2*pi)) * sumc(exp(-0.5*(z*z))/c)) / r
        }
        
        pz<-pz/(.eigraph.bvsmth)^2
        return (pz)
}


###
## betabs = einonp(t,x);
##
## INPUTS: t,x = (p X 1) turnout, race
##
## OUTPUTS: betabs = (p x _Esims) simulations of beta^b,
##                   computed nonparametrically
##
einonp <- function(t,x){
        ##  local nobs,bnds,tt,i,j,betab,betabs,x1,pz,areas,cnts,fa,fb,a,b,d,e,c,ce,cd, betaw,res;
        nobs<-length(t)
        if (.Eprt>=2)
          print( "Nonparametric Density Estimation...")
        
        x<-recode(x,cbind((x<(.EnumTol)),(x>(1-(.EnumTol)))),c(.EnumTol,(1-(.EnumTol))))
        t<-recode(t,cbind((t<(.EnumTol)),(t>(1-(.EnumTol)))),c(.EnumTol,(1-(.EnumTol))))
        tmp<-bounds1(t,x,1)
        bnds<-tmp$aggs
        t<-tmp$bs
        
        betab<-zeros(.EnonEval,nobs)
        
        for (i in 1:nobs)
          betab[,i]<-seqase(bnds[i,1],bnds[i,2],.EnonEval)
        x1<-(1-x)
        print(x/x1)
        print(dim(betab))
        
        betaw<-t(t/x1)-t(x/x1)*betab
        
        if (.Eprt>=2)
          print ("Computing conditional density for each line...")
        pz<-nonbiv(t,x,betab,betaw)
        
        if (.Eprt>=2)
          print ("Drawing simulations...")
        
        ## compute interpolated trapazoidal areas */
        areas<-(betab[2,]-betab[1,])* na.omit(pz+0.5*abs(pz-lag(pz)));
        areas<-cumsumc(areas/t(sumc(areas)))
        
        betabs<-zeros(nobs,.Esims)
        for (i in 1:nobs){
                ## count within categories defined by areas of interpolated trapazoids */
                cnts<-counts(runif(.Esims),areas[,i])
                
                ## use inverse CDF to draw from trapazoidal dist  within categories */
                fa<-trimr(lag(pz[,i]),1,0)
                fb<-trimr(pz[,i],1,0)
                a<-trimr(lag(betab[,i]),1,0)
                b<-trimr(betab[,i],1,0)
                D<-fa*b-fb*a
                E<-fb-fa
                E[(abs(E)<0.000001)]<-0.000001      # @ fix for uniform distribution @
                C<-1/( (b-a)*D+0.5*(b^2-a^2)*E )
                CE<-C*E
                CD<-C*D
                
                res<-c()
                for (j in 1:(.EnonEval-1))
                  if (cnts[j+0]>0)
                    res<-rbind(res,((-CD[j]+sqrt((C[j]^2)*(D[j]^2)+ 2*CE[j]*(CD[j]*a[j]+0.5*CE[j]*(a[j]^2)+rndu(cnts[j],1))))/CE[j]))
                
                res<-t(trimr(res,1,0))
                betabs[i,]<-res[order(runif(.Esims))]
                
        }
        return (betabs)
}
