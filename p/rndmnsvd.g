/*
** Multivariate Normal Sampling when the negative of the Hessian is
** not positive definite. This uses the singular value decomposition to draw
** directly from the singular normal density.
**
   res = rndmnsvd(mean,invvc,sims,bounds,tol);
**
** INPUTS:
**   mean = k x 1 vector of means
**   invvc = k x k matrix, the negative of the hessian
**   sims = scalar, the number of simulations
**   bounds = k x 2 matrix, upper bounds ~ lower bounds
**   tol = scalar, the tolerance level for the diagonal element of s
**         where {u,s,v}=svd1(invvc) (for example 10^(-3))
**
** OUTPUTS:
**   res = sims x k matrix of the multivariate singular normal random draws
**
*/

proc (1)=rndmnsvd(mean,invvc,sims,bounds,tol);
   local u,s,v,res,k,indx,dist,temp,bnds,limsim;
 
   k=rows(mean);
   /* some input checks */
   if sumc(bounds[.,2].>=bounds[.,1])/=0;
     "error(rndmnsvd): Upper bounds must be greater than lower bounds.";
     res={.};
     retp(res);
   endif; 
   if rows(invvc)/=rows(invvc);
     "error(rndmnsvd): the Inverse of variance matrix has to be a square matrix.";
     res={.};
     retp(res);
   endif;
   if k/=rows(invvc);
     "error(rndmnsvd): the dimensions of the Inverse of variance matrix and";
     "mean vector have to be the same.";
     res={.};
     retp(res);
   endif;

   {s,u}=eighv(invvc); @ eigen value decomposition @
   indx=makefac(k,2);
   dist=zeros(2^k,1);
   for i(1,rows(indx),1);
     temp=-mean;
     for j(1,k,1);
       temp[j]=temp[j]+bounds[j,indx[i,j]];
     endfor;
     dist[i]=temp'*temp;
   endfor;
   dist=maxc(dist);
   bnds=dist~(-dist);      @ transformed bounds @
   mean=u'*mean;           @ transformed mean @
   res=zeros(k,sims);
   for i(1,sims,1);
     limsim=1;
     do while (sumc((bounds[.,1] .< res[.,i])+(bounds[.,2] .> res[.,i]))>0)
              or limsim==1;
       for j(1,k,1);
         if (s[j] > tol);
           res[j,i]=rndtni(mean[j],1/s[j],bnds[2]~bnds[1]);
         else;
           res[j,i]=rndu(1,1)*(2*bnds[1])+bnds[2];
         endif;
       endfor;
       res[.,i]=u*res[.,i];
       limsim=limsim+1;
       if limsim==100000;
         "error(rndmnsvd): the sampling method failed. adjust the bounds.";
         res={.};
         retp(res);
       endif;
     endo;
   endfor;
       
   retp(res');
endp;







