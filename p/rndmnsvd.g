/*
** Multivariate Normal Sampling when the inverse of variance matrix is
** not invertible. This uses the singular value decomposition to avoid
** using the variance matrix.
**
   res = rndmnsvd(mean,invvc,sims,bounds,tol);
**
** INPUTS:
**   mean = k x 1 vector of means
**   invvc = k x k matrx of Inverse of variance covariance matrix
**   sims = scalar, the number of simulations
**   bounds = 2 x k matrix, lower bounds | upper bounds
**          or 2 x 1 vector if the same bounds can be used for all the
**          parameters
**   tol = scalar, the tolerance level for the infinite variance 
**
** OUTPUTS:
**   res = sims x k matrix of the multivariate normal random draws
**
*/

proc (1)=rndmnsvd(mean,invvc,sims,bounds,tol);
   local u,s,v,tmpmean,res;

   /* some input checks */
   if sumc(bounds[1,.].>=bounds[2,.])/=0;
     "error(rndmnsvd): Upper bounds must be greater than lower bounds.";
     res={.};
     retp(res);
   endif;
   if rows(invvc)/=rows(invvc);
     "error(rndmnsvd): the Inverse of variance covariance matrix has to be a square matrix.";
     res={.};
     retp(res);
   endif;
   if rows(mean)/=rows(invvc);
     "error(rndmnsvd): the dimensions of the Inverse of variance covariance matrix and the";
     "mean vector have to be the same.";
     res={.};
     retp(res);
   endif;

   {u,s,v}=svd1(invvc);
   tmpmean=u*mean;
   res=zeros(rows(mean),sims);
   for i (1,rows(mean),1);
     if (1/s[i,i])<tol;
       res[i,.]=rndn(1,sims)/sqrt(s[i,i])+tmpmean[i];
     else;
       res[i,.]=rndu(1,sims).*(bounds[2,.]-bounds[1,.])+bounds[1,.];
     endif;
   endfor;
   for i (1,sims,1);
     res[.,i]=inv(u)*res[.,i];
   endfor;
   retp(res');
endp;