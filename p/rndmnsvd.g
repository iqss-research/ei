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
**   bounds = 2 x k matrix, lower bounds | upper bounds
**         or 2 x 1 vector to use the same bounds for all parameters
**   tol = scalar, the tolerance level for the diagonal element of s
**         where {u,s,v}=svd1(invvc) (for example 10^(-3))
**
** OUTPUTS:
**   res = sims x k matrix of the multivariate singular normal random draws
**
*/

proc (1)=rndmnsvd(mean,invvc,sims,bounds,tol);
   local u,s,v,res;

   /* some input checks */
   if sumc(bounds[1,.].>=bounds[2,.])/=0;
     "error(rndmnsvd): Upper bounds must be greater than lower bounds.";
     res={.};
     retp(res);
   endif; 
   if rows(invvc)/=rows(invvc);
     "error(rndmnsvd): the Inverse of variance matrix has to be a square matrix.";
     res={.};
     retp(res);
   endif;
   if rows(mean)/=rows(invvc);
     "error(rndmnsvd): the dimensions of the Inverse of variance matrix and";
     "mean vector have to be the same.";
     res={.};
     retp(res);
   endif;

   {u,s,v}=svd1(invvc);
/* OLD: KOSUKE, PLS CHECK MY CHANGES AND DELETE THE OLD VERSION IF THEY'RE OK.
   v=s+(diag(s) .< tol).*eye(rows(s));
   res=u'*(rndmn(u*mean,invpd(v),sims))';
   res=res-(diag(s) .< tol).*res;       
   res=res+(diag(s) .< tol).*(rndu(rows(mean),sims).*(bounds[2,.]-bounds[1,.])'
       +bounds[1,.]');
*/
   t=(diag(s) .< tol);
   v=s+t.*eye(rows(s));
   res=(1-t).*(u'*(rndmn(u*mean,invpd(v),sims))')
       +t.*((rndu(rows(mean),sims).*(bounds[2,.]-bounds[1,.])'+bounds[1,.]'));
       
   retp(res');
endp;



