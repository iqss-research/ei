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
   clearg u,s,v,res,k,indx,dist,temp,bnds,limsim,tmpmean,tmps,snum,midbounds,
   bnds1;
 
   k=rows(mean);
   /* some input checks */
   if sumc(bounds[.,2].<bounds[.,1])/=0;
     "error(rndmnsvd): Upper bounds must be greater than lower bounds.";
     res={.};
     retp(res);
   endif; 
   if rows(invvc)/=rows(invvc);
     "error(rndmnsvd): the -Hessian must be square.";
     res={.};
     retp(res);
   endif;
   if k/=rows(invvc);
     "error(rndmnsvd): the dimensions of the -Hessian matrix and";
     "mean vector must be the same.";
     res={.};
     retp(res);
   endif;

   {s,u}=eighv(invvc);  @ spectral decomposition: u*diag(s)*u'=invvc @
   s=recode(s, s.<tol, tol);
   snum=sims*10;           @ # of draws at a time @
 
   /* shift to the middle of the bounds */
   midbounds=sumc(bounds')/2;
   tmpmean=u'*(mean-midbounds); @ transformed mean @
   tmpmean=vec(tmpmean.*ones(k,snum));
   tmps=vec(s.*ones(k,snum));
   temp=bounds-midbounds;
   temp=maxc(abs(temp)');
   dist=sqrt(temp'*temp);
   bnds=(-dist)~dist;      @ transformed bounds @

   limsim=1;              
   res=zeros(1,k); 
   v=1./tmps;
   bnds1=bnds.*ones(rows(tmps),2);
   do while rows(res)<=sims;
     temp=rndtni(tmpmean,v,bnds1);
     temp=reshape(temp,rows(temp)/k,k)';
     temp=u*temp+midbounds;
     temp=selif(temp',sumc((bounds[.,1].<temp)+(bounds[.,2].>temp)).==0);
     if scalmiss(temp)/=1;
       res=res|temp;
     endif;
     "(rndmnsvd): trying" limsim "th time...";
     limsim=limsim+1;
     if limsim==50;
       "error(rndmnsvd): the sampling method failed. adjust the bounds.";
       retp(miss(1,1));
     endif;
   endo;
   retp(res[2:sims+1,.]);
endp;


