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
   local u,s,v,res,k,indx,dist,temp,bnds,limsim,tmpmean,tmps,snum;
 
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

   {s,u}=eighv(invvc);  @ eigen value decomposition: u*diag(s)*u'=invvc @
   s=recode(s, s.<tol, tol);
   snum=sims*10;           @ # of draws at a time @
 
/* if you want to shift to the mean */
   @ indx=makefac(k,2);
   dist=zeros(2^k,1);
   for i(1,rows(indx),1);
     temp=-mean;       
     for j(1,k,1);
       temp[j]=temp[j]+bounds[j,indx[i,j]];
     endfor;
     dist[i]=sqrt(temp'*temp);
   endfor;
   tmpmean=zeros(snum*k,1); 
   tmps=vec(s.*ones(k,snum));
   dist=maxc(dist); 
   bnds=(-dist)~dist; @

   /* shift to the middle of the bounds */
   tmpmean=u'*(mean-(bounds[.,1]+bounds[.,2])/2); @ transformed mean @
   tmpmean=vec(tmpmean.*ones(k,snum));
   tmps=vec(s.*ones(k,snum));
   temp=bounds-(bounds[.,1]+bounds[.,2])/2;
   temp=maxc(abs(temp)');
   dist=sqrt(temp'*temp);
   bnds=(-dist)~dist;      @ transformed bounds @

   limsim=1;              
   res=zeros(1,k); 
   do while rows(res)<=sims;
     temp=rndtni(tmpmean,1/tmps,bnds.*ones(rows(tmps),2));
     temp=reshape(temp,rows(temp)/k,k)';
     @ temp=u*temp+mean; @
     temp=u*temp+(bounds[.,1]+bounds[.,2])/2;
     temp=selif(temp',sumc((bounds[.,1].<temp)+(bounds[.,2].>temp)).==0);
     if scalmiss(temp)/=1;
       res=res|temp;
     endif;
     "(rndmnsvd): trying" limsim "th time...";
     limsim=limsim+1;
     if limsim==50;
       "error(rndmnsvd): the sampling method failed. adjust the bounds.";
       res={.};
       retp(res);
     endif;
   endo;
   retp(res[2:sims+1,.]);
endp;


