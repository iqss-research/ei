new;
library ei,cml,pgraph;
rndseed 92232;
fmtt;

/* set parameters */
obs=100;
o=ones(obs,1);
bb0= 0.25*o;  bw0=0.5*o; sb0=.1; sw0=.2; rho0=-.2;
bb=.75*o;     bw=.25*o;   sb=.1;  sw=.2;  rho=.2;
x=seqas(0,1,obs);
n=int(rndu(obs,1)*1000)+100;


/* generate data randomly */
beta=rndbtni(bb0,bw0,sb0,sw0,rho0,0~1~0~1); 
lambda=rndbtni(bb,bw,sb,sw,rho,0~1~0~1); 
t=beta[.,1].*x+beta[.,2].*(1-x);
xx=beta[.,1].*x./t;
v=lambda[.,1].*xx+lambda[.,2].*(1-xx);
"truth lambdas=";;meanc(lambda)';
"truth betas  =";;meanc(beta)';
save beta,lambda;

output file=sample2.asc reset;
v~t~x~n;
output off;

