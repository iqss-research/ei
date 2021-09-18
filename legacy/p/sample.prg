/* set parameters */
new;
library ei;
#include rndp.src;
rndseed 9828;
eiset;
bb=.2; bw=.8; sb=.1; sw=.2; rho=.7;
sims=75;
fmtt;

/* simulate data */
x=rndu(sims,1);
x1=1-x;
bb=bb*ones(sims,1);
bw=bw*ones(sims,1);
d=rndbtni(bb,bw,sb,sw,rho,0~1~0~1);
t=d[.,1].*x+d[.,2].*x1;
n=rndp(sims,1,500);

output file=sample.asc reset;
  screen off;
  t~x~n;
  screen on;
output off;
" sample.asc saved";
