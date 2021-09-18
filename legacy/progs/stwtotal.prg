/*
** This program gives aggregate (state-level) results 
** for the state-wide primary races. 
**
** It uses stwprint.g procedure to print out the results
** in a nice format.
**
*/

new;
library ei;

/** inputs **/
optfile = "stw_pr00p.out"; @ output file name @
n_dist=30;                 @ the number of districts @

/** initializing variables **/
alwbd=zeros(1,3);
aupbd=zeros(1,3);
betaH=zeros(1,100);
betaB=zeros(1,100);
betaW=zeros(1,100);
stwN=zeros(1,3);
stwNt=0;

/*** main program ***/
_Eprt=0;

output file = ^optfile on;
for i(1,n_dist,1);
  name="db"$+ftos(i+0,"*.*lf",1,0);
  loadm ndbuf=^name;


  format/rd 4,0;
  "**************************************************************";
  i "th District";
  "**************************************************************";
  stwprint(ndbuf);
  format/rd 7,3;
  
  /*** read in from the data buffer ***/
  dbuf1=vread(ndbuf,"dbuf1");
  dbuf3=vread(ndbuf,"dbuf3");

  n=sumc(eiread(dbuf1,"N"));
  nH=sumc(eiread(dbuf1,"Nb"));
  nW=sumc(eiread(dbuf3,"Nb"));
  nB=n-nH-nW;

  xH=nH/n;
  xW=nW/n;
  xB=nB/n;

  t=sumc(eiread(dbuf1,"Nt"))/n;

  aboundbH=eiread(dbuf1,"abounds");
  aboundbH=aboundbH[.,1];  
  aboundbW=eiread(dbuf3,"abounds");
  aboundbW=aboundbW[.,1];  
  {temp,temp1}=bounds1(eiread(dbuf1,"t"),1-eiread(dbuf1,"x")-eiread(dbuf3,"x"),eiread(dbuf1,"N"));
  aboundbB=temp1[1,.]';

  betabsH=eiread(dbuf1,"aggs");
  betabsH=betabsH[.,1];  
  betabsW=eiread(dbuf3,"aggs");
  betabsW=betabsW[.,1];
  betabsB=(t-betabsH.*xH-betabsW.*xW)./xB;

  alwbd=alwbd|(aboundbH[1]~aboundbB[1]~aboundbW[1]);
  aupbd=aupbd|(aboundbH[2]~aboundbB[2]~aboundbW[2]);
  betaH=betaH|betabsH';
  betaB=betaB|betabsB';
  betaW=betaW|betabsW';
  stwN=stwN|(nH~nB~nW);
  stwNt=stwNt|sumc(eiread(dbuf1,"Nt"));

endfor;

alwbd=alwbd[2:rows(alwbd),.];
aupbd=aupbd[2:rows(aupbd),.];
betaH=betaH[2:rows(betaH),.];
betaB=betaB[2:rows(betaB),.];
betaW=betaW[2:rows(betaW),.];
stwN=stwN[2:rows(stwN),.];
stwNt=stwNt[2:rows(stwNt),.];
stwX=sumc(stwN)/sumc(sumc(stwN));

stwlwbd=meanwc(alwbd,stwN);
stwupbd=meanwc(aupbd,stwN);
stwlwbde=stwlwbd.*sumc(stwN);
stwupbde=stwupbd.*sumc(stwN);

stwbetab=zeros(100,3);
stwNbt=zeros(100,3);
for i(1,100,1);
  stwbetab[i,.]=meanwc(betaH[.,i],stwN[.,1])~meanwc(betaB[.,i],stwN[.,2])~meanwc(betaW[.,i],stwN[.,3]);
  stwNbt[i,.]=stwbetab[i,.].*stwX'.*sumc(sumc(stwN));
endfor;

format/rd 4,0;
?;
"*******************************************************";
"State-wide results";
"*******************************************************";
?;
"DEMOGRAPHIC STATISTICS AND ELECTION RESULTS";
"Note: These statistics summarize the data"; 
"      available for EI estimation.";
?;
print "N " sumc(sumc(stwN));
?;
format/rd 7,1;
print "Percent Hispa:       " stwX[1]*100;
print "Percent Black:       " stwX[2]*100;
print "Percent White+Other: " stwX[3]*100;
?;
format/rd 10,0;
print "Hispans:       " sumc(stwN[.,1]);
print "Blacks:        " sumc(stwN[.,2]);
print "Whites+Other:  " sumc(stwN[.,3]);
?;
print "Number of Votes Cast: " sumc(stwNt);
format/rd 7,1;
print "Turnout %: " sumc(stwNt)/sumc(sumc(stwN))*100;
?;
"******************************************";
?;
"EI ESTIMATION RESULTS FOR HISPANICS, BLACKS AND WHITES";    
?; 
"** Betas **";	
?;
"Estimates of Aggregate Quantities of Interest";
format/rd 7,3;
"  BETA-H  BETA-B  BETA-W";;
meanc(stwbetab)'|stdc(stwbetab)';
?;
"Composition of electorate";
"  HISPAN  BLACK   WHITES";;
meanc(stwNbt./sumc(stwNbt'))'|stdc(stwNbt./sumc(stwNbt'))';
?;
"Aggregate Bounds";
"        BETA-H  BETA-B  BETA-W";
"Lower:" stwlwbd';
"Upper:" stwupbd';
?;
"Bounds for composition of electorate";
"        Hispan  BLACK   WHITES";
"Lower:" stwlwbde'/sumc(stwNt);
"Upper:" stwupbde'/sumc(stwNt);
?;
"*********************************************************************";
output off;

