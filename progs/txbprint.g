/*
** 
   txbprint(ndbuf); 
**
** This proc prints out the EI turnout results for Hispanics, Blacks and 
** Whites+Others. The results are based on the two EI runs which are 
** Hispanics vs. Non-Hispanics and Whites+Others vs. Hispanics+Blacks.
**
** This proc does not give the estimates for Blacks if the population of 
** Blacks is less than 0.5 percent. It also does not give estimates for 
** Blacks if they are out of bounds etc.
**
** INPUT
**   ndbuf = a nested data buffer which includes two data buffers, dbuf1 
**           (Hispanic vs. Non-Hispanics) and dbuf2 (Whites+Others vs.
**           Hispanics+Blacks). It also has to have the two variables: 
**           "county" (county number) and "cname" ("county name") in dbuf1.
**
** OUTPUT
**   No output except it will print out the results on screen.
**
*/ 

proc (0)=txbprint(ndbuf);
 local p,c,dbuf1,dbuf2,cty,unicty,n,nb,t,x,nt,v,d,dv,betab,betaw,cbetab,
       nx,vx,bnds,bbnds1,bbnds2,paggs1,paggs2,abds1,betabs,betaws,abds2,
       temp,temp1,mask,fmt,prt,eprt,ctynm,goodman,double,a,b;
  
  eprt=_Eprt; _Eprt=0;
  /* initialize the parameters */
  p=rows(eiread(vread(ndbuf,"dbuf1"),"betab"));
  c=rows(unique(eiread(vread(ndbuf,"dbuf1"),"county"),1));
  nb=zeros(p,3); x=zeros(p,3); x2=zeros(p,3);
  betab=zeros(p,3);   betaw=zeros(p,3);  cbetab=zeros(c,3);  
  bbnds1=zeros(c,3);  bbnds2=zeros(c,3); paggs1=zeros(2,6);  
  paggs2=zeros(2,6);  abds1=zeros(2,6);  abds2=zeros(2,6);   
  goodman=zeros(2,6); double=zeros(2,3);

  /* read some parameters from the data buffer */
  cty = eiread(vread(ndbuf,"dbuf1"),"county");
  unicty = unique(cty,1);
  ctynm = eiread(vread(ndbuf,"dbuf1"),"cname");
  n = eiread(vread(ndbuf,"dbuf1"),"n");
  t = eiread(vread(ndbuf,"dbuf1"),"t");
  nt = eiread(vread(ndbuf,"dbuf1"), "Nt");

  /* loop for Hispanic and Whites */
  for i (1,2,1);
    dbuf1=vread(ndbuf,"dbuf"$+ftos(i*2-1,"*.*lf",1,0));

    /*** population calculation ***/
    nb[.,i] = eiread(dbuf1, "Nb");
    x[.,i] = eiread(dbuf1, "x");

    /*** betas ***/
    betab[.,i] = eiread(dbuf1, "betab");
    betaw[.,i] = eiread(dbuf1, "betaw");
    {vx, cbetab[.,i], nx} = cmeanwc(betab[.,i], cty, nb[.,i]); 
    bnds = eiread(dbuf1, "bounds");
    {vx, bbnds1[.,i], nx} = cmeanwc(bnds[.,1], cty, nb[.,i]);
    {vx, bbnds2[.,i], nx} = cmeanwc(bnds[.,2], cty, nb[.,i]);   
    paggs1[.,(i*2-1):i*2] = eiread(dbuf1, "paggs");
    abds1[.,(i*2-1):i*2] = eiread(dbuf1, "abounds");
    goodman[.,(i*2-1):i*2] = eiread(dbuf1, "goodman");
  endfor;

  /* Black from Hispanic and Whites */
  nb[.,3] = n-nb[.,1]-nb[.,2];
  x[.,3] = 1-x[.,1]-x[.,2];

  /*** Betas for Blacks ***/
  @ see eq 5.2 on p.80 and eq 15.4 on p.266 @
  {temp,temp1}=bounds1(t,x[.,3],n);
  bnds = temp[.,1:2];
  {nx, bbnds1[.,3], nx} = cmeanwc(bnds[.,1], cty, nb[.,3]);
  {nx, bbnds2[.,3], nx} = cmeanwc(bnds[.,2], cty, nb[.,3]);
  abds1[.,5:6]=temp1';
  bnds[.,1] = missrv(bnds[.,1],0);
  bnds[.,2] = missrv(bnds[.,2],1);
 
  betabs = (t-eiread(vread(ndbuf,"dbuf1"),"betabs").*x[.,1]
           -eiread(vread(ndbuf,"dbuf3"),"betabs").*(1-x[.,2]-x[.,3]))./x[.,3];
  for i (1,p,1);
    temp = delif(betabs[i,.]',(((betabs[i,.]'.<bnds[i,1])+(betabs[i,.]'.>bnds[i,2])
           +scalmiss(bnds[i,1])+scalmiss(bnds[i,2])).>0));
    betab[i,3]=meanc(temp);
    if rows(temp)/=100 and not(scalmiss(temp));
      temp = rndsmpl(temp,100);
      betabs[i,.]=temp';
    endif;
    if scalmiss(temp);
      for j (1,100,1);
        betabs[i,j]=miss(0,0);
      endfor;
    endif;
  endfor;
  
  for i (1,rows(unicty),1);
    if scalmiss(bbnds1[i,3]) or sumc(selif(nb[.,3],cty.==unicty[i]))/sumc(selif(n,cty.==unicty[i]))<0.005;
      cbetab[i,3]=miss(0,0);
    else;
      cbetab[i,3]= (sumc(selif(nt,cty.==unicty[i]))/sumc(selif(n,cty.==unicty[i]))
                  -cbetab[i,1]*sumc(selif(nb[.,1],cty.==unicty[i]))/sumc(selif(n,cty.==unicty[i]))
                  -cbetab[i,2]*sumc(selif(nb[.,2],cty.==unicty[i]))/sumc(selif(n,cty.==unicty[i])))
                  /(sumc(selif(nb[.,3],cty.==unicty[i]))/sumc(selif(n,cty.==unicty[i])));
      if (cbetab[i,3]>bbnds2[i,3] or cbetab[i,3]<bbnds1[i,3]);
        cbetab[i,3] = miss(1,1);
      endif;
    endif;
  endfor;

  @ see eq 15.6 on p.270 @
  paggs1[1,5] = (sumc(nt)/sumc(n)-paggs1[1,1]*sumc(nb[.,1])/sumc(n)-paggs1[1,3]*sumc(nb[.,2])/sumc(n))
                /(sumc(nb[.,3])/sumc(n));   

  paggs1[2,5] = sqrt(diag(vcx(meanwc(betabs,nb[.,3]))));
  for i (1,rows(unicty),1);
    if (sumc(nb[.,3])/sumc(n) < .05)or(paggs1[1,5]>abds1[2,5]) or (paggs1[1,5] < abds1[1,5]);
      paggs1[1,5]=miss(0,0);
      paggs1[2,5]=miss(0,0);
      cbetab[i,3]=miss(0,0);
    endif;
  endfor;

  {a,b,temp,temp,temp,temp}=reg(x[.,3]~(1-x[.,3]),t);
  goodman[.,5:6]=(a')|(b');

  if paggs1[1,5]==miss(0,0);
    paggs1[2,5]=miss(0,0);
  endif;  

  if rows(unicty)==1;
    paggs1[1,1 3 5]=cbetab;
  endif;

  format/rd 4,0;
  ?;
  "DEMOGRAPHIC STATISTICS AND ELECTION RESULTS";
  "Note: These statistics summarize the data"; 
  "      available for EI estimation.";
  ?;
  print "N " sumc(n);
  ?;
  format/rd 7,1;
  print "Percent Hispa:       " sumc(nb[.,1])/sumc(n)*100;
  print "Percent Black:       " sumc(nb[.,3])/sumc(n)*100;
  print "Percent White+Other: " sumc(nb[.,2])/sumc(n)*100;
  ?;
  format/rd 10,0;
  print "Hispans:       " sumc(nb[.,1]);
  print "Blacks:        " sumc(nb[.,3]);
  print "Whites+Other:  " sumc(nb[.,2]);
  ?;
  print "Number of Votes Cast: " sumc(nt);
  format/rd 7,1;
  print "Turnout %: " sumc(nt)/sumc(n)*100;
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
  temp=paggs1[.,1 5]~paggs1[.,3];
  temp;
  save temp;
  ?;
  "Aggregate Bounds";
  "        BETA-H  BETA-B  BETA-W";
  "Lower:" abds1[1,1 5]~abds1[1,3];
  "Upper:" abds1[2,1 5]~abds1[2,3];
  ?;
  "cnty#   County  betah betab betaw bh-lo bh-up bb-lo bb-up bw-lo bw-up";;
  mask = 1~0~ones(1,9);
  ?;
  fmt = {"*.*lf" 5 0, "*.*lG" 10 10, "*.*lf" 6 3, "*.*lf" 6 3,
  "*.*lf" 6 3, "*.*lf" 6 3, "*.*lf" 6 3, "*.*lf" 6 3, "*.*lf" 6 3,
  "*.*lf" 6 3, "*.*lf" 6 3};
  for k (1,rows(vx),1);
    temp = vx[k]~ctynm[vx[k]]~cbetab[k,1 3]~cbetab[k,2]~bbnds1[k,1]~bbnds2[k,1]~ 
           bbnds1[k,3]~bbnds2[k,3]~bbnds1[k,2]~bbnds2[k,2];
    prt = printfm(temp,mask,fmt);
    ?;
  endfor;
  ?;"Goodman's Regression";
  "  BETA-H  BETA-B  BETA-W";;
  goodman[.,1 5]~goodman[.,3];
  "*********************************************************************";
  _Eprt=eprt;
endp;
