/*
**
    compelct(dists,filenm);
**
**  This proc reaggregates the EI turnout results based on 
**  different redistricting plans. This proc uses another 
**  proc, cmp.g
**
**
**  INPUT
**     dists = redistricting plan indicator variable.
**          1=PLAN1001C, 2=PLAN1021C, 3=PLAN1025C, 4=PLAN1034C, 
**          5=PLAN1040C, 6=PLAN1043C, 7=PLAN1044C, 8=PLAN1045C, 
**          9=PLAN1046C, 10=PLAN1047C, 11=PLAN1048C 
**     filenm = output file name.
**
**  OUTPUT
**     It prints out the results in specified file.
**
*/

proc 0=compelct(dists,filenm);
  local p,c,dbuf,n,nb,nt,dv,betabs,cty,plans,name,ndbuf,t,x,
        ctynm,i,upbnds,temp,temp1,temp2,mask,wgts,
        fmt,prt,a,b,lwbnds,plan_nm;
  
  let plan_nm= P1001C, P1021C, P1025C, P1034C, 
               P1040C, P1043C, P1044C, P1045C, 
               P1046C, P1047C, P1048C, P1000C; 
  _Eprt=0;
  output file = ^filenm on;
  "               ESTIMATES";
  "plan#  dist# HIS   BLA   WHI";;
  output off;
  
for distplan(1,rows(dists),1); @ loop for plans @
  temp2=zeros(1,3);
  lwbnds=zeros(1,3);
  upbnds=zeros(1,3);
  wgts=0;
  
  for h (1,2,1);  @ loop for primary races @  
    i=dists[distplan];
    
    if distplan<rows(dists);
    ndbuf="";
    /* loop for Hispanic and Whites */
    for j (1,2,1);
     /*** Betas ***/
      if h==1;
        name="c:/gauss36/kosuke/ac_98p/db"$+ftos(1+0,"*.*lf",1,0);
      elseif h==2;
        name="c:/gauss36/kosuke/ag_98p/db"$+ftos(1+0,"*.*lf",1,0);
      elseif h==3;
        name="c:/gauss36/kosuke/pr_96p/db"$+ftos(1+0,"*.*lf",1,0);
      elseif h==4;
        name="c:/gauss36/kosuke/s_96p/db"$+ftos(1+0,"*.*lf",1,0);
      elseif h==5;
        name="c:/gauss36/kosuke/s_00p/db"$+ftos(1+0,"*.*lf",1,0);
      elseif h==6;
        name="c:/gauss36/kosuke/pr_00p/db"$+ftos(1+0,"*.*lf",1,0);
      else;
        name="c:/gauss36/kosuke/s_96r/db"$+ftos(1+0,"*.*lf",1,0);
      endif;
    
      loadm dbuf=^name;
      dbuf=vread(dbuf,"dbuf"$+ftos(j*2-1,"*.*lf",1,0));
      betabs=eiread(dbuf,"betabs");  
      t=eiread(dbuf,"t");  
      x=eiread(dbuf,"x");  
      n=eiread(dbuf,"n");  
      cty=eiread(dbuf,"county");
      ctynm=eiread(dbuf,"cname");
      plans=eiread(dbuf,"plans");

      for k (2,30,1);
        if h==1;
          name="c:/gauss36/kosuke/ac_98p/db"$+ftos(k+0,"*.*lf",1,0);
        elseif h==2;
          name="c:/gauss36/kosuke/ag_98p/db"$+ftos(k+0,"*.*lf",1,0);
        elseif h==3;
          name="c:/gauss36/kosuke/pr_96p/db"$+ftos(k+0,"*.*lf",1,0);
        elseif h==4;
          name="c:/gauss36/kosuke/s_96p/db"$+ftos(k+0,"*.*lf",1,0);
        elseif h==5;
          name="c:/gauss36/kosuke/s_00p/db"$+ftos(k+0,"*.*lf",1,0);
        elseif h==6;
          name="c:/gauss36/kosuke/pr_00p/db"$+ftos(k+0,"*.*lf",1,0);
        elseif h==7;
          name="c:/gauss36/kosuke/s_96r/db"$+ftos(k+0,"*.*lf",1,0);
        endif;

        loadm dbuf=^name;
        dbuf=vread(dbuf,"dbuf"$+ftos(j*2-1,"*.*lf",1,0));
        betabs=betabs|eiread(dbuf,"betabs");
        t=t|eiread(dbuf,"t");
        x=x|eiread(dbuf,"x");
        n=n|eiread(dbuf,"n");
        cty=cty|eiread(dbuf,"county");
        plans=plans|eiread(dbuf,"plans");
      endfor;

      dbuf=vput("",selif(betabs,plans[.,distplan].==i),"betabs");
      dbuf=vput(dbuf,selif(t,plans[.,distplan].==i),"t");
      dbuf=vput(dbuf,selif(x,plans[.,distplan].==i),"x");
      dbuf=vput(dbuf,selif(n,plans[.,distplan].==i),"n");
      dbuf=vput(dbuf,selif(cty,plans[.,distplan].==i),"county");
      dbuf=vput(dbuf,ctynm,"cname");
      dbuf=vput(dbuf,selif(plans,plans[.,distplan].==i),"plans");

      ndbuf=vput(ndbuf,dbuf,"dbuf"$+ftos(j*2-1,"*.*lf",1,0));
    endfor;
 
    else;
        if h==1;
          name="c:/gauss36/kosuke/ac_98p/db"$+ftos(i+0,"*.*lf",1,0);
        elseif h==2;
          name="c:/gauss36/kosuke/ag_98p/db"$+ftos(i+0,"*.*lf",1,0);
        elseif h==3;
          name="c:/gauss36/kosuke/pr_96p/db"$+ftos(i+0,"*.*lf",1,0);
        elseif h==4;
          name="c:/gauss36/kosuke/s_96p/db"$+ftos(i+0,"*.*lf",1,0);
        elseif h==5;
          name="c:/gauss36/kosuke/s_00p/db"$+ftos(i+0,"*.*lf",1,0);
        elseif h==6;
          name="c:/gauss36/kosuke/pr_00p/db"$+ftos(i+0,"*.*lf",1,0);
        elseif h==7;
          name="c:/gauss36/kosuke/s_96r/db"$+ftos(i+0,"*.*lf",1,0);  
        endif;  
        loadm ndbuf=^name;
    endif;
    
    {a,b,c}=cmp(ndbuf);
    temp2=temp2|a[1,.];
    upbnds=upbnds|b[1,.];
    lwbnds=lwbnds|b[2,.];
    wgts=wgts|c;
  endfor;

  temp2=temp2[2:rows(temp2),.];
  upbnds=upbnds[2:rows(upbnds),.];
  lwbnds=lwbnds[2:rows(lwbnds),.];
  wgts=wgts[2:rows(wgts)];
  ?;
  temp2=temp2~ones(2,1)~wgts;
  temp2;
  @ temp2=selif(temp2,ones(2,1)|zeros(1,1)|ones(3,1)); @
  @ temp2=selif(temp2,zeros(1,1)|ones(1,1)|zeros(4,1)); @
  @ temp2=packr(temp2); @
  rows(temp2);
  a=zeros(1,3);
  b=zeros(1,3);
  c=zeros(1,3);
  for k (1,3,1);
    {temp,a[k],temp1}=cmeanwc(temp2[.,k],temp2[.,4],temp2[.,5]);
    {temp,b[k],temp1}=cmeanwc(upbnds[.,k],ones(6,1),wgts);
    {temp,c[k],temp1}=cmeanwc(lwbnds[.,k],ones(6,1),wgts);
  endfor;
  sumc(a');
  
  output file = ^filenm on;  
  /** if you want bounds **/
  mask = 0~ones(1,10);
  ?;
  fmt = {"*.*lG" 6 6, "*.*lf" 6 0, "*.*lf" 6 3, "*.*lf" 6 3,"*.*lf" 6 3, 
  "*.*lf" 6 3, "*.*lf" 6 3, "*.*lf" 6 3, "*.*lf" 6 3, "*.*lf" 6 3, "*.*lf" 6 3};
  temp = plan_nm[distplan]~i~a~b~c;
  ?;
  /** if you don't want bounds **/
  mask = 0~ones(1,4);
  fmt = {"*.*lG" 6 6, "*.*lf" 6 0, "*.*lf" 6 3, "*.*lf" 6 3, "*.*lf" 6 3};
  temp = plan_nm[distplan]~i~a;
  prt = printfm(temp,mask,fmt);
  output off;
endfor;
endp;


/*
** Support proc for compelct.g 
**
*/


proc (3)=cmp(ndbuf);
 local p,c,dbuf1,dbuf2,cty,n,n2,nb,t,t2,x,x2,nt,v,d,dv,betab,betaw,cbetab,
       nx,vx,bnds,bbnds1,bbnds2,lbnds1,lbnds2,lambdab,lambdaw,clambdab,
       paggs1,paggs2,abds1,betabs,betaws,abds2,temp,temp1,mask,fmt,prt,
       eprt,ctynm,lambdas,lambdabs,lambdaws,goodman,double,a,b,unicty;
  
  /* initialize the parameters */
  p=rows(eiread(vread(ndbuf,"dbuf1"),"betab"));
  c=rows(unique(eiread(vread(ndbuf,"dbuf1"),"county"),1));
  nb=zeros(p,3); x=zeros(p,3); x2=zeros(p,3);
  betab=zeros(p,3);   betaw=zeros(p,3);
  lambdab=zeros(p,3); lambdaw=zeros(p,3);
  cbetab=zeros(c,3);  clambdab=zeros(c,3);
  bbnds1=zeros(c,3);  bbnds2=zeros(c,3);
  lbnds1=zeros(c,3);  lbnds2=zeros(c,3);
  paggs1=zeros(2,6);  paggs2=zeros(2,6);
  abds1=zeros(2,6);   abds2=zeros(2,6);
  goodman=zeros(2,6); double=zeros(2,3);

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
    temp = eiread(dbuf1,"aggs");
    paggs2[2,(i*2-1)] = stdc(temp[.,1].*sumc(nb[.,i])/sumc(nt)); 
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
        cbetab[i,3] = miss(0,0);
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


  paggs2[1,1 3 5]=paggs1[1,1 3 5].*sumc(nb)'/sumc(nt);
  paggs2[2,5]=sqrt(diag(vcx(meanwc(betabs,nb[.,3])*sumc(nb[.,3])/sumc(nt))));
  if scalmiss(paggs2[1,5]);
    paggs2[2,5]=miss(0,0);
  endif;
 
  {a,b,temp,temp,temp,temp}=reg(x[.,3]~(1-x[.,3]),t);
  goodman[.,5:6]=(a')|(b');

  if paggs1[1,5]==miss(0,0);
    paggs1[2,5]=miss(0,0);
  endif;  

  if rows(unicty)==1;
    paggs1[1,1 3 5]=cbetab;
  endif;

  temp=paggs2[.,1 5]~paggs2[.,3];
  temp1=abds1[1,1]*sumc(nb[.,1])/sumc(nt)~abds1[1,5]*sumc(nb[.,3])/sumc(nt)~abds1[1,3]*sumc(nb[.,2])/sumc(nt);
  temp1=temp1|(abds1[2,1]*sumc(nb[.,1])/sumc(nt)~abds1[2,5]*sumc(nb[.,3])/sumc(nt)~abds1[2,3]*sumc(nb[.,2])/sumc(nt));
  retp(temp,temp1,sumc(nt));
endp;







