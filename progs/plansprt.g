/*  
**   This procedure produces the output file for EI and
**   EI2 for different redistricting plans using the EI
**   estimates from the current district.  This is based 
**   on 3 seperate EI and EI2 runs, one for each of 3 
**   ethnic groups specified as X. This proc corresponds 
**   to txprint1.g
**   
**   Inputs:
**    distplan = redistricting plan indicator variable.
**            1=PLAN1001, 2=PLAN1021, 3=PLAN1025, 4=PLAN1034, 
**            5=PLAN1040, 6=PLAN1043, 7=PLAN1044, 8=PLAN1045, 
**            9=PLAN1046, 10=PLAN1047, 11=PLAN1048 
**    filenm = output file name.
**   
*/

proc (0)=plansprt(distplan,plannm);
  local p,c,dbuf,n,nb,nt,dv,betabs,cty,plans,name,ndbuf,t,x,ctynm,filenm;
  
  filenm=plannm$+".out";
  /* loop for each district */
  for i (1,32,1); 

    output file = ^filenm on;
    format/rd 4,0;
    ?;?;
    "**************************************************************";
    $plannm; 
    i "th District";
    "**************************************************************";
    format/rd 7,3;

    ndbuf="";
    /* loop for Hispanic and Whites */
    for j (1,2,1);
      /*** Betas ***/
      name="db"$+ftos(1+0,"*.*lf",1,0);
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
        name="db"$+ftos(k+0,"*.*lf",1,0);
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
    stwprint(ndbuf);
    output off;
  endfor;
endp;







