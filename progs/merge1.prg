/* 
** This program merges the primary data set with redistricting 
** plan data sets using a key variable common in all the data 
** sets.
**
** Note: The program will ignore the precincts in the plan
**       data sets which do not appear in the primary data set.
**
*/


new;
/*** inputs to be specified ***/
prmrds = "0010a";                @ primary dataset name @
let plands = plan1001 plan1021 plan1025 plan1034 plan1040 plan1043 
             plan1044 plan1045 plan1046 plan1047 plan1048;          
             @ plan data sets @
let keyvar = VTDKEY01 VTDKEY02; @ 2 key variabels @
needvar = "plan_cd";            @ variables to merge @
newdata = "0010b";              @ new data set name @

/*** main program ***/
varname1=getname(prmrds);
mstrdta=subdatd(prmrds,varname1);
kyvar1a=subdatd(prmrds,keyvar[1]);
kyvar2a=subdatd(prmrds,keyvar[2]);
nobs1=rows(kyvar1a);


for i(1,rows(plands),1);

  mergeds=subdatd(plands[i],needvar);
  kyvar1b=subdatd(plands[i],keyvar[1]);
  kyvar2b=subdatd(plands[i],keyvar[2]);
  nobs2=rows(kyvar1b);

  /*** merging the two data sets ***/
  temp=zeros(nobs1,rows(needvar));
  indx=zeros(nobs2,1);

  for j(1,nobs1,1);
    "Now merging data " $prmrds "and" $plands[i];
    k=1;
    do while (kyvar1a[j]$/=kyvar1b[k] or kyvar2a[j]$/=kyvar2b[k]) and k<nobs2;
      k=k+1;
    endo;
    if kyvar1a[j]$==kyvar1b[k] and kyvar2a[j]$==kyvar2b[k];
      j "th obs matched with" k "th obs!";
      temp[j,.]=mergeds[k,.];
      indx[k]=1;
    else;
      j "th obs didn't match!";
      for l(1,cols(temp),1);
        temp[j,l]=miss(0,0);
      endfor;
    endif;
  endfor;
  mstrdta=mstrdta~temp;
  saved(mstrdta,newdata,varname1|plands[1:i]);
endfor;
