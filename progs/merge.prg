/* 
** This program merges the primary data set with demographic data
** set using a key variable common in the two data sets.
**
** Note: The program will ignore the precincts in the demographic
**       data set which do not appear in the primary data set.
**
*/


new;

/*** inputs to be specified ***/
prmrds = "0010";                 @ primary dataset name @
let plands = txdata5e;           @ demog data sets @
let keyvar = VTDKEY01 VTDKEY02;  @ 2 key variabels @
let needvar = CNTY_NAM CD103 CD107 VNTOTAL VNHISPAN VNNONHIS
VNWHITE VNBLACK VNAMINAK VNASIAN VNHAWOP VNOTHER POP90 W90
B90 NA90 AP90 O90 VAW90 VAB90 VANA90 VAAP90 VAO90 H90 NHW90
NHB90 NHNA90 NHAP90 NHO90 VAH90 VANHW90 VANHB90 VANHNA90
VANHAP90 VANHO90;            @ variables to merge from demog data set @
newdata = "0010a";           @ new dataset name @

/*** main program ***/
varname1=getname(prmrds);
mstrdta=subdatd(prmrds,varname1);
kyvar1a=subdatd(prmrds,keyvar[1]);
kyvar2a=subdatd(prmrds,keyvar[2]);
nobs1=rows(kyvar1a);

mergeds=subdatd(plands,needvar);
kyvar1b=subdatd(plands,keyvar[1]);
kyvar2b=subdatd(plands,keyvar[2]);
nobs2=rows(kyvar1b);

/*** merging the two data sets ***/
temp=zeros(nobs1,rows(needvar));
indx=zeros(nobs2,1);
for j(1,nobs1,1);
  k=1;
  do while (kyvar1a[j]$/=kyvar1b[k] or kyvar2a[j]$/=kyvar2b[k]) and k<nobs2;
    k=k+1;
  endo;
  if kyvar1a[j]$==kyvar1b[k] and kyvar2a[j]$==kyvar2b[k];
    j "th obs matched!";
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

saved(mstrdta,newdata,varname1|needvar);
