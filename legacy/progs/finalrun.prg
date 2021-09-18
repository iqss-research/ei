/*
**  A program for Texas election data
**  (C) Copyright 2001 Gary King
**  All Rights Reserved.
**
**  This program uses EI and EI2 to estimate Turnout and
**  Democrat vote for Hispanics, Blacks, and Whites in
**  each district of Texas for each election 1992-2000
**
**  Inputs: Data set requires precinct level data for:
**    County Name, district, Census data from 1990 and 2000 for
**    voting age population, hispanic v.a. population, white v.a.pop,
**    black v.a.pop & Election results coded as the number of
**    votes cast for the Democrats and Republicans.
**  VN abbreviation denotes 2000 Census data.
**  VA abbreviation denotes 1990 Census data.
**  USC92D: United States Congressional election 1992, Democratic votes.
**
**  Creates 2 output files:
**     results.txt computes estimates for Blacks from estimates of
**          Whites and Hispanics.
**     error.txt lists EI estimations where the Return code is non-zero.
**          EI Results in these districts are less reliable.
**
**  Outputs:
**     Betas - Generated from EI estimation.
**     Beta-H, Beta-B, Beta-W (Aggregate): District-level estimates of turnout
**          for Hispanics, Blacks, and Whites.  (And Standard Errors)
**     Aggregate Bounds for Beta-H, Beta-B, and Beta-W:
**          District-level Lowest and highest proportion of Hispanics,
**          Blacks, and Whites that could have voted.
**         ("Bounds" are the logical limits of how high and how low
**          the fractions of interest COULD be.  Estimates of BETAH, BETAB,
**          and BETAW must logically and always will fall between the high
**          and low bounds.)
**     County-level results:
**          Beta-H, Beta-B, and Beta-W and bounds for each county in the district.
**     Goodman's Regression: Results of Goodman's estimation technique.
**          Note that non-sensical results for Goodman's technique are due
**          to the technique, not our analyses.
**
**     Lambdas - Generated from EI2 Estimation.
**     Lambda-H, Lambda-B, Lambda-W: District-level estimates of the proportion
**          of Hispanics, Balcks, and Whites that voted Democrat.
**     Aggregate Bounds for Lambda-H, Lambda-B, Lambda-W:
**          District-level lowest and highest proportions of Hispanics,
**          Blacks, and Whites that could have voted Democrat (see note above).
**     County-level results:
**          Lambda-H, Lambda-B, Lambda-W and bounds for each county in the district.
**     Double Regression: results of Double Regression Technique.
**          Note that non-sensical results here are due to the problems
**          of the Double Regression technique, not due to the data or
**          programming error.
**
**  Procedures needed to run this program:
**    cln.g - eliminates missing data.
**    intpl.g - interpolates data for 1992, 1994, 1996, and 1998 from Census
**        figures of 1990 and 2000.
**    globals.g - for specifying 4 globals for EI.
**    txprint.g - Produces file of results with Estimates for Blacks derived
**        from estimates of Whitess and Hispanics.
**    minr1.g. - needed to generate bounds in txprint2.
**/

new;
/*** inputs to be specified ***/
ds = "txdata5e";             @ input dataset name @

otpfile = "result.txt";      @ output file name @
errfile = "reserr.txt";   @ error file name @
yrstart = 3;              @ Start with year#____: 1=1992, 2=1994, 3=1996, 4=1998, 5=2000 @
n_year = 3;                @ last election year @
diststar = 3;              @ Start with District#____@
n_dist = 30;               @ # of districts in Texas @

/*** clearing variables ***/
clear plan, vtdkey00, county, cnty_nam, plan_cd, IS_SPLIT, PCTOFVTD, cd103, cd107, 
PRS92D, PRS92R, PRS92I, USC92D, USC92R, RRC92D, RRC92R, SCP3_92D, 
SCP3_92R, CCAP4_92, CCAP4_91, STSEN92D, STHSE92D, STSEN92R, STHSE92R, 
USN94D, USN94R, USC94D, USC94R, GOV94D, GOV94R, LG94D, LG94R, ATG94D, 
ATG94R, COMP94D, COMP94R, TRS94D, TRS94R, LC94D, LC94R, AGRI94D,AGRI94R, 
RRC94D, RRC94R,SCP3_94D, SCP3_94R, CCAPP_94, CCAP7_94, CCAP7_91, CCAP8_94, 
CCAP8_91, STSEN94D, STSEN94R, STHSE94D, STHSE94R, PRS96R, PRS96D, PRS96I, 
USN96R, USN96D, RRC96R, RRC96D, SCCJ_96R, SCCJ_96D, SCP3_96R, SCP3_96D, 
CCAP4_96, CCAP4_92, STSEN96R, STSEN96D,STHSE96R, STHSE96D, USC96R, USC96D, 
usc96r2, usc96d2,  @ <------------ added to fix 96 redistricting problem @
USC98R, USC98D, GOV98R, GOV98D,LG98R, LG98D, ATG98R, ATG98D, COMP98R,
COMP98D,LC98R,LC98D,AGRI98R,AGRI98D,RRC98R,RRC98D,SCP3_98R,SCP3_98D,CCAP4_98,
CCAP4_91,STHSE98R,STSEN98R,STHSE98D,STSEN98D,STHSE98,PRS00R,PRS00D,
USN00R,USN00D,USC00R,USC00D,RRC00R,SCP3_00R,CCAPP_00,CCAPP_01,
CCAP7_00,CCAP8_00,CCAP8_01,STSEN00R,STSEN00D,STHSE00R,STHSE00D,TAPERSON,
TA1RACE,TAWHITE,TABLACK,TAAMERIN,TAASIAN,TAHAWPAC,TA1OTHER,TA2MRACE,
TNTOTAL,TNHISPAN,TNNONHIS,TN1RACE,TNWHITE,TNBLACK,TNAMINAK,TNASIAN,
TNHAWOP,TNOTHER,TN2MRACE,VAPERSON,VA1RACE,VAWHITE,VABLACK,VAAMERIN,
VAASIAN,VAHAWPAC,VA1OTHER,VA2MRACE,VNTOTAL,VNHISPAN,VNNONHIS,VN1RACE,
VNWHITE,VNBLACK,VNAMINAK,VNASIAN,VNHAWOP,VNOTHER,VN2MRACE,POP90,
W90,B90,NA90,AP90,O90,VAW90,VAB90,VANA90,VAAP90,VAO90,H90,NHW90,
NHB90,NHNA90,NHAP90,NHO90,VAH90,VANHW90,VANHB90,VANHNA90,
VANHAP90,VANHO90,VAPOP90; 

/*** reading in variables ***/
call subdatv(ds,0);
nobs=rows(vtdkey00);

vtdkey00[1:100];


end;

/* renaming the variables */
VAHISP90 = VAH90;
VANHWHIT = VANHW90;
VANHBLAC = VANHB90;
VAPOP90 = VAW90+VAB90+VANA90+VAAP90+VAO90;
@VAPOP90 = VAH90+VAW90+VAB90+VANA90+VAAP90+VAO90;@

/*** create a new county var, n_county ***/
ctynm = unique(cnty_nam, 0);   @ county name @
n_county = ones(nobs, 1);    @ county number @

@ works if counties are out of order @
ctyi=seqa(1,1,rows(ctynm)); 
for i (1,nobs,1);
  n_county[i]=selif(ctyi,ctynm.$==cnty_nam[i]);
endfor;

/*** Creating EI inputs ***/
/* N's:  voting age population through interpolation */
nttl = intpl(VAPOP90,VNTOTAL,4);              @ total @
nhis = intpl(VAHISP90,VNHISPAN,4);            @ hispanics @
nnonhis = intpl(VAPOP90-VAHISP90,VNNONHIS,4); @ non-hispanics @
nbla = intpl(VANHBLAC,VNBLACK,4);             @ black @
nnonbla = intpl(VAPOP90-VANHBLAC,VNTOTAL-VNBLACK,4); @ non-black @
@ whites + others @
nwhi = intpl(VAPOP90-VAHISP90-VANHBLAC,VNTOTAL-VNHISPAN-VNBLACK,4);
nnonwhi = nhis+nbla;
@ whites only @
@nwhi = intpl(VANHWHIT,VNWHITE,4);@             
@nnonwhi = intpl(VAPOP90-VANHWHIT,VNTOTAL-VNWHITE,4);@  @ non-whites @

nthis=nhis+nnonhis;
ntbla=nbla+nnonbla;
ntwhi=nwhi+nnonwhi;

@ X's: ratio of racial groups @
xhis = nhis./nthis;
xbla = nbla./ntbla;
xwhi = nwhi./ntwhi;

end;

@ ** fix 96 redistricting problem @
cd103=stof(cd103);
cd107=stof(cd107);
ks=seqa(1,1,nobs);
let kchgd = 3 5 6 7 8 9 18 22 24 25 26 29 30;
kt=sumc((cd107.==kchgd')');
kchg=selif(ks,kt); @ indices of changed districts @
usc96d[kchg]=usc96d2[kchg];
usc96r[kchg]=usc96r2[kchg];

@ T's: turnout @
VOTES = (USC92D+USC92R)~(USC94D+USC94R)~(USC96D+USC96R)~(USC98D+USC98R)~(USC00D+USC00R);
T = VOTES./nttl[.,2:6];
this= VOTES./nthis[.,2:6];
tbla= VOTES./ntbla[.,2:6];
twhi= VOTES./ntwhi[.,2:6];

@ V's: Democratic vote share @
DEMV = USC92D~USC94D~USC96D~USC98D~USC00D;
V = DEMV./VOTES;

/*** EI inputs for His vs. NonHis, Bla vs. NonBla, and Whi vs. NonWhi ***/
hisdta = this~nthis[.,2:6]~xhis[.,2:6]~V~cd103~n_county~cd107;
hisdta = delif(hisdta, hisdta[.,21].==0); @ drop cd103 = 0 @
whidta = twhi~ntwhi[.,2:6]~xwhi[.,2:6]~V~cd103~n_county~cd107;
whidta = delif(whidta, whidta[.,21].==0);

/*** run ei for each district in every year ***/
library @/tmp/ei/p/@ei;

for i (yrstart, n_year, 1); @ loop for election @
   hisdta1 = cln(hisdta[.,i i+5 i+10 i+15 21:23]); @ i think 21 begins with cd @
   whidta1 = cln(whidta[.,i i+5 i+10 i+15 21:23]);

   /*** initialize betabs for JudgeIt ***/
   betabs=zeros(1,5);

   /* if before 96, then use cd103. otherwise use cd107 */
   if i<3; col=5; 
   else; col=7; 
   endif;

   for j (diststar, n_dist, 1); @ loop for district @
      output file = ^otpfile on;
      format/rd 4,0;
      "**************************************************************";
      "Year:" i*2+90 j "th District";
      "**************************************************************";
      format/rd 7,3;
      output off;

      hisdta2 = selif(hisdta1, hisdta1[.,col].==j);
      whidta2 = selif(whidta1, whidta1[.,col].==j);

      /** skip uncontested districts  ***/
      demshr=meanwc(hisdta2[.,4],hisdta2[.,2]);
      if ((0.05 < demshr) and (demshr < 0.95));


         @EI HISPANICS VS. NON-HISPANICS@
         eiset;
         _Eeta=3;
         _EalphaB=0~0.2;
         _EalphaW=0~0.2;
         _Eres=vput(_Eres,"Hispanics(EI):"$+ftos(((90+i*2)*100+j),"*.*lf",4,0),"titl");
         _Eres = vput(_Eres, hisdta2[.,6], "county");
         _Eres = vput(_Eres, ctynm, "cname");
         _Erho=.1;
	 _EisFac=-2;
         _Ebounds=(-5~5)|(-2~2)|(-5~5)|(-2~2)|(-4~3)|(-4~3)|(-2~2);
         dbuf1 = ei(hisdta2[.,1], hisdta2[.,3], hisdta2[.,2], 1, 1);
         ndbuf=vput("",dbuf1,"dbuf1");


         eiset;
         _Eeta=3;
         _EalphaB=0~0.2;
         _EalphaW=0~0.2;
         _Ei2_m = -1;
         _Eres = vput(_Eres,"Hispanics(EI2):"$+ftos(((90+i*2)*100+j),"*.*lf",4,0),"titl");
         _Erho=.1;
	 _EisFac=-2;
         _Ebounds=(-5~5)|(-2~2)|(-5~5)|(-2~2)|(-4~3)|(-4~3)|(-2~2);
         dbuf2 = ei2(hisdta2[.,4], dbuf1, 1, 1);
         ndbuf=vput(ndbuf,dbuf2,"dbuf2");


         @EI WHITES VS. NON-WHITES@
         eiset;
         _Eeta=3;
         _EalphaB=0~0.2;
         _EalphaW=0~0.2;
         _Eres = vput(_Eres,"Whites(EI):"$+ftos(((90+i*2)*100+j),"*.*lf",4,0),"titl");
         _Eres = vput(_Eres, whidta2[.,6], "county");
         _Eres = vput(_Eres, ctynm, "cname");
         _Erho=.1;
	 _EisFac=-2;
         _Ebounds=(-5~5)|(-2~2)|(-5~5)|(-2~2)|(-4~3)|(-4~3)|(-2~2);
         dbuf3 = ei(whidta2[.,1], whidta2[.,3], whidta2[.,2], 1, 1);
         ndbuf=vput(ndbuf,dbuf3,"dbuf3");


         eiset;
         _Eeta=3;
         _EalphaB=0~0.2;
         _EalphaW=0~0.2;
         _Ei2_m = -1;
         _Eres = vput(_Eres,"Whites(EI2):"$+ftos(((90+i*2)*100+j),"*.*lf",4,0),"titl");
         _Erho=.1;
         _Ebounds=(-5~5)|(-2~2)|(-5~5)|(-2~2)|(-4~3)|(-4~3)|(-2~2);
	 _EisFac=-2;
         dbuf4 = ei2(whidta2[.,4], dbuf3, 1, 1);
         ndbuf=vput(ndbuf,dbuf4,"dbuf4");
         
         

         name="db"$+ftos(((90+i*2)*100+j),"*.*lf",4,0);
         save ^name=ndbuf;

         /*** Print Results and Error File ***/
         output file = ^otpfile on;
         txprint(ndbuf);
         output off;
         loadm temp;
         betabs=betabs|(i~j~temp[1,.]);
         for k (1,4,1);
           dbuf=vread(ndbuf,"dbuf"$+ftos(k+0,"*.*lf",1,0));
           if eiread(dbuf,"retcode")/=0;
             output file = ^errfile on;
             format/rd 4,0;
             "Year:" i*2+90 j "th District"; ?;
             "WARNING: retcode is not zero for "$+eiread(dbuf,"titl");
             "         retcode =" eiread(dbuf,"retcode");
             output off;
           endif;
         endfor;
      else;
         output file = ^otpfile on; ?;
         "This district is not contested. No EI estimation is conducted."; ?;
         output off;
      endif;
   endfor;
   name="betas"$+ftos(((90+i*2)*100),"*.*lf",2,0);
   betabs=betabs[2:rows(betabs),.];
   save ^name=betabs;
endfor;





