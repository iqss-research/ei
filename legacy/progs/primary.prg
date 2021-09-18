/*
**  A program for Texas election data
**  (C) Copyright 2001 Gary King
**  All Rights Reserved.
**
**  This program uses EI to estimate Turnout and
**  voting by race of candidate and voter in
**  Texas primaries between 1992-2000
**
**  Inputs: Data set requires precinct level data for:
**    County Name, district, Census data from 1990 and 2000 for
**    voting age population, hispanic v.a. population, white v.a.pop,
**    black v.a.pop & Election results coded as the number of
**    votes cast for the Democrats and Republicans.
**  VN abbreviation denotes 2000 Census data.
**  VA abbreviation denotes 1990 Census data.
**  Votes for Hispanic, Black, and White candidates in the primary
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
**    txbprint.g - Produces file of results with Estimates for Blacks derived
**                 from estimates of Blacks and Hispanics.
**    minr1.g. - needed to generate bounds in txprint2.
**/

new;
library pgraph, cml;
library /home/kimai/ei/p/ei;

/*** inputs to be specified ***/
ds = "9810b";          @ input dataset name @
otpfile = "ac98p.out";    @ output file name @
errfile = "ac98p.err";    @ error file name @

/*** clearing variables ***/ 
clear vtdkey00, vtdkey01, vtdkey02, DISTRICT, CNTYKEY, E_VAP,
E_ANGVAP, E_BLAKVA, E_HSPVAP, CD103, CD107, CNTY_NAM, VNTOTAL,
VNHISPAN, VNNONHIS, VNWHITE, VNBLACK, VNAMINAK, VNASIAN, VNHAWOP,
VNOTHER, POP90, W90, B90, NA90, AP90, O90, VAW90, VAB90, VANA90,
VAAP90, VAO90, H90, NHW90, NHB90, NHNA90, NHAP90, NHO90, VAH90,
VANHW90, VANHB90, VANHNA90, VANHAP90, VANHO90, PLAN1001, PLAN1021, 
PLAN1025, PLAN1034, PLAN1040, PLAN1043, PLAN1044, PLAN1045, PLAN1046, 
PLAN1047, PLAN1048, MATTOX, OVERSTRE, KELLY, PATTERSO, DE_LEON;

/*** reading in variables ***/
call subdatv(ds,0);

/* Year read in from primary datafile  */
yr = 4;                    @ year: 1=92,2=94,3=96,4=98,5=00 @
diststar = 30;              @ Start with District#____@
n_dist = 30;             @ # of districts in Texas @

/* renaming the variables */
VOTES = PATTERSO+DE_LEON;
VNTOTAL=E_VAP;
VNHISPAN=E_HSPVAP;
VNNONHIS=VNTOTAL-VNHISPAN;
VNBLACK=E_BLAKVA;
VNNONHIS = VNTOTAL - VNHISPAN;
VAHISP90 = VAH90;
VANHWHIT = VANHW90;
VANHBLAC = VANHB90;
VAPOP90 = VAH90+VANHW90+VANHB90+VANHNA90+VANHAP90+VANHO90;

nobs=rows(vtdkey00);

/*** create a new county var, n_county ***/
ctynm = unique(cnty_nam, 0);        @ county name @
n_county = ones(nobs, 1);           @ county number @

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
/* whites only 
nwhi = intpl(VANHWHIT,VNWHITE,4);             
nnonwhi = intpl(VAPOP90-VANHWHIT,VNTOTAL-VNWHITE,4); */

nthis=nhis+nnonhis;
ntbla=nbla+nnonbla;
ntwhi=nwhi+nnonwhi;

@ X's: ratio of racial groups @
xhis = nhis./nthis;
xbla = nbla./ntbla;
xwhi = nwhi./ntwhi;

@ T's: turnout @
T = VOTES./nttl[.,2:6];
this= VOTES./nthis[.,2:6];
tbla= VOTES./ntbla[.,2:6];
twhi= VOTES./ntwhi[.,2:6];

@ plans @
plans=PLAN1001~PLAN1021~PLAN1025~PLAN1034~PLAN1040~PLAN1043~PLAN1044~PLAN1045~PLAN1046~PLAN1047~PLAN1048;


/*** EI inputs for His vs. NonHis, Bla vs. NonBla, and Whi vs. NonWhi ***/
hisdta = this[.,yr]~nthis[.,yr+1]~xhis[.,yr+1]~n_county~district~stof(plans);
whidta = twhi[.,yr]~ntwhi[.,yr+1]~xwhi[.,yr+1]~n_county~district~stof(plans);

/*** run ei for each district in every year ***/
/*** initialize betabs for JudgeIt ***/
betabs=zeros(1,5);

for j (diststar,n_dist,1); @ loop for district @

    output file = ^otpfile on;
    format/rd 4,0;
    "**************************************************************";
    j "th District";
    "**************************************************************";
    format/rd 7,3;
    output off;

    hisdta1=selif(hisdta,hisdta[.,5].==j);
    whidta1=selif(whidta,whidta[.,5].==j);

    /** do the same thing as cln.g **/  
    hisdta1 = packr(hisdta1);
    whidta1 = packr(whidta1);
    hisdta2 = delif(hisdta1,(hisdta1[.,1].>1).or(hisdta1[.,1].<=0).or(hisdta1[.,2].<= 0) 
                   .or(hisdta1[.,3].>1).or(hisdta1[.,3].< 0).or(whidta1[.,1].>1).or(whidta1[.,1].<=0).or(whidta1[.,2].<= 0) 
                   .or(whidta1[.,3].>1).or(whidta1[.,3].< 0));
    whidta2 = delif(whidta1,(hisdta1[.,1].>1).or(hisdta1[.,1].<=0).or(hisdta1[.,2].<= 0) 
                   .or(hisdta1[.,3].>1).or(hisdta1[.,3].< 0).or(whidta1[.,1].>1).or(whidta1[.,1].<=0).or(whidta1[.,2].<= 0) 
                   .or(whidta1[.,3].>1).or(whidta1[.,3].< 0));


    @EI HISPANICS VS. NON-HISPANICS @
    eiset;
    _Eeta=3;
    _EalphaB=0~.3;
    _EalphaW=0~.3;
    _Eres = vput(_Eres,"Hispanics(EI): "$+ftos(j+0,"*.*lf",2,0),"titl");
    _Eres = vput(_Eres,hisdta2[.,4],"county");
    _Eres = vput(_Eres,hisdta2[.,6:16],"plans");
    _Eres = vput(_Eres,ctynm,"cname");
    _Erho=1;
    _EisFac=-2;
    _Ebounds=(-5~5)|(-2~2)|(-5~5)|(-2~2)|(-4~3)|(-4~3)|(-2~2);
    dbuf1 = ei(hisdta2[.,1],hisdta2[.,3],hisdta2[.,2],1,1);
    ndbuf=vput("",dbuf1,"dbuf1");

    @EI WHITES VS. NON-WHITES @
    eiset;
    _Eeta=3;
    _EalphaB=0~.3;
    _EalphaW=0~.3;
    _Eres = vput(_Eres,"Whites(EI):"$+ftos(j+0,"*.*lf",2,0),"titl");
    _Eres = vput(_Eres,whidta2[.,4],"county");
    _Eres = vput(_Eres,whidta2[.,6:16],"plans");
    _Eres = vput(_Eres,ctynm,"cname");
    _Erho=1;
    _EisFac=-2;
    _Ebounds=(-5~5)|(-2~2)|(-5~5)|(-2~2)|(-4~3)|(-4~3)|(-2~2);
    dbuf3 = ei(whidta2[.,1],whidta2[.,3],whidta2[.,2],1,1);
    ndbuf=vput(ndbuf,dbuf3,"dbuf3");

    name="db"$+ftos(j+0,"*.*lf",1,0),;
    save ^name=ndbuf; 

    /*** Print Results and Error File ***/
    output file = ^otpfile on;
    txbprint(ndbuf);
    output off;
    loadm temp;
    betabs=betabs|(yr~j~temp[1,.]);
    for k (1,3,2);
      dbuf=vread(ndbuf,"dbuf"$+ftos(k+0,"*.*lf",1,0));
      if eiread(dbuf,"retcode")/=0;
        output file = ^errfile on;
        format/rd 4,0;
        j "th District"; ?;
        "WARNING: retcode is not zero for "$+eiread(dbuf,"titl");
        "         retcode =" eiread(dbuf,"retcode");
        output off;
      endif;
    endfor;
endfor;
name="betas";
betabs=betabs[2:rows(betabs),.];
save ^name=betabs;

