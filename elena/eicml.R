##/*
##  This archive is part of the program EI
##  (C) Copyright 1995-2001 Gary King
##  All Rights Reserved.
##*/
##/*
##    {b,vc}=quadcml(x,Zb,Zw,y);
##  Ecological Inference Likelihood Maximization via CML
##  With variance-covariance computed via global methods.
## 
##  MODEL:
##  mean function:  E(Y)=Bb*X + Bw*(1-x)
##  var  function:  V(Y)=vb*X^2 + c(b,w)*2*X*(1-X) + vw*(1-X)^2
##                  var params reparam'd: vars>0 and p.d.
##  distribution:   distribution implied on Y if Bb,Bw are bivariate 
##                  truncated normal
##  INPUT:
##  x = explanatory variable 
##  Zb = 1 or covariates for Bb (constant term included automatically)
##  Zw = 1 or covariates for Bw (constant term included automatically)
##  y = dependent variable
**
##  OUTPUT:
##  b = {Bb, Bw, sb, sw, rho}, where Bb, Bw are vectors if Zb,Zw have vars
##  vc = estimation global var cov matrix of b
##**
##  where params of interest are reparameterized as eirepar.g
##**
##  GLOBALS
##  _Erho[1]=standard deviation of normal prior on phi_5 (default=0.5)
##        0 fix it to _Erho[2] and don't estimate
##       <0 estimate without prior
##  _Ebounds=1 set bounds automatically unless z's are included
##           0 don't use bounds
##           kx2 or 1x2 matrix to indicate upper~lower bounds 
##  _Eprt=0 print nothing
##        1 print only final output from each stage
##        2 print everything useful plus friendly iteration numbers etc
##        3 print everything useful, iterations, and all sorts of checks
##  _Estval = 1 use default starting values or set to vector
##  _Eselect = scalar 1 to use all observations
##             vector of 1's to select and 0's to discard observations during
##             likelihood estimation.
*/
#include ei.ext;
quadcml <- function(x,Zb,Zw,y) {
  local vars,stval,ret,stb,t,bb,bw,sb,sbw,sw,vc,bnds,nbnds,tdiag,logl,
  b,tt1,tt2,tt,dataset,vrs,Dvrs,r,rho,pb,mlogl,grds,se,mask,fmt,eta,
  gridl,gr1;

  if _Eprt>=2;
    printfl "Likelihood estimation...";
  endif;

  /* prepare data */
  dataset = packdta(x,Zb,Zw,y);
  if _Eprt>=2;
    if _EselRnd/=1;
      fmtt;
      "_EselRnd: Modifying _Eselect with random selection = ";;_EselRnd;
    endif;
    tt=sumc(1-_Eselect);
    if (rows(_Eselect)/=1) and (tt/=0);
      fmtt;
      "_Eselect: deleting " tt " observations "\
      "from estimation stage;";
      "                   " sumc(_Eselect) "observations remain";
    endif;
  endif;
  
  /* starting values OR grid search */
  if rows(_Estval)==1 and not(scalone(_Estval));  @ grid search @
    if _Estval==0;
      gridl=5;      @ default number of gridlines for grid search @
    else;
      gridl=_Estval;
    endif;
  else;    
    gridl=0;        @ no grid search @
  endif;

  if scalone(_Estval) or gridl/=0;
    stval=zeros(sumc(_Ez),1)|-1.2|-1.2|0;
  else;
    stval=_Estval;
  endif;
  if rows(_Eeta)==4;
    stval=stval|_Eeta[1 2];
  elseif _Eeta[1]==4;
    stval=stval|0|_Eeta[2];
  elseif _Eeta[1]==5;
    stval=stval|_Eeta[2]|0;
  else;
    stval=stval|0|0;
  endif;
  
  /* cml globals */
  __title="EI Likelihood Maximization: "$+eiread(_Eres,"titl");
  _cml_Active=ones(rows(stval)-2,1)|0|0;
  _cml_MaxIters=_Emaxiter;
  _cml_DirTol=_EdirTol;
  _cml_CovPar=0;
  _cml_ParNames=(0+"Zb"$+ftosm(seqa(0,1,_Ez[1]),0,0))|
           (0+"Zw"$+ftosm(seqa(0,1,_Ez[2]),0,0))|(0+"sigB")|
	   (0+"sigW")|(0+" rho")|(0+"etaB")|(0+"etaW");

  @ if change this code, change also eiread @
  if scalzero(_Ebounds);    @ don't use bounds  @
    _cml_bounds={-1e256 1e256};
  elseif cols(_Ebounds)==2;
    _cml_bounds=_Ebounds|(0~.0001)|(0~.0001);
  elseif scalone(_Ebounds);     @ automatic bounds calculation @
    bnds={-10 10};
    nbnds={-20 20};
    if _Ez[1]==1;
      _cml_Bounds=bnds;
    else;
      _cml_Bounds=nbnds.*ones(_Ez[1],1);
    endif;
    if _Ez[2]==1;
      _cml_Bounds=_cml_Bounds|bnds;
    else;
      _cml_Bounds=_cml_Bounds|(nbnds.*ones(_Ez[2],1));
    endif;
    _cml_Bounds=_cml_Bounds|(-6~3)|(-6~3)|(-2~2);
  else;
    "quadcml: problem with _Ebounds";
    end;
  endif;

  if gridl==0 and cols(_Ebounds)/=2;  @ sneak past CML checks on parameters @
    tt=rows(stval);
    tt=tt-1|tt;
    _cml_bounds=_cml_bounds|(stval[tt]~stval[tt]);
    tt=rows(_cml_bounds);
    tt=tt-1|tt;
    _cml_bounds[tt,2]=_cml_bounds[tt,2]+_Edirtol;
  endif;
  
  if _Eprt==0;
    __output=0;
  elseif _Eprt<=2;
    #ifdos;
      _cml_Diagnostic={.};
      __output=2;
    #else;
      _cml_Diagnostic=2;
      __output=1;
    #endif;
  elseif _Eprt==3;
    _cml_Diagnostic=3;
    __output=2;
  else;
    "quadcml: problem with _Eprt";
    end;
  endif;
  
  if _Erho[1]==0;
    r=rows(stval)-2;
    stval[r]=_Erho[2];
    _cml_bounds[r,.]=(_Erho[2]-__macheps)~(_Erho[2]+__macheps);
    _cml_active[r]=0;
  endif;

  if gridl==0;         /* run CML */
    {b,mlogl,grds,vc,ret}=cml(dataset,0,&eiloglik,stval);
    if in(ret,3|4|6|7|10|99,1);
      ?;
      "*********************************************";
      "CML return code: ";;ret;
      "Restarting iterations with trust algorithm on";
      "*********************************************";
      _cml_options={ trust };
      {b,mlogl,grds,vc,ret}=cml(dataset,0,&eiloglik,stval);
      _cml_options=0;
    endif;
    _Eres=vput(_Eres,ret,"retcode");
    logl=mlogl*_cml_NumObs;
    _Eres=vput(_Eres,logl,"loglik");
    if _Eprt>=2;
      printfl "CML converged; Computing variance-covariance matrix...";
      ?;
    endif;
  
  else;                /* run GRID search */
    if _Eprt>=2;
      "Preliminary mean grid search (on 2 parameters)...";
    endif;
    tt=rows(_cml_bounds);
    gr1=meanc(_cml_bounds[3:tt,.]');
    gr1=_cml_bounds[1 2,.]|(gr1~gr1);
    {b,logl}=eigrid(dataset,gr1,53,_Edirtol);
    _cml_bounds[1 2,.]=b[1 2]+((-1~1)|(-1~1));
    if _Eprt>=2;
      ?;"Main grid search (on all parameters)...";
    endif;
    {b,logl}=eigrid(dataset,_cml_bounds,gridl,_Edirtol);
    ret=3333;
    _Eres=vput(_Eres,ret,"retcode");
    _Eres=vput(_Eres,logl,"loglik");
    grds=b*0+_Edirtol;
  endif;
  
  /* compute var cov matrix */
  if _EdoSim/=-1;
    clearg _GhFix;
    _GhFix=ones(rows(stval)-2,1)|0|0;
    if _Erho[1]==0;
      _GhFix[r]=0;
    endif;
    vc=gvc(&eiloglik,b,dataset);
    if _Erho[1]==0;
      vc=vc~zeros(rows(vc),1);
      vc=vc|zeros(1,cols(vc));
      vc[rows(vc),cols(vc)]=0.00000000001;
    endif;
    _Eres=vput(_Eres,_GhActual,"GhActual");
    if scalmiss(vc);
      "quadcml: couldn't compute positive definite variance matrix";
      retp(miss(1,1),miss(1,1));
    endif;
  else;
    vc=zeros(rows(b)-2,rows(b)-2);
  endif;
  
  /* print likelihood results */
  if _Eprt>=1 or ret/=0;
    if ret/=0 and ret/=3333;
      "Overriding no print instructions.  See CML return code";
    endif;
    if ret==20;
      ret=1;
    endif;
    format/ro 7,4;
    ?;
    __title;
    "CML Return code:          ";;ret;
    "log-likelihood:           ";;logl;
    "Number of observations:   ";;_cml_NumObs;
    "Number of iterations:     ";;_cml_IterData[1];
    if _ei_vc[_Ghactual,1]/=-1;
      "Variance computed from:   ";;_ei_vc[_Ghactual,.];
    else;
      "-Hessian computed from:   ";;_ei_vc[_Ghactual,.];
    endif;
    ?;
    "*## Parameter values in scale of estimation ***";
    mask=0~1~1;
    fmt=("*.*s"~  1~8)|
        ("*.*lf"~ 10~4)|
	("*.*lf"~ 10~4);
    if _ei_vc[_Ghactual,.]/=-1;
      se=sqrt(diag(vc))|0|0;
    else;
      se=miss(0,0).*ones(rows(b),1);
    endif;	
    call printfm(_cml_ParNames~b~se,mask,fmt);
    ?;
    {bb,bw,sb,sw,rho}=eirepar(b,Zb,Zw,x);
    let vrs=bb bw sb sw rho;
    "Reparameterization back to the truncated scale, parameterized";
    "        according to the underlying untruncated distribution.";
    $vrs';
    pb=meanc(bb)|meanc(bw)|sb|sw|rho;
    pb';
    if _Eprt>=3;
      ?;
      "Reparameterization back to ultimate truncated scale";
      pb=eirepart(b,Zb,Zw,x);
      $vrs';
      pb';
    endif;
  endif;  

  _Eres=vput(_Eres,_cml_parnames,"parnames");
  _Eres=vput(_Eres,b,"phi");
  _Eres=vput(_Eres,vc,"vcphi");
  
  retp(b,vc);
endp;

