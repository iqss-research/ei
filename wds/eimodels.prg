/*
** A sample code from the manual (_Esims = 100 so that it runs quickly...)
*/

new;                                     @ clear workspace @
library ei;                              @ initialize libraries @
clear t,x,n;                             @ clear all variables in dataset @
loadvars sample.asc t x n;               @ load variables from disk file @

eiset;                                   @ clear for Model 1 @
_Eres = vput(_Eres, "Model 1", "titl");  @ print out model number @
_Eeta = 1;                               @ zb=x, zw=1 @
_Ealpha_B = {0 2};                       @ prior for betaB @
_Esims = 100;                           @ number of simulation @
_EiLlikS = 1;                            @ store log-likelihood at each simulation @
dbufdef = eimodels_def("",1,t,x,n,1,1);  @ save Model 1, the first model @

eiset;                                   @ clear for Model 2 @
_Eres = vput(_Eres, "Model 2", "titl");
_Eeta = 2;                               @ zb=1, zw=x @
_Ealpha_W = {0 2};                       @ prior for betaW @
_Esims = 100;
_EiLlikS = 1;
dbufdef = eimodels_def(dbufdef,2,t,x,n,1,1); @ save Model 2 @

eiset;                                   @ clear for Model 3 @
_Eres = vput(_Eres, "Model 3", "titl");
_Eeta = 3;                               @ zb=x, zw=x @
_Ealpha_B = {0 2};                       @ prior for betaB @
_Ealpha_W = {0 2};                       @ prior for betaW @
_Esims = 100;
_EiLlikS = 1;
dbufdef = eimodels_def(dbufdef,3,t,x,n,1,1); @ save Model 3 @

save rvdef = dbufdef;                    @ save ndbuf in file rvdef.fmt @
dbufrun = eimodels_run(dbufdef);         @ run ei on all the models @
save rvrun = dbufrun;                    @ save ndbres in file rvrun.fmt @

_EIMetaR = 1;
call eiread(dbufrun, "sum");             @ summary of ei run for Model 1 @

_EIMetaR = 3;
graphon;
call eigraph(dbufrun, "tomog");          @ tomography plot for Model 3 @

_EI_bma_prior = {0.2,0.4,0.4};           @ set prior model probabilities @
dbufavg = eimodels_avg(dbufrun);         @ Bayesian Model Averaging @
save rvavg = dbufavg;                    @ save dbres in file rvavg.fmt @

call eiread(dbufavg,"sum");              @ summary for dbres @
graphclr;                                @ clear graphics screen @
call eigraph(dbufavg,"post");              @ draw district level posterior @
graphoff;                                @ close graphics window @


