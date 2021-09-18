/*
**  This proc sets 4 globals for individual runs of EI 
**   and EI2.  See the EI manual for explanation of 
**   how each global functions.
**  
*/

proc 0=globals(eeta,ealphab,ealphaw,erho);
   _Eeta = eeta;
   _EalphaB = ealphab; 
   _EalphaW = ealphaw;
   _Erho = erho;
endp;