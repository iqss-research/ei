/* 
** This proc interpolates population data for years between censuses.
** Census data are available for 1990 and 2000.
** INPUTS:
**   var1: the value of the variable at the beginning of the period.
**   var2: the value of the variable at the end of the period.
**   num:  The number of values to be interpolated.
** OUTPUT:
**   res: a matrix of original and interpolated values.
**
*/

proc (1)=intpl(var1,var2,num);
  local res;
  res=var1;
  for i(1,num,1);
    res=res~(var1+(i/(num+1)).*(var2-var1));
  endfor;
  res=res~var2;
  retp(res);
endp;