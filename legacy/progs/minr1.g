/*   This proc. chooses the lowest value in a row.
**   The data set MASTER was created by concatenating
**   three data sets for 3 seperate redistricting plans.
**   In the very rare case where a precinct was split
**   differently in two plans, the population figures for 
**   each plan in the precinct double counted some individuals.
**   Using this proc ensured that no double counting occurs in 
**   our analysis.
*/
 
proc 1=minr1(a);
  local res;
  res=minc(a');
  retp(res);
endp;
