  
/*
**  usage:    outz = cln(inz);
**  This procedure discards precincts with unusable data
**   prior to EI estimation.  Note, the precinct will only be
**   thrown out for this election and year and may be used in other
**   analyses if data are present.  (for example, some precincts have 
**   election results for some years but not others).
**  Problems that can create such precincts are:
**     Missing or incorrect census information.
**     Missing election returns.
**     Election returns that add up to a higher number than
**       people listed as living in the precinct.
**  Precincts are thrown out if: 
**     Turnout is greater than 100% or is 0% or less than 0%.
**     Number of people in the precinct is 0%.
**     Proportion of ethnic group is higher than 100% or less than 0%.
**     Proportion of democratic vote is higher than 100% or less than 0%.
**  Inputs: temp1
**  Outputs: temp1
*/

proc cln(temp1);
  temp1 = packr(temp1); 
  temp1 = delif(temp1,(temp1[.,1].>1).or(temp1[.,1].<=0).or(temp1[.,2].<= 0) 
                   .or(temp1[.,3].>1).or(temp1[.,3].< 0).or(temp1[.,4].>  1) 
                   .or (temp1[.,4] .< 0));
  retp(temp1);
endp;  



