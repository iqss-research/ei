/*
**
   compare.prg 
**
** This program shows how to compare different redistricinting
** plans using the procedure, compelct.g
**
** distplan = redistricting plan indicator variable.
**            1=PLAN1001, 2=PLAN1021, 3=PLAN1025, 4=PLAN1034, 
**            5=PLAN1040, 6=PLAN1043, 7=PLAN1044, 8=PLAN1045, 
**            9=PLAN1046, 10=PLAN1047, 11=PLAN1048, 
**            12=PLAN1000 
*/

dists={10,10,10,32,10,10,31,10,10,10,10,10};
compelct(dists,"res10b.out");
dists={12,24,31,24,24,24,32,5,24,24,24,24};
compelct(dists,"res24b.out");
dists={25,25,25,9,25,25,25,25,25,25,25,25};
compelct(dists,"res25b.out");
