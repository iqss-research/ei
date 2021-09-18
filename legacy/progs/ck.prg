/*
**
** This program provides an easy way to examine the two
** graphs: eigraphs(dbuf,"fit") and eigraph(dbuf,"tomogs").
**
*/

new;
library /home/kimai/ei/p/ei;

keyword kd(db);
  clearg dbuf1,dbuf2,dbuf3,dbuf4,ndbuf;
  clearg s;
  s="/home/kimai/texas/ac_98p/db"$+db;
  loadm ndbuf=^s;
  dbuf1=vread(ndbuf,"dbuf1");
  dbuf2=vread(ndbuf,"dbuf2"); 
  dbuf3=vread(ndbuf,"dbuf3");
  dbuf4=vread(ndbuf,"dbuf4"); 
  graphon;
  graphno;
endp;
  
keyword kf(n);
  local s;
  _eigraph_bwhi=0.4;
  s="dbuf"$+n;
  graphclr;
  eigraph(varget(s),"fit");
  graphno;
endp;

keyword kt(n);
  local s;
  s="dbuf"$+n;
  graphclr;
  _eigraph_bwhi=0.4;
  eigraph(varget(s),"tomogs");
  graphno;
endp;
