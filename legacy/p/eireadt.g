/* quick way to run eiread(dbuf,"sum") and save output in file called t
*/
proc 0=eireadt;
  eiset;
  format/rd 7,4;
  output file=t reset;
     call eiread(dbuf,"sum");
     ?;
     call eiread(dbuf,"neighbor");
     ?;
     call eiread(dbuf,"goodman");
     output off;
endp;

