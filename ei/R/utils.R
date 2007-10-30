user.prompt <- function (verbose=TRUE){
  verb <- T
 if(length(verbose) <= 0 || verbose){
   silent <- readline("\nPress <return> to continue or Ctrl-c Ctrl-c to quit: ")
   
 }else{
   
   answer <- readline("\nPress 'Y' for verbose output or enter 'N' otherwise. \nPress Ctrl-c Ctrl-c to quit: " )
   
   if(substr(answer,1,1)=='N')
     verb <- F
   else if(substr(answer,1,1)=='Y')
     verb <- T
   
 }
   return(verb)
     
}
messout <- function(str, verbose=T, obj=NULL){
  if(verbose)
  message(str);
  if(length(obj) >0 && verbose)
    print(obj)
}
 

