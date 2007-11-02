###DESCRIPTION : trims leading and trailing blanks for any string but not
###              inter-words blanks.  For example, trim.blanks("   abc  ") = "abc";
###              but trim.blanks("   abc    def ")= "abc    def".
### AUTHOR: Elena Villalon
##          evillalon@iq.harvard.edu
###
trim.blanks <- function(x) {
### at the beginning of string"^" gets anny number (+) of white spaces
  f <- x
  if(length(x))
    f <- na.omit(x)
  
  if(length(f) <= 0)
    return(x)
  if(length(f)>1)
    print(f)
  if(f=="" )
    return(x)
  x <- sub("^[[:space:]]*(.*)", "\\1",x) ###get \n\t
  x <- sub('^ +', '', x) ###get white spaces
  
### at the ending of string"$" gets anny number (+) of white spaces
  
  x <- sub("(.*)[[:space:]]*$", "\\1", x)
  x <- sub(' +$', '', x)
  return(x)
}
###Request input from user 
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
### Used instead of message to control verbose output
messout <- function(str, verbose=T, obj=NULL){
  if(verbose)
  message(str);
  if(length(obj) >0 && verbose)
    print(obj)
}
 

