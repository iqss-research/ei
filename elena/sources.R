vec <- c("eibounds.R","eicml.R","eidata.R", "eienv.R","eigrid.R","eiloglik.R","einonpar.R","einormal.R",  
"eiproc.R", "eiread.R", "eirepar.R", "eitnorm.R","eiutils.R","g2rutils.R", "probs.R","eisims.R", 
"testcases.R","utils.R")
print(getwd())
message("You need to be in directory ~/GAUSSCODE/ei/R")
lapply(vec,source)
### cml.Rdata contains stval, cml.bounds, dataset
### quacml.Rdata contains b, vc from R and bg, vcg from Gauss
getSample()
###dbuf <- ei(t0,x0,tvap,1,1,dbug=TRUE,EnonPar=1)
###dbuf$betabs
### qq3 <- quadcml(x0,1,1,t0,evglobal)
### lst <- eirepar(qq3$b,1,1,x0,get("Ez", env=evglobal),evglobal)
### pb <- eirepart(qq3$b,1,1,x0,get("Ez", evglobal)
### load("cml.Rdata")
### llk <- eiloglik(b, dataset,evbase=evglobal)
