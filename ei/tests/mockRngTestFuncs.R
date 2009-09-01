#####File Paths#####
#print(Sys.getenv("RCMDCHECK", unset="unset"))
#basepath <- ifelse(Sys.getenv("RCMDCHECK", unset="unset") == "FALSE", ".", file.path("..","ei","unitTests"))

#####mock RNGs#####
mock_rnorm <- function(n,mean=0,sd=1){
    return(sapply(1:n, function(x) (KR_rnorm()+mean)*sd))
}

KR_rnorm <- function(){
    xi = 2.216035867166471
    u <- mock_runif(1)
    if(u < 0.884070402298758){
        v <- mock_runif(1)
        return(xi*(1.131131635444180*u+v-1))
    }
    else if(u < 0.973310954173898){
        if(u < 0.958720824790463){
            if(u < 0.911312780288703){
                while(1){
                    vw <- mock_runif(2)
                    z <- vw[2]-vw[1]
                    t <- 0.479727404222441 - 0.595507138015940 * min(vw)
                    if(max(vw) <= 0.805577924423817 || 0.053377549506886*abs(z) <= KR_f(t))
                        return(ifelse(z<0, t, -t))
                }
            }
            else{
                while(1){
                    vw <- mock_runif(2)
                    z <- vw[2]-vw[1]
                    t <- 0.479727404222441 + 1.105473661022070 * min(vw)
                    if(max(vw) <= 0.872834976671790 || 0.049264496373128*abs(z) <= KR_f(t))
                        return(ifelse(z<0, t, -t))
                }
            }
        }
        else{
            while(1){
                vw <- mock_runif(2)
                z <- vw[2]-vw[1]
                t <- xi - 0.630834801921960*min(vw)
                if(max(vw) < 0.755591531667601 || 0.034240503750111*abs(z) <= KR_f(t))
                    return(ifelse(z<0, t, -t))
            }
        }
    }
    else{
        vw <- c(10,10)
        t <- 10
        while(vw[2]^2*t > xi^2/2){
            vw <- mock_runif(2)   
            t <- xi^2/2 - log(vw[1])
        }
        return(ifelse(u < 0.986655477086949, sqrt(2*t), -sqrt(2*t)))
    }
}

KR_f <- function(t){
    xi <- 2.216035867166471
    dnorm(t)-0.180025191068563*(xi-abs(t))
}

mock_runif <- function(n){
    old_cursor <<- cursor+1
    cursor <<- cursor+n

    return(rndu_list[old_cursor:cursor])
}

reset_mock <- function(){
    cursor <<- 0
}

#####Setup Calls#####
.setUp <- function(){
    require(ei)
    .mockRNGState <<- new.env(parent=.GlobalEnv)
    .mockRNGState$cursor <<- 0
    load("rndu_1.RData", envir=.mockRNGState)
    #load("rndu_1.RData", envir=.mockRNGState) 
    environment(mock_rnorm) <<- .mockRNGState
    environment(mock_runif) <<- .mockRNGState
    environment(reset_mock) <<- .mockRNGState
}
######Global parameters#####
qis <- c('aggs',
        'beta',
        'betaBs',
        'betaWs',
        'CI50b',
        'CI50w',
        'CI80bw',
        'CI95bw',
        'CsbetaB',
        'CsbetaW',
        'EaggBias',
        'ExpVarCI',
        'ExpVarCI0',
#        'ExpVarCIs',
#        'GEbw',
#        'GEbwa',
#        'GEwb',
#        'GEwba',
        'Maggs',
        #'mpPsiu',
        'Paggs',
        'psi',
        'R',
        'sbetaB',
        'sbetaW',
        'tsims',
        'tsims0',
        'VCaggs',
        'ABounds',
        'bounds',
        'phi',
        #'psiu',
        'VCphi')

good_datasets <- c("CENS1910", "KYCK88", "MATPROII", "NJ")
bad_datasets <- c("LAVOTEALL", "SCSP")
ugly_datasets <- c("FULTONGEN")

datasets <- c(good_datasets)

#####Actual test functions#####
test.cens1910 <- function(){
    print(mock_rnorm)
    run_all("CENS1910")
}

test.kyck88 <- function(){
    run_all("KYCK88")
}

test.matproII <- function(){
    run_all("MATPROII")
}

test.nj <- function(){
    run_all("NJ")
}

run_all <- function(dataset){
#    format_gauss_all(dataset, qis)
    format_r_all(dataset, qis)
    qi_test(dataset)
}

qi_test <- function(dataset){
    load(file.path(dataset, "r_res.RData"))
    #load(sprintf("%s/%s/r_res.RData", basepath, dataset))
    load(file.path(dataset, "g_res.RData"))
    #load(sprintf("%s/%s/g_res.RData", basepath, dataset))

    res <- list()

    for(q in qis){
        print(q)
        x <- all.equal(g_res[[q]], round(r_res[[q]],4), tol=1e-3, check.attributes=F)
    }
}    

#####GAUSS formatting functions#####

format_gauss_all <- function(datasets, qis){
    for(d in datasets){
        print(d)
        g_res <- format_gauss(d, qis)
        save(g_res, file=file.path(d,"g_res.RData"))
        #save(g_res, file=paste(basepath, '/', d, "/g_res.RData", sep=''))
    }
}

format_gauss <- function(dataset, qis){
    g_res <- list()
    for(q in qis){
        print(q)
        table <- as.matrix(read.table(file.path(dataset, sprintf("test%s.tab", q))))
        #table <- as.matrix(read.table(paste(basepath, '/', dataset, "/test", q, ".tab", sep='')))
        print(dim(table))
        d <- qidims(dataset,q)
        if(q=="VCphi") print(table)
        g_res[[q]] <- matrix(as.numeric(table), byrow=T, nr=d[1], nc=d[2])
    }

    return(g_res)
}

#####R formatting functions#####
format_r_all <- function(datasets, qis, ...){
    for(d in datasets){
        print(d)
        r_res <- format_r(d, qis, ...)
        print(.mockRNGState$cursor)
        save(r_res, file=file.path(d,"r_res.RData"))
        #save(r_res, file=paste(basepath, '/', d, "/r_res.RData", sep=''))
    }
}

format_r <- function(dataset, qis, ...){
    r_res <- list()
    #.setUp()
    dat <- read.table(file.path(dataset, sprintf("%s.tab", dataset)))
    #dat <- read.table(paste(dataset,'/',dataset,".tab",sep=''))
    t <- dat[,1]
    x <- dat[,2]
    n <- dat[,3]
    if(dim(dat)[2] > 3){
        Zb <- dat[,4]
        Zw <- dat[,5]
    }
    else{
        Zb <- 1
        Zw <- 1
    }
    
    load(file.path(dataset, "g_res.RData"))
    #load(paste(dataset, "/g_res.RData", sep=''))
    vcphi <- g_res[["VCphi"]]
    phi <- g_res[["phi"]]
    rm(g_res)
    dbuf <- ei(t,x,n,Zb,Zw,EdoML=0,EdoML.phi=phi,EdoML.vcphi=vcphi,cml.bounds=NA, dbug=T)
    #dbuf <- ei(t,x,n,Zb,Zw, ...)
    for(q in qis){
        reset_mock()
        print(q)
        table <- as.matrix(eval(parse(text=sprintf("summary(dbuf, %s)", q))))
        d <- qidims(dataset,q)
        r_res[[q]] <- matrix(table, nr=d[1], nc=d[2])
    }

    return(r_res)
}

#####Formatting parameters#####
get_p <- function(dataset="sample"){
    switch(dataset,
        CENS1910=1040,
        FULTONGEN=289,
        KYCK88=118,
        LAVOTEALL=3262,
        MATPROII=275,
        NJ=567,
        SCSP=3187,
        75)
}

qidims <- function(dataset,qi){
    p<-get_p(dataset)
    Esims<-100
    switch(qi,
        aggs=c(Esims,2),
        beta=c(p,2),
        betaBs=c(p,Esims),
        betaWs=c(p,Esims),
        CI50b=c(p,2),
        CI50w=c(p,2),
        CI80bw=c(p,4),
        CI95bw=c(p,4),
        CsbetaB=c(p,1),
        CsbetaW=c(p,1),
        EaggBias=c(4,2),
        ExpVarCI=c(100,4),
        ExpVarCI0=c(p,4),
        ExpVarCIs=c(100,4),
        GEbw=c(p,3),
        GEbwa=c(p,2),
        GEwb=c(p,3),
        GEwba=c(p,2),
        Maggs=c(2,1),
        mpPsiu=c(5,1),
        Paggs=c(2,2),
        psi=c(5,1),
        R=c(1,1),
        sbetaB=c(p,1),
        sbetaW=c(p,1),
        tsims=c(100,Esims+1),
        tsims0=c(p,Esims+1),
        VCaggs=c(2,2),

        ABounds=c(2,2),
        bounds=c(p,4),
        phi=c(7,1),
        psiu=c(5,1),
        VCphi=c(5,5))
}
