###
##  here are some function that exists in Gauss itself but nor in R
##



###
##    m <- zeros(r,c)
##    r - nr of rows
##    c - nr of columns
##    m - matrix of zeros
##
zeros<-function(r,c){
        if(c==1)
          res<-rep(0,r)
        else
          res<-matrix(0,r,c)
        return(res)
}


###
##    m <- ones(r,c)
##    r - nr of rows
##    c - nr of columns
##    m - matrix of ones
## 

ones<-function(r,c){
        if(c==1)
            res<-rep(1,r)
        else
            res<-matrix(1,r,c)
        return(res)
}

###
##    res <- menac(x)
##    x   - matrix
##    res - vector with mean of each col. in matrix x 
## 
meanc<-function(x){
        return(colMeans(x))
}


###
## counts the number of elemetns of a vector that fall into a specific rage
## c = count(x,v)
## x = Nx1 vector containing the numbers to be counting
## v = Px1 vector (sorted in asc order) containing the ranges within which counts are to be made

counts<-function(x,v){
        res<-c()
        v<-sort(v)
        for (i in 1:length(v)){
                tmp<-x[(x<=v[i])]
                res<-c(res,length(tmp))
                x<-setdiff(x,tmp)
        }
        return (res)
}


## remove t first and b last rows from x

trimr<-function(x,t,b){
        return(x[(t+1):(nrow(x)-b),])
}

###
## removes the last row of the original matrix AND add a
## row with "missing values" (as the first row of the original matrix)

lag<-function(x){
        if(is.matrix(x))
          res<-rbind(rep(NA,ncol(x)),x[1:(nrow(x)-1),])
        else
          res<-c(NA,x[1:(length(x)-1)])
        return (res)
}
