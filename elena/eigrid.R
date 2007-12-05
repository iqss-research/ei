##/*  
##**  This archive is part of the program EI
##**  (C) Copyright 1995-2001 Gary King
##**  All Rights Reserved.
##**
##**  Grid search procedures
##**  by Eric Dickson & Gary King
##**
##**  {b,logl}=eigrid(dataset,gridlohi,gridline,toler);
##**
##**  INPUTS:
##**
##**  dataset = input to loglik proc defined as required by CML
##**        (loglik must be called eiloglik(b,dataset) and defined separately)
##**
##**  gridlohi = k x 2 vector.  Columns designate lower~upper bounds of the
##**           search region scale of estimation for each of the parameters.
##**           In EI, _cml_bounds, as constructed by _Ebounds, is used.
##**           To hold a parameter fixed, set both entries equal to each other.
##**
##**  gridline = scalar; number of gridlines to compute per dimension.
##**             Is a constant for all dimensions over which variation
##**             is desired.   So far only 3 and 5 have been tested.
##**
##**  toler = scalar, specifies the accuracy to which the 
##**          optimal gridpoint must be calculated.  Once we have
##**          zoomed in sufficiently, and the gridlines are within
##**          toler of one another, the routine ceases.
##**
##**  OUTPUTS:
##**
##**  b = the gridpoint at which the log-likelihood function is maximized.
##**      It is a kx1 vector.
##**
##**  logl = the value of the log-likelihood function at b.
##*/

#include ei.ext;
eigrid <- function(dataset,gridlohi,gridline,toler,evbase=get("evbase", parent.frame()))
  {
    dta <<- dataset
    logrid <<- gridlohi
    hggrid <<- gridline
    
    gridlo <- gridlohi[,1];
    gridhi <- gridlohi[,2];
    str <- "********************************************************************";

    gridvar <- gridhi - gridlo;

    ##/* number of zoom-ins to get to _Edirtol level of accuracy */
    k <- floor((log(toler/max(gridvar)))/(log(2/(gridline + 1))))+ 1;

    grid <- t(makegrid(gridlo,gridhi,gridline)) ### /* get initial grid */

###/* calculate the likelihood of each grid point */
    for (j in (1:k)){ ###               /* zoom-in number */
      jj <- j
      print(str)
###  format/ldn 4,0;
      print(paste("Beginning zoom-in ", jj))
      print(paste(" of ", k))
      rr <- cols(grid)
      loglk <- matrix(0,nrow=rr,ncol=1)
      et1 <- hsec <- proc.time()[3]
      for (f in 1:rr) { ###          /* regulates grid-line to be analyzed */
        ff <- f
        if ((jj==1) && (ff==200)){
          et2 <- hsec <- proc.time()[3]
          et3 <- (et2 - et1)*rr*k/ff
###  format/ldn 1,1;
          print("The grid search will end approximately")
          print(paste((et3-(et2-et1))/6000," minutes from now (", proc.time()[3], ")"))
        }
        loglk[ff] <- colSums(as.matrix(eiloglik(grid[,ff],dataset,evbase)))
    
      }

###/* Find grid point with highest likelihood and zoom in */
      best <- maxindc(loglk) ###        /* max like row number */
### format/ldn 7,2;
      print(paste("Best grid point, loglik: ", loglk[best]))
      gridopt <- grid[,best]
### format/ldn 6,3;
      print(paste("Coordinates:             ", as.vector(gridopt)))
      print(paste("Scale explored:          ", as.vector(gridvar)))
      if (jj <= k-1){
        gridvar <- gridvar * 2 / (gridline + 1)
        gridhi <- gridopt + gridvar/2
        gridlo <- gridopt - gridvar/2
        diffgrid <- (gridvar <toler)
        for (m in (1:rows(diffgrid))){
          mm <- m
          if (scalone(diffgrid[mm])){  
            gridlo[mm] <- (gridlo[mm]+gridhi[mm])/2
            gridhi[mm] <- gridlo[mm]
          }
        }
        grid <- makegrid(gridlo,gridhi,gridline)
        grid <- as.vector(grid)
      }
    }
    print(str)
    lst <- c(list(gridopt=gridopt),list(loglk=loglk[best]))
    return(lst)
  }

##/**************************************************************************
##**   makegrid.src
##**   version of 9 June 1998
##**   by Eric S. Dickson
##**
##**   This proc, given the boundaries and level of detail, produces grids 
##**   which are used in EIGRID()
##**
##**
##**          {grid} = makegrid(gridlo,gridhi,gridline);
##**
##**
##**   INPUTS:
##**
##**   gridlo = _Egridlo, a k x 1 vector analogous to _Estval.  Each entry
##**           designates the lower bound of the search region in the
##**           scale of estimation for one of the parameters.
##**
##**   gridhi = _Egridhi, a k x 1 vector analogous to _Estval.  Each entry
##**           designates the upper bound of the search region in the
##**           scale of estimation for one of the parameters.
##**
##**     NOTE ON USAGE: If one wishes to hold a particular parameter fixed
##**      at a certain value, then one should set gridlo = gridhi for
##**      that parameter.
##**
##**   gridline = _Egrline, the scalar number of gridlines per dimension.
##**             Is a constant for all dimensions over which variation
##**             is desired.   So far only 3 and 5 have been tested.
##**
##**   OUTPUTS:
##**
##**   grid = a matrix whose rows are the various gridpoints.  The order
##**          of the columns is the same as the order of the rows in
##**          gridlo or gridhi, and cols(grid) = rows(gridlo). 
##**          rows(grid) = gridline^(gridrows), where gridrows is the 
##**          number of parameters which are varied in the grid search.
##**
##**   In the following example, gridline = 7, and various variables of
##**   interest are labelled:
##**
##**   gridlo                                                         gridhi 
##**        }---->|     *     *     *     *     *     *     *     |<----{   
##**                    ^                                   ^  
##**                    |                                   |
##**                truedown                             trueup
##**
##**
##*************************************************************************/

makegrid <- function(gridlo,gridhi,gridline, evbase=get("evbase", parent.frame()))
  {


    glo <- gridlo;##      /* vector containing lower ranges in scale of estimn */
    ghi <- gridhi;##      /* vector containing upper ranges in scale of estimn */
    j <- gridline;##    /* the number of gridlines per dimension */
    l <- rows(ghi);##   /* the number of dimensions */

    griddiff <- ghi - glo; ##                /* difference vector: the range */
    griddum <- (griddiff !=  0); ###         /* 1 if ghi /= glo, else 0 */
    gridsum <- colSums(as.matrix(griddum)); ###            /* how many dim's will be varied */
    gridnums <- seq(from=1, by=1,length.out=rows(griddum))
###/* convenient labelling device */
    trueup <- glo+(1/(1+j))*griddiff;  ### /* vector with highest gridpoint coords */
    truedown <- ghi-(1/(1+j))*griddiff; ###  /* vector with lowest gridpoint coords */

###/* The way the third (main) if-loop is written at present, it cannot handle 
###   gridum = 0 or gridsum = l because of the use of trimr.  I am sure there is
###   a more elegant workaround than adding the two additional loops, but these
###   work okay for right now.  There is probably also a better way of doing the
###   main loop in general, but this is the best I've tried up to now.  */

    if (gridsum == l){ ###       /* i.e., if we are to vary over ALL dimensions */
      gridpts <- matrix(0,nrow=j,ncol=l)         
      for (k in(1:l)){
        kk <- k
        
        gridpts[,kk] <- seq(from=truedown[kk],to=trueup[kk],length.out=j)
      }
      
      grid <- makefacn(l,gridpts); 
    }

    if (gridsum == 0)###  /* to vary over NO dimensions: i.e., ghi = glo = grid' */
      grid <- t(ghi);
  
    if (gridsum !=0 && gridsum!=l)
      {### /* the complicated case where only some vary */
      
        gridpack <- cbind(griddum,gridnums, glo, ghi, truedown, trueup)
        gridpack <- sortc(gridpack,1)
        gridcut <- gridpack[1:(rows(gridpack)-gridsum),2:4]
        gridpack <- trimr(gridpack,rows(gridpack)-gridsum,0)
        gridpack <- sortc(gridpack,c=2)
        gridrows <- rows(gridpack)
        gridpts <- matrix(0,nrow=j,ncol=gridrows)
        grid <- matrix(0, nrow=j^(gridrows),ncol=l)
  
        for( i in (1:gridrows)){
          ii <- i
     
          gridpts[,ii] <- seq(from=gridpack[ii,5],to= gridpack[ii,6],length.out=j)
        }

        gridint <- makefacn(gridrows,gridpts)

        for (m in (1:rows(gridcut))){
          mm <- m
          grid[,gridcut[mm,1]] <- gridcut[mm,2] * matrix(1,nrow=rows(grid),ncol=1)
        }

        for (n in (1:gridrows)){
          nn <- n
          grid[,gridpack[nn,2]] <- gridint[,nn]
        }
      }
   
    return(grid)
  }
