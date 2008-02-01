\name{ei}
\alias{ei}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Time-series cross-sectional Forecasting}
\description{
  Uses aggregate data to infer individual levels relations It uses
  non-parametric and bayesian model averaging procedures}
\usage{ei(t,x,n,Zb,Zw,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{Vector. Dependent variable, e.g. fraction of people who votes.}

  \item{x}{Vector.Explanatory variable, e.g, fraction of voting age people who are black.}
  
  \item{n}{Vector.Explanatory variable, e.g. number of voting people.}
  
  \item{Zb}{Optional covariate, which may be set equal to 1 for no
    covariate. Default 1.}
  
  \item{Zb}{Optional covariate, which may be set equal to 1 for no covariate. Default 1.}
  \item{"..."} {Additional inputs.}
}
\details{
  ~~ If necessary, more details than the description above ~~}
\value{
  Returns a list containing the following components: 

  \item{"eiversion"}{The curent version of the package}
  \item{"t"}{Vector with the dependent variable}
  \item{"x"}{Vector of explanatory variable}
  \item{"n"}{Vector of explanatory variable}             
  \item{"Zb"}{Optional Covariate or 1}
  \item{"Zw"}{Optional Covariate or 1}
  \item{"EalphaB"}{{\eqn{(ncol(Zb) \times 2)}} matrix of means (in the
    first column) and standard deviations (in the second) of an
    independent normal prior distribution on elements of
    \eqn{\alpha^b}. If you specify Zb, you should probably specify a
    prior, at least with mean zero and some variance (default=NA; which
    indicates no prior). (See Equation 9.2, page 170, to interpret
    \eqn{\alpha^b$}).}
  
\item{"EalphaW"}{\eqn{(ncols(Zw) \times 2)} matrix of means (in the
  first column) and standard deviations (in the second) of an
  independent normal prior distribution on elements of
  \eqn{\alpha^w}. If you specify Zw, you should probably specify a
  prior, at least with mean zero and some variance (default=NA; which
  indicates no prior). (See Equation 9.2, page 170, to interpret
  \eqn{\alpha^w$})}

\item{"Ebeta"}{Standard deviation of the ``flat normal'' prior on 
  \eqn{\breve{\mathfrak{B}}^b} and \eqn{\breve{\mathfrak{B}}^w}. The flat normal
  prior is uniform within the unit square and dropping outside the
  square according to the normal distribution. Set to zero for no prior
  (default). Setting to positive values probabilistically keeps the
  estimated mode within the unit square. 0.25 is a reasonable value to
  experiment with at first.}

  \item{"Ebounds"}{If equal to 1 sets CML bounds on parameters
    automatically unless z's are included; if 0 don't use bounds; 
    \eqn{k\times 2} (where \eqn{k} is the number of starting values) or 
    \eqn{1\times 2} matrix to indicate upper \eqn{\sim} lower bounds. (Do not
    confuse the bounds referred to here with the bounds on the
    quantities of interest.) Default=1.}
  
  \item{"Ecdfbvn"}{Determines which procedure to use for computing the
    area of the bivariate normal distribution above the unit square:
    \item{1}{based on the Gauss function CDFBVN;}
    \item{2}{Martin van der Ende's method (based on D.R. Divgi,
	 ``Calculation of the univariate and bivariate
    Normal integral,'' Annals of Statistics, 1979, 903-910, with
    additional options available for this method in the proc
    cdfbvn_div);}
  \item{3} {Integration of log of the unit square;}
  \item{4}{ Direct integration on unit square;}
  \item{5}{ Fairly accurate and fast, based on direct integration on the
    unit square from a new Gauss internal
    procedure (DEFAULT);}
   \item{6}{ Most accurate but slow, based on a cdfbvn
    procedure by Alan Genz (using results from Drezner, Z. and
    G.O. Wesolowsky, 1989. ``On the computation of the bivariate normal
    integral,'' Journal of Statist. Comput. Simul. 35: 101-107). See
    Appendix F.}
  
  
    Option 5 (the default) appears to be the best tradeoff between speed
    and accuracy currently available (and so this global should not be
    changed to anything other than 6, which is more accurate but much
    slower, unless you have a good reason to do so). However,
    fundamental progress remains to be made on methods of integrating
    the bivariate normal, as all currently available methods are
    innacurate and jump
  }
    \item{"EdirTol"}{Direction tolerance for CML
    convergence. Default=0.0001. Set to smaller values if most of your
    values of \eqn{T_i} or \eqn{X_i} are very small.
  }       
  \item{"EcdfTol"}{Tolerance for the lncdfbvn function (when Ecdfbvn=5,
    its default), with smaller calculated values truncated at the value
    of this global (DEFAULT=2.220446e-11). This can be any positive
    number, although lncdfbvn gets imprecise for small values. Only set
    to smaller values if you think you need the precision, such as if
    most of your values of \eqn{T_i} or \eqn{X_i} are very small.
  }
  \item{"EvTol"}{Numerical tolerance for the conditional variance
    calculation. Must be greater than 0; Default is \eqn{1e-322}
  }
  \item{"EdoML"}{do maximum likelihood (default); if 0 don't do maximum
    likelihood, using instead the values of \eqn{\phi} stored in EdoML.phi
    and vcphi in EdoML.vcphi.
  }
  \item{"doml.phi"}{If EdoML=1, this should include a vector of
    values of \eqn{\phi} and will be used instead of the output of the
    likelihood maximization. (This global is ignored unless EdoML=1.)
  }      
  \item{"doml.vc"}{EdoML=1, this should include a matrix of values of
    estimated variance matrix \eqn{V(\phi)} and will be used instead of the
    output of the likelihood maximization procedure. (This global is
    ignored unless EdoML=1.)
  }
  \item{"EdoSim"}{If 1 do simulations (default); if 0 don't do simulations; 
    if \eqn{-1} don't do simulations or compute the maxlik variance (use this
    option for computing conditional log-likelihood of eta's).
  }

    \item{"Eeta"}{Automatically includes \eqn{X_i} in the inputs Zb  and/or
    Zw. The actual inputs Zb and Zw  must be set to 1 if the default is
    changed. Using this global is better than explicitly including 
    \eqn{X_i} in the inputs, because eiread and eigraph will be ``aware'' of
    the contents of Zb and Zw. If you set this global, it is generally
    best to also set the priors EalphaB and EalphaW. See Chapter 9,
    and the parameterization in Equation 9.2 (page 170). Options
    include:}  

    \item{"Eeta=0"}{ excludes \eqn{X_i}, which is equivalent to setting
	 \eqn{\alpha^b=\alpha^w=0} (default).}
    
    \item{"Eeta=1"}{ sets zb=x, zw=1, which estimates \eqn{\alpha^b} and
	 fixes \eqn{\alpha^w=0}.}
    
    \item{"Eeta=2"}{ sets zb=1, zw=x, which estimates \eqn{\alpha^w} and
	 fixes \eqn{\alpha^b=0}}
    
    \item{"Eeta=3"}{sets zb=zw=x, which estimates \eqn{\alpha^b} and
	 \eqn{\alpha^w}}
    
    \item{"?"}{Set to a \eqn{4\times 1} vector with elements Eeta =
	 \eqn{\alpha^b\vert\alpha^w\vert} se \eqn{(\alpha^b)\vert} se
	 \eqn{(\alpha^w)} to fix \eqn{\alpha^b} and \eqn{\alpha^w},
	 and their standard deviations, during estimation.}
    
    \item{"?"}{Finally, set to a \eqn{3\times 1} vector \eqn{4\vert a\vert b} to set
	 zb=x and zw=1, to estimate \eqn{\alpha^b}, and fix \eqn{\alpha^w=a} and
    its standard error to b. Set to a \eqn{3\times 1} vector \eqn{5\vert
	 a\vert b} to set zb=1 and zw=x, to estimate \eqn{\alpha^w}, and fix 
    \eqn{\alpha^b=a} and its standard error to b.}

   \item{"eigraph.bvsmth"}{Smoothing parameter for nonparametric
    estimation; used only if Enonpar=1. Default=0.08. (The same
    parameter controls the nonparametric bivariate density estimation
    for diagnostic purposes in eigraph.) See Section 9.3.2.
  }
  \item{"EisChk"}{If 0 to do nothing (default); if 1 change lnir from the
    scalar mean importance ratio to a (\Esims\eqn{\star}Eisn)
    \eqn{\times rows(\eqn{\phi})+1} matrix containing the log of the importance ratio as the
    first column and normal simulations of \eqn{\tilde{\phi}} as the
    remaining columns. Also changes PhiSims from the mean and standard
    deviation of the posterior phi's to a \eqn{Esims \times}rows(\eqn{\phi})
    matrix of normal simulations of phi.}
  
  \item{"EiLliks"}{If 1 save (Esims\eqn{\times 1}) the log-likelihoods
    evaluated for each simulation; if 0 saves only the means of these
    likelihoods (default). These can be used for computing the marginal
    likelihood.
  }
  \item{"EisFac"}{Factor to multiply by estimated variance matrix in the
    normal approximation for use in importance sampling, or set to \eqn{-1}
    to use normal approximation only or \eqn{-2} to condition on the
    maximum posterior estimates. Adjust this, Eisn, or Eist if
    eiread's resamp larger than 15 or 20. If this is set too low,
    estimation variability will not be sufficient and your confidence
    intervals may be too narrow; it must be greater than zero and should
    probably be at least one. See Section 7.5. (Default=4).
  }
  \item{"Eisn"}{factor to multiply by Esims to compute the number of
    normals to draw before resampling. This is used to to try to get
    Esims samples from exact posterior. Increase this or change EisFac
    or _EisT if resamp is larger than 15 or 20. Default=10. See Section
    7.5.
  }
  \item{"Eist"}{If 0 (default) to use multivariate normal density to draw
    random numbers for initial approximation for importance sampling; or
    if greater than 2, use the multivariate Student \eqn{t} density, with
    degrees of freedom Eist. Use this, EisFac, or Eisn if resamp is
    larger than 15 or 20. See Section 7.5.
  }
  \item{"EmaxIter"}{Maximum number of iterations for CML. Default=500.}
  \item{"EnonEva"}{Number of nonparametric density evaluations for each
    tomography line (default=11). Only used if EnonPar=1.
  }
  \item{"EnonNumInt"}{Number of points to evaluate for numerical
    integration in computing the denominator for the bivariate kernel
    density (default=50). Only used if EnonPar=1.
  }
  \item{"EnonPar"}{If 0 do not run nonparametric model (default); if 1 run
    nonparametric model. (When choosing nonparametric estimation, only
    relevant options will be available under eigraph and eiread.) See
    Section 9.3.2.
  }
  \item{"EnumTol"}{Numerical tolerance. A homogeneous precinct is one
    for which \eqn{X_i<}EnumTol or \eqn{X_i>}(1- EnumTol). Default is
    0.0001. Set to smaller values if most of your values of \eqn{T_i} or 
    \eqn{X_i} are very small.
  }
  \item{"Erho"}{The first element is the standard deviation of normal
    prior on \eqn{\phi_5} for the correlation; set to 0 to fix \eqn{\phi_5} to
    a second element, Erho[2]; set to \eqn{-1} to estimate without a
    prior. Default=0.5. _Erho should be a scalar unless the first
    element is 0, in which case it should be a \eqn{2\times 1}  vector,
    where the second element is the value at which the \eqn{\phi_5}  is
    fixed (and not estimated). See Section 7.4.
  }          
  \item{"Eselect"}{Controls which observations are included in the
    estimation stage, including both likelihood maximization and
    importance sampling. All observations are included in the simulation
    stage unless you delete them from the data set before starting 
    \eqn{\mathfrak{E}I}. This allows users to base the truncated bivariate
    normal contours on a subset of observations that might be more
    representative (such as those for which \eqn{T_i} is not 0 or 1). Set
    to \eqn{p\times 1} vector to of 1's to include and 0's to exclude
    individual observations.
  }
  \item{"EselRnd"}{Set to scalar 1 to include all observations not
    already deleted by Eselect (default), or a scalar greater than 0
    and less than 1 to randomly select this fraction of observations in
    the estimation stage. This global is especially useful for speeding
    up estimation in very large datasets, since thousands of
    observations are not always needed for estimating $ \phi$. Since all
    observations will still be included in the simulation stage,
    precinct-level estimates of all quantities of interest will still be
    available. (If used with EI2, each iteration of EI includes a
    different randomly selected set of observations)}
  
  \item{"Esigma"}{Standard deviation of an underlying normal
    distribution, from which a half normal is constructed as a prior for
    both \eqn{\breve{\sigma}_b} and \eqn{\breve{\sigma}_w}. Note: the expected
    value under this prior is Esigma \eqn{\sqrt{2/\pi}\approx}
    Esigma0.8. Set to zero or negative for no prior. Default =0.5. See Section 7
  }
  \item{"Esims"}{Number of simulations. Default is 100.}
  
  \item{"Estval"}{For gradient methods: Scalar 1, use best guess
    starting values (default); or set to \eqn{k\times 1} vector of starting
    values. If Eeta[1]=0 (its default), \eqn{k=5} with elements guesses of
    \eqn{\phi}, that is on the scale of estimation. If you have starting
    values on the untruncated normal scale, 
    \eqn{\breve{\psi}=\{{\mathfrak{B}}^b,{\mathfrak{B}}^w,\sigma_b,\sigma_w,\rho\}},
    you can reparameterize as in this example:
    Estval=eireparinv(c(.5,.5,.2,.2,-.1)).If Eeta[1]=1, 2, 4, or 5, \eqn{k=6}; if Eeta=3, \eqn{k=7};
    and if covariates are used and nrow(Eeta)=4, then \eqn{k} is 5
    plus the number of covariates included, with Zb coming before Zw.

    For a grid search: Set Estval to scalar 0 (with 5 divisions per
    zoom), or to a scalar integer greater than or equal to 3 for a grid
    search with this number of divisions per zoom. (That is, the grid
    search procedure divides the parameter space into a number of
    divisions, evaluates the likelihood for every combination of values
    on all the parameters, chooses the region of highest likelihood,
    zooms in and repeats the procedure on the narrower parameter
    region. This continues until differences in the parameters differ by
    the global Edirtol.)}
   \item{"ei.vc"}{\eqn{M\times 2} matrix (\eqn{M\geq 1}), each row of which
    represents instructions for one attempt to compute an estimated
    positive definite variance matrix of \eqn{\phi}. The procedure exits
    after the first positive definite hessian is found. Options to
    include in various rows are: c(1, 0) the usual numerical hessian
    computation; c(1 , d) use usual hessian procedure and then adjust
    eigenvalues together so they are greater than d; c(2, f) use wide
    step lengths at f fraction falloff in the likelihood function;c(3, f) use quadradic
    approximation with falloff in likelihood function set at  f; c(4, 0)
    use a generalized inverse (to deal with singularity) and a
    generalized cholesky (to deal with non-positive definiteness) based
    on work in progress by Jeff Gill and Gary King; c(5, 0) use wide step
    lengths but check that the gradients for each are correct (and if
    necessary search for better ones); c(-1, 0) avoid the computation of
    the variance covariance matrix in case of non-positive definiteness
    and use the singular value decomposition for the multinomial normal
    sampling (i.e. EisT has to be set to 0). In order to use this
    option, also make sure to define relatively narrow upper and lower
    bounds of the parameters by using Ebounds. DEFAULT=matrix(c(1, 0, 4, 0, 2,
	 0.1, 2, 0.05, 3, 0.1, 1, 0.1, 1, 0.2),ncol=2,byrow=TRUE). The variance computation only
    very rarely gets beyond the second try.
  

    When the likelihood surface is normal (i.e., quadratic), which is
    true asymptotically, all options produce identical results. In
    practice, this procedure is useful for ensuring that a positive
    definite variance matrix can be found due to numerical, rather than
    theoretical or empirical, difficulties, as can happen when the mode
    of the truncated normal is far from the unit square due to
    imprecision in the function that computes the bivariate normal
    CDF. (Another, sometimes better, way to fix these numerical problems
    is to reduce the variances of the priors in Erho and Esigma.)
    Because importance sampling is used after this procedure, different
    values of the variance matrix can produce identical estimates of the
    quantities of interest. Be sure to verify that the simulations are
    being appropriately drawn from the estimated contours (see compare
    the right two figures in eigraph's tomogS).
  }
  \item{"betabs"}{Matrix of point estimates. Number of rows equal to the lenght of 
    either input vectors, x,t, or n; number of columns equal to number of
    simulation \code{Esims}. The prescint estimates from the aggregate
    data.}
  \item{"retcode"}{~~ If necessary, add description  ~~}
  
  \item{"phi"}{if EdoML =1, this should include a vector of values of
    \eqn{\phi} and will be used instead of the output of the likelihood
    maximization. (This global is ignored unless EdoML=1.)
  }
  \item{"loglik"}{if 1 save (Esims \eqn{ \times 1}) the log-likelihoods
    evaluated for each simulation; if 0 saves only the means of these
    likelihoods (default). These can be used for computing the marginal
    likelihood.
  }
  \item{"ghactual"}{~~ If necessary, add description  ~~}
  \item{"vcphi"}{EdoML=1, this should include a matrix of values of
    estimated variance matrix \eqn{V(\phi)} and will be used instead of the
    output of the likelihood maximization procedure. (This global is
    ignored unless EdoML=1.)
  }         
  
  }
  \item{"date"}{Current date}          

  
\references{\url{http://gking.harvard.edu/ei} }  

\author{ 
Gary King<\email{king@harvard.edu}}
\note{ 
See the documentation on line: \url{http://gking.harvard.edu/ei/docs}
 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

function(t,x,tvap,Zb, Zw,...)
{
 ###  local res,et,MLpsi,MLvc,betaBs,betaWs,tst,Eselect,flat;
  evbase <- eiset(t,x,tvap,Zb,Zw,...)  ##environment 
  drvdot <- match.call(expand.dots=TRUE)
  drv  <-  match.call(expand.dots=FALSE)

  n <- tvap
 
  evbase <- expanddots(drvdot,drv,evbase)
  
  ###	@ timing start @
  et   <- proc.time()
  param  <- ls(env=evbase)

  if(length(grep("Eprt", param)) && length(grep("Eversion", param))){
    Eversion <- get("Eversion", env=evbase)
    message(Eversion)
  }
##  print(EnonPar) 
  ### copy variables from evbase to local environment evei
  evei <- getEnvVar(evbase, environment())  ##environment 
 
  if(!scalzero(Eres) && vin(Eres,"titl") && Eprt >0)
    message(vread(Eres,"titl"))
 
  Eres <- add.to.Eres(Eres, round=1, evbase)
 
###testing and debugging            
###  lapply(param, function(x){
###    print(x)
###    print(get(x, env=evbase))})
###    checkinputs(t,x,n,Zb,Zw);
###   return(evbase)

  
###  /* verify inputs */
  if(as.logical(Echeck)){
 
    tst <- ""
    tst <- checkinputs(t,x,tvap,Zb, Zw,evbase)
    if(tst != ""){
      print(tst)
  ###    "----- EI Aborted -----";
      return(Eres)
    }
    if(Eprt>0){
      if(EnonPar)
        print("Inputs ok, beginning nonparametric estimation...")
      else
      print("Inputs ok, beginning preliminary estimation...")
    }
  }
  
 ###  /* augment _Eselect if _EselRnd<1 */
  assign("Eselect", Eselect, env=evbase);  ###$@ save existing value @
  Eselect0 <- Eselect
  Eselect <- matrix(1,nrow=rows(x),ncol=1)* as.vector(Eselect);
 
  if(EselRnd<1){
    vec <- runif(rows(x),min=0,max=1)
    vec <- matrix(vec)
    Eselect <- Eselect & (vec < EselRnd);
  }
 
 
###  /* nonparametric estimation */
  
  if(dbug==T) evglobal <<- evbase 
  if(EnonPar>=1){
       
    betaBs <- einonp(t,x, evbase);
    assign("betaBs", betaBs, env=evbase); 
  

    Eres <- add.to.Eres(Eres, round=2, evbase)  
    return(timing(et,Eprt,Eres,Eselect0))
  }
  ### /* parametric estimation: */
  
  ###/* eta influence on Zb,Zw */
  if(rows(Eeta)!=4){
  if(Eeta[1]==1 ||  Eeta[1]==4){
    Zb <- x;
    Zw <- 1;
  }else if(Eeta[1]==2 || Eeta[1]==5){
    Zb <- 1;
    Zw <- x;
  }else if(Eeta[1]==3){
    Zb <- Zw <- x;
 
  }
}

  ### /* set internal global */


 ### clearg _Ez;	@ n of covariates, incl. implied constant term for Zb|Zw @
  Ez <- as.matrix(c((cols(Zb)+1-scalone(Zb)),(cols(Zw)+1-scalone(Zw))));
  assign("Ez", Ez, env=evbase)

###  /* likelihood estimation */
  if(EdoML==1){
   
     lst  <- quadcml(x,Zb,Zw,t);
     MLpsi <- lst$Mlpsi
     MLvc <- lst$MLvc
     
     if(length(MLvc) <= 1 && is.na(MLvc)){
       Eres <- vput(Eres,MLpsi,"phi");
       Eres <- vput(Eres,MLvc,"vcphi");
      return(Eres);
     }
  }else{
    message("Skipping likelihood estimation..");
     Eres <- add.to.Eres(Eres, round=3, evbase)  
  
    MLpsi <- mlpsi <- EdoML.phi;
    MLvc <-  mlvc <- EdoML.vcphi;
    if(rows(mlvc)!= rows(mlpsi))
      stop("ei: EdoML.phi or EdoML.vcphi input error");
    
  }

 ###/* simulation */
  if(EdoSim==1){
  ###  {betaBs,betaWs} = psim1(T,X,tvap,Zb,Zw,MLpsi,MLvc);
     lst <- psim1(t,x,tvap,Zb,Zw,MLpsi,MLvc);
     betaBs <- lst$betaBs
     betaWs <- lst$betaWs
     Eres   <- vput(Eres,betaBs,"betaBs"); ###@ no need to save betaWs; see eiread @
     assign("Eres", Eres,env=evbase)
   }


 return(timing(et,Eprt,Eres,Eselect0));
 

}
}
\keyword{file}% at least one, from doc/KEYWORDS
\keyword{datasets}% __ONLY ONE__ keyword per line
