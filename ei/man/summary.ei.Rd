\name{summary.ei}
\alias{summary.ei}

\title{Summary and analysis of relevant computational results for the procedure ei.}

\description{
  Extract information from \sQuote{dbuf} the output of running ei and
  further calculate parameters of interest. The list \sQuote{dbuf} is
  output from procedure ei or ei2 and, if a second argument
  (\sQuote{name}) is provided, it will return the item you desire
  to be extracted or computed from \sQuote{dbuf} (see the list below),
  It is a wrapper around the internal procedure \code{eiread}. 
}

\usage{summary.ei(object,...)}

\arguments{
  Choose from the following items for name. The possibilities include
  items stored as is in dbuf (and so could also be retrieved as
  dbuf[["name"]],special items that can be added to
  \sQuote{dbuf} for extra information. In fact, other items can be 
  computed by required some additional information. For example,
  those items that use true values from individual-level data will only work
  if you have added this information into the list \sQuote{dbuf}; if 
  \sQuote{truth} is a matrix of p (number of precincts rows)
  and two columns corresponding to betaB and betaW, append it to the final computation:
  as dbuf[["truth"]] <- truth, so that those additional items  
  may be computed. From the perspective of the user, items in each of these categories
  can be treated identically, as all calculations are automatic.
  These are useful for verifying the method if individual-level data are
  available, but they are not useful for most ecological inference
  applications, where this information is not available.
  The options are included here only because I needed them while writing
  the book.)
  
  When used with an EI2 output data buffer (i.e., for  2\eqn{\times} C
  tables), summary.ei (or the internal procedure eiread) uses use
  the mean posterior estimate for $ X_i$ for some items;
  the correct multiply imputed values of $ x_i$ (x2) are used only where specifically noted.

  \item{object}{String or scalar to add to the prompt-messages. Default is NA
    that will not change the prompt-messages.}
  
  \item{...}{
    By default, if no second argument is provided, the
    function call \code{summary(dbuf)} will return \sQuote{Paggs}.
    If a second argument is provided then it needs to be one of the
    following options (note that the arguments are case insensitive):
  
\itemize{ 

   \item{EalphaB}{
    value of this global (priors on \eqn{\alpha^b})}

   \item{EalphaW}{
    value of this global (priors on \eqn{\alpha^w})}

  \item{Ebeta}{
   %\mathfrak{B}
    value of this global (priors on \eqn{\breve{\strong{B}}}).}

  \item{Ebounds}{
    value of this global (bounds for \code{nmlinb}, not quantities of interest).}

  \item{Ecdfbvn}{
    value of this global (method of calculating CDF of the bivariate
    normal).}
  
  \item{EdirTol}{
    value of this global (\sQuote{optim} convergence tolerance)}

  \item{EdoML}{
    value of this global (do maxlik)}

  \item{EdoML.phi}{
    value of this global (input \eqn{\phi}'s). Only relevant if \code{EdoML=0}.}

  \item{EdoML.vcphi}{
    value of this global (input variance-covariance of \eqn{\phi}'s). Only relevant if \code{EdoML=0}.}

  \item{Eeta}{
    value of this global (expanded model specifications; see EIBAD Chapter 9).}

  \item{EI.bma.prior}{
    value of this global (prior model probabilities for Bayesian Model Averaging).}

  \item{EI.vc}{
    value of this global (variance matrix computation).}

  \item{EImodels.save}{
    value of this global (file name for eimodels.run).}

  \item{EIgraph.bvsmth}{
    value of this global (smoothing parameter for nonparametric density
    estimation).}
  
  \item{EisChk}{
    value of this global (check importance sampling).}
  
  \item{EisFac}{
    value of this global (number to multiply by variance matrix
    by in importance sampling or $ -1$ for normal approximation).}

  \item{EisN}{
    value of this global (first stage importance sampling factor).}

  \item{EisT}{
    value of this global (multivariate $ t$ or normal for importance sampling).}

  \item{EmaxIter}{
    value of this global (maximum iterations for \sQuote{nmlinb}).}

  \item{EnonEval}{
    value of this global (number of nonparametric density evaluations
    for each tomography line). Only relevant if \code{Enonpar=1}.}
  
  \item{EnonNumInt}{
    value of this global (number of points to evaluate for numerical
    integration in computing the denominator for the bivariate kernel
    density). Only relevant if \code{Enonpar=1}.}
  
  \item{EnonPar}{
    value of this global (0, parametric or 1, nonparametric, estimation).}

  \item{EnumTol}{
    value of this global (numerical tolerance for homogeneous and unanimous precincts).}

  \item{Erho}{
    value of this global (prior on \eqn{\rho}). See EIBAD Section 7.4.}

  \item{Eselect}{
    value of this global (vector of 0/1 to delete/select observations during likelihood stage, or 1 to select all).}

  \item{EselRnd}{
    value of this global (fraction of observations to select randomly).}

  \item{Esigma}{
    value of this global (priors on \eqn{\breve{\sigma}_b} and \eqn{ 
    \breve{\sigma}_w}). See EIBAD Section 7.4.}

  \item{Esims}{
    value of this global (number of simulations). See EIBAD Appendix F.}

  \item{Estval}{
    value of this global when ei was run (starting values).}

  \item{undern}{
    p \eqn{\times} 1: the original variable $ N_i$ from the first stage EI analysis. This works only for EI2.}

  \item{undert}{
    p \eqn{\times} 1: the original variable $ T_i$ from the first stage EI analysis. This works only for EI2.}

  \item{EvTol}{
    value of this global (numerical tolerance for the conditional variance).}

  \item{underx}{
    p \eqn{\times} 1: the original variable $ X_i$ from the first stage EI analysis. This works only for EI2.}

  \item{Ez}{
    2 \eqn{\times} 1: number of covariates (see EIBAD Chapter 9), including implied constant term for Zb\eqn{\vert}Zw; also clears and sets this in global memory.}
  
  \item{ABounds}{
     2 \eqn{\times} 2: aggregate bounds rows:lower,upper;
    columns:betab,betaw. See EIBAD Chapter 5.}

  \item{ABounds2}{
     2 \eqn{\times} 2: aggregate bounds for the \eqn{\lambda}'s from EI2,
    rows:lower,upper; columns:lambdaB,lambdaW. See Chapter 5.}

  \item{AggBias}{
    4 \eqn{\times} 2: regressions of true \eqn{\beta_i^b} and \eqn{\beta_i^w} on a
    constant and $ X_i$. Requires truth to have been added to
    \sQuote{dbuf}, the output of running ei. See EIBAD Table 11.1, page
    219.}

  \item{aggs}{
    Esims \eqn{\times} 2: simulations of district-level quantities of
    interest, \eqn{\hat{B}^b\sim\hat{B}}. See Section 8.3.
%(Uses x2 simulations when used with ei2).
}

  \item{AggTruth}{
    2 \eqn{\times} 1: true district-level \eqn{B^b} and \eqn{B^w}.
    Requires truth to have been addded to \sQuote{dbuf} the output of
    running ei.See EIBAD Chapter 2.}

  \item{beta}{
    p \eqn{\times} 2 point estimates of \eqn{\beta_i^b} (in the first column)
    and \eqn{\beta_i^w} (in the second) based on their mean posterior.
    See Section 8.2. 
%(Uses x2 simulations when used with ei2.)
}

  \item{betaB}{
    p \eqn{\times} 1 point estimates of \eqn{\beta_i^b} based on its mean
    posterior. See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2.)
}

  \item{betaBs}{
    p \eqn{\times} Esims: simulations of \eqn{\beta_i^b}.
    See EIBAD Chapter 8. 
%(Uses x2 simulations when used with ei2).
    }

  \item{betaW}{
    p \eqn{\times} 1 point estimate of \eqn{\beta_i^w} based on its mean
    posterior. See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
   }

  \item{betaWs}{
    p \eqn{\times} Esims: simulations of \eqn{\beta_i^w}.
    See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
  }

  \item{bounds}{
    p \eqn{\times} 4: bounds on \eqn{\beta_i^b} and \eqn{\beta_i^w},
    lowerB, upperB, lowerW, upperW. See EIBAD Chapter 5.}

 \item{bounds2}{
   p \eqn{\times} 4: bounds on \eqn{\lambda_i^b} and \eqn{\lambda_i^w},
   lowerB, upperB, lowerW, upperW. See EIBAD Chapter 5, eqns 5.4.}

  \item{checkR}{
    nrow(\eqn{\phi}) \eqn{\times} 2: matrix with rows corresponding to \eqn{\phi},
    columns corresponding to slightly less, more (by amount EdirTol)
    than the MLEs,and each element indicating that the CDFBVN function
    is sufficiently precise (when 1) and insufficiently precise (when 0).
    This function calculates the R portion of the likelihood function
    to make the comparison (Equation 6.15, p.104). See option R.}

  \item{CI50b}{
    p \eqn{\times} 2: lower, upper 50\% confidence intervals for \eqn{\beta_i^b}.
	 See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
 }

  \item{CI50w}{
    p \eqn{\times} 2: lower, upper 50\% confidence intervals for 
    \eqn{\beta_i^w}. See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
}

  \item{CI80b}{
    p \eqn{\times} 2: lower, upper 80\% confidence intervals for \eqn{\beta_i^b}. 
	See EIBAD Section 8.2. %(Uses x2 simulations when used with ei2).
}

  \item{CI80bw}{
    p \eqn{\times} 4: lowerB, upperB, lowerW, upperW 80\%
    confidence intervals for \eqn{\beta_i^b} and \eqn{\beta_i^w}.
    See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
}

  \item{CI80w}{
    p \eqn{\times} 2: lower, upper 80\% confidence intervals
    for \eqn{\beta_i^w}. See Section 8.2. 
%(Uses x2 simulations when used with ei2).
}

  \item{CI95b}{
    p\eqn{\times} 2: lower, upper 95\% confidence intervals for \eqn{\beta_i^b}. See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
}

  \item{CI95bw}{
    p\eqn{\times} 4: lowerB, upperB, lowerW, upperW
    95\% confidence intervals for \eqn{\beta_i^b} and \eqn{\beta_i^w}.
    See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
}

  \item{CI95w}{
    p \eqn{\times} 2: lower, upper 95\% confidence intervals
    for \eqn{\beta_i^w}. See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
}

  \item{coverage}{
    2 \eqn{\times} 4: confidence interval coverage; percent of true values
    within the 50\% and 80\% confidence intervals:
    50b, 80b, 50w, 80w (1st row = means, 2nd = weighted
    means).Requires truth to have been added to \sQuote{dbuf} the output
    of ei. 
%(Uses x2 simulations when used with ei2).
}

  \item{CsbetaB}{
    p \eqn{\times} 1: confidence interval-based standard error of \eqn{\beta_i^b}.
%(Uses x2 simulations when used with ei2).
}

  \item{CsbetaW}{
    p \eqn{\times} 1 confidence interval-based standard error of
    \eqn{\beta_i^w}. 
%(Uses x2 simulations when used with ei2).
}

  \item{DataSet}{
    Zb, Zw, x, t, used for input to eiloglik.
    If Eselect is a scalar less than 1 (for random selection of cases).
    The global EselRnd is ignored.}

  \item{date}{
    %\mathfrak{E}
    a string containing the date and time at which execution completed,
    as well as the version number and date of the program
    {\strong{E}I} that created the input data buffer.}

  \item{double}{
    2 \eqn{\times} 1 coefficients from a double regression.
    Requires an EI2 data buffer as input. See EIBAD Section 4.3.}

  \item{EaggBias}{
    4 \eqn{\times} 2: regressions of estimated \eqn{\beta_i^b}
    and \eqn{\beta_i^w} on a constant term and $ X_i$.
    First row: coefficients, second row: standard errors. See EIBAD Section 9.2.4.}

  \item{etaC}{
    2\eqn{\times} 1 coefficients implied by global Eeta.}

  \item{etaS}{
    2\eqn{\times} 1 standard errors implied by global Eeta.}

  \item{ExpVarCI}{
    100\eqn{\times} 4: 80\% confidence intervals for $ T_i$
    given $ X_i$. X, 20\%CI, mean, 80\%CI of sims from
    P(T_i\eqn{\vert} X), where $ X$ in this context is 100 numbers
    equally spaced between 0 and 1. See EIBAD Section 8.5.
%(Uses x2 simulations when used with ei2).
}

  \item{ExpVarCI0}{
    p\eqn{\times} 4: 80\% confidence intervals for $ T_i$ given
    the observed values of $ X_i$ and/or Zb and Zw when applicable.
    T, 20\%CI, mean, 80\%CI of sims from
    P(T_i\eqn{\vert} X, Zb, Zw) 
%(Uses x2 simulations when used with ei2).
}

  \item{ExpVarCIs}{
    100\eqn{\times} 4: same as expvarci, but smoothed with LOESS.
    Used for eigraph's xtfit. See EIBAD Section 8.5. 
%(Uses x2 simulations when used with ei2).
}

  \item{GEbw}{
    p\eqn{\times} 3: \eqn{\beta_i^b, \beta_i^w, }Nsims.
    Point estimates of \eqn{\beta_i^b} and \eqn{\beta_i^w},
    based on the mean posterior under the prior that \eqn{\beta_i^b\geq\beta_i^w}.
     Use this option (instead of beta)
    if you are reasonably certain that \eqn{\beta_i^b\geq\beta_i^w}.
    This procedure is based on simulations for which the inequality holds,
    and thus also reports Nsims, the number of simulations on which
    each estimate is based (If Nsims is not close to Esims,
    you may wish to question your assumption or increase Esims).}

  \item{GEbwa}{
    p\eqn{\times} 2: \eqn{B^{b}\sim B^{w}}. Aggregate level mean posterior
    estimates, based on simulations where
    \eqn{\beta_{i}^{b}\geq\beta_{i}^{w}}.See GEbw for more information.}

  \item{GEwb}{
    p\eqn{\times} 3: \eqn{\beta_i^b,\beta_i^w,}Nsims.
    Point estimates of \eqn{\beta_i^b} and \eqn{\beta_i^w}, based on the
    mean posterior under the prior that \eqn{\beta_i^b\leq\beta_i^w}.
    Use this option (instead of beta) if you are reasonably certain that
    \eqn{\beta_i^b\leq\beta_i^w}. This procedure is based on simulations for
    which the inequality holds, and thus also reports Nsims, the number
    of simulations on which each estimate is based
    (If Nsims is not close to Esims, you may wish to question your assumption or increase Esims).}

  \item{GEwba}{
    p\eqn{\times} 2: \eqn{B^b\sim B^w}. Aggregate level mean posterior
    estimates, based on simulations where \eqn{\beta_i^b\leq\beta_i^w}.
    See GEwb for more information.}

  \item{GhActual}{
    value of the output global GhActual, the row of ei.vc for
    which a positive definite variance matrix was found.}

  \item{Goodman}{
    2\eqn{\times} 2: row 1: Goodman's Regression coefficients,
    row 2: standard errors. See EIBAD Section 3.1.}

  \item{LogLik}{
   value of log-likelihood at the maximum (unnormalized).}
 
 \item{LogLikS}{
   value of log-likelihood at the maximum (unnormalized) for each observation.}

  \item{lnIR}{
    If EisChk=1, this is a \eqn{Esims\ast Eisn \times nrow(\phi)+1} matrix ,
    containing the log of the importance ratio as the first column
    and normal simulations of \eqn{t(\tilde{\phi})} as the remaining columns.
    If EisChk=0, this is the scalar mean importance ratio (equivalent to
    meanIR below).}
 
   \item{LLikSims}{
    If EiLlikS=1, this is Esims\eqn{\times}1 vector of
    log-likelihood values for each simulation;
    otherwise it is a scalar mean of these.}

  \item{Maggs}{
    2\eqn{\times} 1: point estimate of 2 district-level parameters,
    \eqn{\hat{B}^b} and \eqn{\hat{B}^w}: meanc(aggs).
    See EIBAD Section 8.3. 
%(Uses x2 simulations when used with ei2).
}

  \item{MeanIR}{
    scalar log of the mean importance ratio.}

  \item{mpPsiu}{
    Mean Posterior of \eqn{\breve{\psi}} (rather than MLEs).}

  \item{N}{
    p\eqn{\times} 1: number of individual elements in precinct i, N_i, an input to ei.}

  \item{Nb}{
    p\eqn{\times} 1: denominator of x and t, equal to \eqn{x\ast n}, \eqn{N_i^{b}}
    (e.g., number of blacks in the voting age population).}

%  \item{Nb2}{
%    for ei2 data buffers only, p\eqn{\times} Esims: denominator of x2 and V,
%    equal to \eqn{x2\ast n}, \eqn{N_i^{b}}
%    (e.g., number of blacks in the voting age population).}

  \item{Nt}{
    p\eqn{\times} 1: numerator of t, equal to \eqn{t\ast n},
    \eqn{N^t_i} (e.g., number of people who turnout).}

  \item{NbN}{
    p\eqn{\times} 1: \eqn{{N_i^{bN}}} (e.g., number of blacks who don't vote).
    Requires truth to have been added to \sQuote{dbuf}, the output of ei.
    See EIBAD Chapter 2.}

  \item{NbT}{
    p\eqn{\times} 1: \eqn{N_i^{bT}} (e.g., number of blacks who vote).
    Requires truth to have been added to \sQuote{dbuf}, the output of ei.
    See EIBAD Chapter 2.}

  \item{Neighbor}{
    Freedman et al.'s neighborhood model point estimates
    (i.e., assuming \eqn{ \beta_i^b=\beta_i^w}, with implied standard errors of zero).}

  \item{nobs}{
    scalar: number of observations,  p.}

  \item{Nw}{
    p\eqn{\times} 1: equal to (1-x)\eqn{\ast} n, \eqn{N_i^{w}}
    (e.g., number of whites in the voting age population).}

%\item{Nw2}{
%  for ei2 data buffers only, p\eqn{\times} Esims:
%  equal to \eqn{(1-x2)\ast n}, \eqn{N_i^{w}} 
% (e.g., number of whites in the voting age population).}

  \item{NwN}{
    p\eqn{\times} 1: \eqn{N_i^{wN}} (e.g., number of whites who don't vote).
    Requires truth to have been added to \sQuote{dbuf} the output of ei.
    See EIBAD Chapter 2.}

  \item{NwT}{
    p\eqn{\times} 1: \eqn{N_i^{wT}} (e.g., number of whites who Turnout).
    Requires truth to have been added to \sQuote{dbuf} the output of ei.
    See EIBAD Chapter 2.}

  \item{Paggs}{
    2\eqn{\times} 2: row 1: \eqn{\hat{B}^b} and \eqn{\hat{B}^w}; row 2: standard
    errors. See EIBAD Section 8.3. 
%(Uses x2 simulations when used with ei2).
}

  \item{Palmquist}{
    scalar: Palmquist's Inflation Factor. See EIBAD Equation 3.14, page 52.}

  \item{ParNames}{
    character vector of names for \eqn{\phi}.}

  \item{phi}{
    maximum posterior estimates from \sQuote{nmlinb}.}

  \item{PhiSims}{
    If EisChk=1, PhiSims is a Esims \eqn{\times}nrow(\eqn{\phi}) matrix of random
    simulations of \eqn{\phi}; otherwise, it is a rows(\eqn{\phi}) \eqn{\times} 2
    matrix of the means (in the first column) and standard deviations
    (in the second column) of the simulations
    (which are the mean and standard deviations of the posterior
    distribution of \eqn{\phi}). See EIBAD Section 8.2.}

  \item{Pphi}{
    2\eqn{\times} 5 maximum posterior estimates. row 1: \eqn{\phi}, row 2:
    standard errors. See EIBAD Chapter 7.}

  \item{psi}{
    reparameterized \eqn{\phi} into ultimate truncated scale.
    See EIBAD Section 6.2.2.}

  \item{PsiTruth}{
    5\eqn{\times} 1: true values of \eqn{\psi}
    (i.e., on truncated scale). Requires truth to have been added to
    \sQuote{dbuf}, which is the output ei. See EIBAD Table 10.3, page 207.}

  \item{psiu}{
    \eqn{\breve{\psi}}, which was reparameterized from
    \eqn{\phi} into untruncated scale. See EIBAD Equation 7.4, page 136.
  }
  \item{R}{
    %\mathfrak{B} 
    The sum of the log of the volume above the unit square under the
    bivariate normal, \eqn{R(\breve{\strong{B}},\breve{\Sigma})}.
    This is the last piece of the likelihood function. See checkR.}

  \item{Ri}{
  %\mathfrak{B} 
  The log of the volume above the unit square under
  the bivariate normal, \eqn{R(\breve{\strong{B}},\breve{\Sigma})}
  for each $ i$. This is the last piece of the likelihood function,
  and will differ over $ i$ only if covariates are included. See checkR.}

\item{resamp}{
  number of resampling tries. This number will range between 1
  and Esims and is better if small. If it is greater than 15-20,
  you can try adjusting Eisn or EisFac and rerunning ei.}

\item{RetCode}{Non applicable.}

\item{RNbetaBs}{
  p\eqn{\times} Esims: randomly horizontally permuted
  simulations of \eqn{\beta_{i}^{b}}. This is essentially equivalent
  to betaBs except that it randomly permutes estimation
  variation also. 
%(Uses x2 simulations when used with ei2).
}

\item{RNbetaWs}{
   p\eqn{\times} Esims: randomly horizontally permuted simulations of
  \eqn{\beta_i^w}. This is essentially equivalent to betaWs except
  that it randomly permutes estimation variation also.
  
%(Uses x2 simulations when used with ei2).
}

\item{sbetaB}{
  p\eqn{\times} 1 standard error for the estimate of \eqn{\beta_i^b},
  based on the standard deviation of its posterior.
  See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
}

\item{sbetaW}{
  p\eqn{\times} 1 standard error for the estimate of \eqn{\beta_i^w},
  based on the standard deviation of its posterior.
  See EIBAD Section 8.2. 
%(Uses x2 simulations when used with ei2).
}

\item{STbetaBs}{
  p\eqn{\times} Esims: SORTED simulations of \eqn{\beta_i^b}
  (e.g., the 80\% confidence interval lower bound is
  STbetaBs[int(0.1* Esims)]). See EIBAD Section 8.2.
  
%(Uses x2 simulations when used with ei2).
}

\item{STbetaWs}{
  p\eqn{\times} Esims: SORTED simulations of \eqn{\beta_i^w}
  (e.g., the 80\% confidence interval upper bound is
  STbetaWs[int(0.9* Esims)]).
  See EIBAD Section 8.2. (Uses x2 simulations when used with ei2).
}

\item{sum}{
  prints a summary of district-level information.}

\item{t}{
  p\eqn{\times} 1: outcome variable proportion, $ T_i$
  (e.g., turnout); input to ei. If dbuf is the output from ei2,
  then this option gives $ V_i$ (e.g., Democratic fraction of the major
  party vote), an input to ei2.}

\item{Thomsen}{
  2\eqn{\times} 1: Estimates of \eqn{B^b\vert B^w} from Thomsen's Ecological Logit Model.}

\item{titl}{
  string: a title with descriptive information.
  Must have been added to \sQuote{dbuf} the output of running  ei.}

\item{truPtile}{
  p x 2: percentile of sorted simulates at which the true value
  falls for \eqn{\beta_i^b} and \eqn{\beta_i^w}.
  Requires truth to have been added  to \sQuote{dbuf}, which is the
  output of ei. See EIBAD Figure 10.7, page 213.
% (Uses x2 simulations when used with ei2).
}


\item{truth}{
  p\eqn{\times} 2: true values of the precinct-level quantities of interest
  \eqn{\beta_i^b\sim\beta_i^w}. Must have been added to \sQuote{dbuf}, the output of running ei.}

\item{truthB}{
  truth[,1], the true \eqn{\beta_i^b}.
  Requires truth to have been added to \sQuote{dbuf} the result of running ei.}

\item{truthW}{
  truth[,2], the true \eqn{\beta_i^w}.
  Requires truth to have been added to \sQuote{dbuf} the result of running ei.}

\item{tsims}{
  100\eqn{\times}Esims +1: simulations of $ T_i$ given $ X_i$.
  rows correspond to 100 values of $ X_i$ equally spaced between zero
  and one, columns are $ X_i$ and sorted simulations of \eqn{T_i\vert X_i}.
  See EIBAD Section 8.5. 
%(Uses x2 simulations when used with ei2).
}

\item{tsims0}{
  p\eqn{\times} Esims +1: simulations of $ T_i$ given observed
  values of $ X_i$ and/or Zb and Zw. The first column is the actual
  values of $ T_i$ 
%(Uses x2 simulations when used with ei2).
}

\item{VCaggs}{
  2\eqn{\times} 2: variance matrix of 2 district-level parameters,
  \eqn{\hat{B}^b} and \eqn{\hat{B}^w}.
  See EIBAD Section 8.3. 
%(Uses x2 simulations when used with ei2).
}

\item{VCphi}{
  global variance matrix of coefficients phi from gvc().
  If \eqn{\texttt{EI.vc}} is set to {-1 0 }, then it returns
  the inverse of the variance matrix.}

\item{x}{
  p\eqn{\times} 1: explanatory variable proportion,
  $ X_i$ (e.g., black voting age population); input to ei.
  (for ei2, this is the mean posterior estimate of $ x_i$,
  as, e.g., the black fraction of those voting).}

%\item{x2}{
%  for ei2 data buffers only, p\eqn{\times}\texttt{Esims}:
%  simulations of the explanatory variable proportion,
%  $ x_i$ (e.g., black fraction of voters).}

%\item{x2rn}{
%  for ei2 data buffers only, p\eqn{\times}\texttt{Esims}:
%  horizontally randomly permuted simulations of the explanatory variable
%  proportion, $ x_i$ (e.g., black fraction of voters).
%  (This is useful because x2 has only EI2.m unique columns.)} 

\item{Zb}{
  matrix of covariates for \eqn{\beta_i^b} or 1 for none;
  as affected by Eeta; input to ei. See EIBAD Section 9.2.1.}

\item{Zw}{
  matrix of covariates for \eqn{\beta_i^w} or 1 for none;
  as affected by Eeta; input to ei. See EIBAD Section 9.2.1.} 
}
}
}

\details{Allows users to obtain relevant parameters and information from
  running ei, by inputing the \sQuote{dbuf} output of ei and the name of
  the parameter or method to derive from the analysis of \sQuote{dbuf}.
  Same as \code{v = ei:::eiread(dbuf,"name")} but eiread is internal to
  the R-package.\\
  Global variables relevant to this procedure are:
\itemize{
\item{Eprt}{
  if \eqn{\geq} 1, in addition to outputting results in v,
  print a display of the item to the screen
  (available for options whose description begins
  with an asterisk, below) (default); 0 do not print.}

\item{EIMetaR}{
  If dbuf is a meta-data buffer (i.e., from the ei2.mta
  output global of ei2; the output buffer from eimodels.def
  and eimodels.run), this global denotes which of the imputed
  data buffers (in case of ei2.mta) or which of the stored model
  (in case of eimodels.def and eimodels.run) in dbuf should be
  accessed when running summary (i.e. eiread) and plot (i.e. eigraph).
  Default=1.} 
}
}

\value{Can be a list,  a vector, matrix or string}

\references{\url{http://gking.harvard.edu/ei}
 
  \emph{A Solution to the Ecological Inference Problem: Reconstructing
  Individual Behavior from Aggregate Data}, Gary King, Princeton
University Press, 1997 (EIBAD)
}

\author{Gary King, \email{king@harvard.edu}}

\note{Elena Villalon, \email{evillalon@iq.harvard.edu},
  translated the original Gauss code into R language.}

\seealso{See demos; for example: demo(eip.default)}

