#' Ecological Inference Estimation
#'
#' \code{ei} is the main command in the package \code{EI}.  It gives
#' observation-level estimates (and various related statistics) of
#' \eqn{\beta_i^b} and \eqn{\beta_i^w} given variables \eqn{T_i} and \eqn{X_i}
#' (\eqn{i=1,...,n}) in this accounting identity: \eqn{T_i=\beta_i^b*X_i +
#' \beta_i^w*(1-X_i)}.  Results are stored in an \code{ei} object, that can be
#' read with \code{summary()} or \code{eiread()} and graphed in \code{plot()}.
#'
#' The \code{EI} algorithm is run using the \code{ei} command.  A summary of
#' the results can be seen graphically using \code{plot(ei.object)} or
#' numerically using \code{summary(ei.object)}.  Quantities of interest can be
#' calculated using \code{eiread(ei.object)}.
#'
#' @aliases EI ei
#' @param formula A formula of the form \eqn{t ~x} in the \eqn{2x2} case and
#' \eqn{cbind(col1,col2,...) ~ cbind(row1,row2,...)} in the RxC case.
#' @param total `total' is the name of the variable in the dataset that
#' contains the number of individuals in each unit
#' @param Zb \eqn{p} x \eqn{k^b} matrix of covariates or the name of covariates
#' in the dataset
#' @param Zw \eqn{p} x \eqn{k^w} matrix of covariates or the name of covariates
#' in the dataset
#' @param id `id' is the name of the variable in the dataset that identifies
#' the precinct. Used for `movie' and `movieD' plot functions.
#' @param data data frame that contains the variables that correspond to
#' formula.  If using covariates and data is specified, data should also
#' contain \code{Zb} and \code{Zw}.
#' @param erho The standard deviation of the normal prior on \eqn{\phi_5} for
#' the correlation. Numeric vector, used one at a time, in order. Default `c(.5, 3, 5)`.
#' @param esigma The standard deviation of an underlying normal distribution,
#' from which a half normal is constructed as a prior for both
#' \eqn{\breve{\sigma}_b} and \eqn{\breve{\sigma}_w}. Default \eqn{= 0.5}
#' @param ebeta Standard deviation of the "flat normal" prior on
#' \eqn{\breve{B}^b} and \eqn{\breve{B}^w}.  The flat normal prior is uniform
#' within the unit square and dropping outside the square according to the
#' normal distribution.  Set to zero for no prior. Setting to positive values
#' probabilistically keeps the estimated mode within the unit square.
#' Default\eqn{=0.5}
#' @param ealphab cols(Zb) x 2 matrix of means (in the first column) and
#' standard deviations (in the second) of an independent normal prior
#' distribution on elements of \eqn{\alpha^b}.  If you specify Zb, you should
#' probably specify a prior, at least with mean zero and some variance (default
#' is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^b}).
#' @param ealphaw cols(Zw) x 2 matrix of means (in the first column) and
#' standard deviations (in the second) of an independent normal prior
#' distribution on elements of \eqn{\alpha^w}.  If you specify Zw, you should
#' probably specify a prior, at least with mean zero and some variance (default
#' is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^w}).
#' @param truth A length(t) x 2 matrix of the true values of the quantities of
#' interest.
#' @param simulate default = TRUE:see documentation in \code{eiPack} for
#' options for RxC ei.
#' @param covariate see documentation in \code{eiPack} for options for RxC ei.
#' @param lambda1 default = 4:see documentation in \code{eiPack} for options
#' for RxC ei.
#' @param lambda2 default = 2:see documentation in \code{eiPack} for options
#' for RxC ei.
#' @param covariate.prior.list see documentation in \code{eiPack} for options
#' for RxC ei.
#' @param tune.list see documentation in \code{eiPack} for options for RxC ei.
#' @param start.list see documentation in \code{eiPack} for options for RxC ei.
#' @param sample default = 1000
#' @param thin default = 1
#' @param burnin default = 1000
#' @param verbose default = 0:see documentation in \code{eiPack} for options
#' for RxC ei.
#' @param ret.beta default = "r": see documentation in \code{eiPack} for
#' options for RxC ei.
#' @param ret.mcmc default = TRUE: see documentation in \code{eiPack} for
#' options for RxC ei.
#' @param usrfun see documentation in \code{eiPack} for options for RxC ei.
#'
#' @author Gary King <<email: king@@harvard.edu>> and Molly Roberts <<email:
#' molly.e.roberts@@gmail.com>>
#' @references Gary King (1997). A Solution to the Ecological Inference
#' Problem.  Princeton: Princeton University Press.
#'
#' @export
#' @return TODO
#'
#' @examples
#' data(sample_ei)
#' form <- t ~ x
#' dbuf <- ei(form, total = "n", data = sample_ei)
#' summary(dbuf)
ei <- function(formula, total = NULL, Zb = 1, Zw = 1, id = NA, data = NA,
               erho = .5, esigma = .5, ebeta = .5, ealphab = NA, ealphaw = NA,
               truth = NA, simulate = TRUE, covariate = NULL, lambda1 = 4,
               lambda2 = 2, covariate.prior.list = NULL, tune.list = NULL,
               start.list = NULL, sample = 1000, thin = 1, burnin = 1000,
               verbose = 0, ret.beta = "r", ret.mcmc = TRUE, usrfun = NULL) {
  # Extract formula
  dv <- terms.formula(formula)[[2]]
  iv <- terms.formula(formula)[[3]]
  t <- as.character(dv)
  x <- as.character(iv)
  n <- as.character(total)
  id <- as.character(id)

  if (length(dv) == 1) {
    cli::cli_progress_step("Running 2x2 ei")
    if (!simulate) {
      dbuf <- ei.estimate(t, x, n,
        id = id, data = data, Zb = Zb, Zw = Zw,
        erho = erho, esigma = esigma, ebeta = ebeta,
        ealphab = ealphab, ealphaw = ealphaw, truth = truth
      )
      return(dbuf)
    }
    if (simulate) {
      dbuf <- NULL
      i <- 1
      while (i <= length(erho) & is.null(dbuf)) {
        try(
          {
            dbuf <- ei.estimate(t, x, n,
              id = id,
              data = data, Zb = Zb, Zw = Zw, erho = erho[i],
              esigma = esigma, ebeta = ebeta,
              ealphab = ealphab, ealphaw = ealphaw,
              truth = truth
            )
          },
          silent = TRUE
        )
        i <- i + 1
      }
      if (is.null(dbuf)) {
        cli::cli_abort("{.fn ei.estimate} did not converge. Try a different value of {.arg erho}.")
      }
      cli::cli_progress_done()
      dbuf.sim <- ei.sim(dbuf)
      return(dbuf.sim)
    }
  }

  if (length(dv) > 1) {
    cli::cli_progress_step("Running eiRxC")
    # If the table is RxC use eiRxC
    dbuf <- ei.MD.bayes(formula,
      data = data, total = total, covariate = covariate,
      lambda1 = lambda1, lambda2 = lambda2,
      covariate.prior.list = covariate.prior.list,
      tune.list = tune.list, start.list = start.list,
      sample = sample, thin = thin, burnin = burnin,
      verbose = verbose, ret.beta = ret.beta, ret.mcmc = ret.mcmc,
      usrfun = usrfun
    )
    dbuf$data <- data
    dbuf$total <- n
    dbuf$formula <- formula
    class(dbuf) <- "ei"
    cli::cli_progress_done()
    return(dbuf)
  }
}

ei.estimate <- function(t, x, n, id, Zb = 1, Zw = 1, data = NA, erho = .5,
                        esigma = .5, ebeta = .5, ealphab = NA, ealphaw = NA,
                        truth = NA, Rfun = 2, precision = 4) {

  # Check to make sure data is not null
  if (!missing(data)) {
    if (is.character(t)) t <- data[[t]]
    if (is.character(x)) x <- data[[x]]
    if (is.character(n)) n <- data[[n]]
    if (is.character(Zb)) Zb <- data[[Zb]]
    if (is.character(Zw)) Zw <- data[[Zw]]
    if (is.character(id)) id <- data[[id]]
  }

  Zb <- as.matrix(Zb)
  Zw <- as.matrix(Zw)

  # If there are no covariates, run simplest R function
  if (dim(Zb)[1] == 1 & Zb[1, 1] == 1 & dim(Zw)[1] == 1 & Zw[1, 1] == 1) {
    Rfun <- 5
  }
  if (dim(Zb)[1] == 1 & Zb[1, 1] == 1) {
    Zb <- as.matrix(rep(1, length(x)))
  }
  if (dim(Zw)[1] == 1 & Zw[1, 1] == 1) {
    Zw <- as.matrix(rep(1, length(x)))
  }

  # Extract the number of covariates
  numb <- dim(Zb)[2]
  numw <- dim(Zw)[2]

  # Starting values
  start <- c(0, 0, -1.2, -1.2, 0, rep(0, numb + numw))

  cli::cli_alert_info("Maximizing likelihood")

  solution <- optim(start, like,
    y = t, x = x, n = n, Zb = Zb,
    Zw = Zw, numb = numb, erho = erho, esigma = esigma,
    ebeta = ebeta, ealphab = ealphab, ealphaw = ealphaw, Rfun = Rfun,
    hessian = TRUE,
    method = "BFGS"
  )

  # Find values of the Hessian that are 0 or 1.
  covs <- as.logical(ifelse(diag(solution$hessian) == 0 |
    diag(solution$hessian) == 1, 0, 1))
  hessian <- solution$hessian[covs, covs]
  output <- list(
    phi = solution$par,
    hessian = solution$hessian, hessianC = hessian,
    erho = erho, esigma = esigma,
    ebeta = ebeta, ealphab = ealphab, ealphaw = ealphaw, numb = numb,
    x = x, t = t, n = n,
    Zb = Zb, Zw = Zw,
    truth = truth, precision = precision, covs = covs, Rfun = Rfun, id = id
  )

  class(output) <- "ei"
  output
}


#' @noRd
#' @export
print.ei <- function(x, ...)
{
  cat("ei object \n")
}