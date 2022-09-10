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
#' the correlation. Numeric vector, used one at a time, in order. Default `c(.5, 3, 5, .1, 10)`.
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
#' @param ndraws integer. The number of draws. Default is 99.
#' @param nsims integer. The number of simulations within each draw. Default is 100.
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
#' @return `ei` object
#' @examples
#' data(sample_ei)
#' form <- t ~ x
#' dbuf <- ei(form, total = "n", data = sample_ei)
#' summary(dbuf)
ei <- function(formula, total = NULL, Zb = 1, Zw = 1, id = NA, data,
               erho = c(.5, 3, 5, .1, 10), esigma = .5, ebeta = .5, ealphab = NA, ealphaw = NA,
               truth = NA, simulate = TRUE, ndraws = 99, nsims = 100,
               covariate = NULL, lambda1 = 4,
               lambda2 = 2, covariate.prior.list = NULL, tune.list = NULL,
               start.list = NULL, sample = 1000, thin = 1, burnin = 1000,
               verbose = 0, ret.beta = "r", ret.mcmc = TRUE, usrfun = NULL) {

  if (missing(formula)) {
    cli::cli_abort('{.arg formula} is required.')
  }

  # Extract formula
  dv <- terms.formula(formula)[[2]]
  iv <- terms.formula(formula)[[3]]
  t <- as.character(dv)
  x <- as.character(iv)
  n <- as.character(total)
  id <- as.character(id)

  if (missing(data)) {
    cli::cli_abort('{.arg data} is required.')
  }
  data <- as.data.frame(data)

  if (simulate) {
    if (!is.numeric(nsims)) cli::cli_abort('{.arg nsims} must be {.cls numeric}.')
    if (!is.numeric(ndraws)) cli::cli_abort('{.arg ndraws} must be {.cls numeric}.')
  }

  if (length(dv) == 1) {
    cli::cli_progress_step("Running 2x2 ei")
    dbuf <- NULL
    i <- 1
    while (i <= length(erho) && is.null(dbuf)) {
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
      cli::cli_abort(c("{.fn ei.estimate} did not converge. Try a different value of {.arg erho}.",
                       "i" = "Values tried: {erho}."))
    }
    cli::cli_progress_done()
    if (simulate) {
      dbuf <- ei.sim(dbuf, ndraws = ndraws, nsims = nsims)

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
    class(dbuf) <- c("ei", "eiRxC")
    cli::cli_progress_done()
  }
  dbuf
}

ei.estimate <- function(t, x, n, id, Zb = 1, Zw = 1, data = NA, erho = .5,
                        esigma = .5, ebeta = .5, ealphab = NA, ealphaw = NA,
                        truth = NA, Rfun = 2, precision = 4) {

  if (missing(t)) {
    cli::cli_abort('{.arg t} is required for {.fn ei.estimate}.')
  }
  if (missing(x)) {
    cli::cli_abort('{.arg x} is required for {.fn ei.estimate}.')
  }
  if (missing(n)) {
    cli::cli_abort('{.arg n} is required for {.fn ei.estimate}.')
  }
  if (missing(id)) {
    cli::cli_abort('{.arg id} is required for {.fn ei.estimate}.')
  }

  # Check to make sure data is not null
  if (!missing(data)) {
    if (is.character(t)) t <- data[[t]]
    if (is.character(x)) x <- data[[x]]
    if (is.character(n)) n <- data[[n]]
    if (is.character(Zb)) Zb <- data[[Zb]]
    if (is.character(Zw)) Zw <- data[[Zw]]
    if (is.character(id)) id <- data[[id]]
  }

  if (any(is.na(t)) | any(is.na(x)) | any(is.na(n))) {
    drops <- union(union(which(is.na(t)), which(is.na(x))), which(is.na(n)))
    t <- t[-drops]
    x <- x[-drops]
    n <- n[-drops]
    if (!is.na(id)) id <- id[-drops]
    cli::cli_warn('Found and removed {length(drops)} NA{?s}.')
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

  cli::cli_alert_info("Maximizing likelihood for {.arg erho} = {erho}.")

  solution <- optim(
    par = start, fn = like,
    y = t, x = x, n = n, Zb = Zb,
    Zw = Zw, numb = numb, erho = erho, esigma = esigma,
    ebeta = ebeta, ealphab = ealphab, ealphaw = ealphaw, Rfun = Rfun,
    hessian = TRUE,
    control = list(factr = 1e6, pgtol = 0.001),
    method = "L-BFGS-B"
  )
  cli::cli_process_done()

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

  class(output) <- c("ei", "ei2x2", class(output))
  output
}


#' @noRd
#' @export
print.ei <- function(x, ...) {

  cli::cli_text("An ei object with{ifelse('betab' %in% names(x), ' ', 'out ')} simulated QoIs:")

  if ('betab' %in% names(x)) {
    magg <- matrix(round(.maggs(x), digits = x$precision), nrow = 2)
    rownames(magg) <- c("Bb", "Bw")
    colnames(magg) <- c("mean", "sd")

    cli::cli_text('--- Estimates of Aggregate Quantities of Interest ---')
    cli::cat_print(magg)
  } else {
    ab <- matrix(.abounds(x), nrow = 2)
    rownames(ab) <- c("lower", "upper")
    colnames(ab) <- c("betab", "betaw")
    cli::cli_text('--- Aggregate Bounds ---')
    cli::cat_print(ab)
  }

  invisible()
}

#' Returning an element in the ei object
#' @param object An \code{ei} object from the function \code{ei}.
#' @param name The name of the element to extract from the \code{ei} object.
#' @export
values_ei <- function(object, name) {
  if (! "ei" %in% class(object)) {
    cli::cli_abort("{.arg object} is not an ei object")
  }
  if (! name %in% names(object)) {
    cli::cli_abort("{.arg name} is not an element of this ei object.")
  }
  object[[name]]
}

