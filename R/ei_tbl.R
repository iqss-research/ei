## Designed by Christopher T. Kenny
## Based on Cory McCartan's redist objects

# helpers ----
new_ei_tbl <- function(data, x = NULL, t = NULL, n = NULL, phi = NULL,
                       hessian = NULL, hessianC = NULL, psi = NULL,
                       betab = NULL, betaw = NULL, sbetab = NULL,
                       sbetaw = NULL, betabs = NULL, betaws = NULL, resamp = NULL,
                       erho = NULL, esigma = NULL, ebeta = NULL,
                       ealphab = NULL, ealphaw = NULL, numb = NULL,
                       Zb = NULL, Zw = NULL, truth = NULL, precision = NULL,
                       id = NULL, Rfun = NULL, covs = NULL) {
  if (missing(data)) cli::cli_abort("`data` required for {.fn new_ei_tbl}.")

  data <- reconstruct.ei_tbl(data)
  data <- add_ei_attr(
    data, x, t, n, phi, hessian, hessianC, psi,
    betab, betaw, sbetab, sbetaw, betabs, betaws, resamp,
    erho, esigma, ebeta, ealphab, ealphaw, numb,
    Zb, Zw, truth, precision, id, Rfun, covs
  )

  data
}

#' attr control
#' @keywords internal
#' @noRd
add_ei_attr <- function(data, x, t, n, phi, hessian, hessianC, psi,
                        betab, betaw, sbetab, sbetaw, betabs, betaws, resamp,
                        erho, esigma, ebeta, ealphab, ealphaw, numb,
                        Zb, Zw, truth, precision, id, Rfun, covs) {
  if (!is.null(x)) attr(data, "x") <- x
  if (!is.null(t)) attr(data, "t") <- t
  if (!is.null(n)) attr(data, "n") <- n
  if (!is.null(phi)) attr(data, "phi") <- phi
  if (!is.null(hessian)) attr(data, "hessian") <- hessian
  if (!is.null(hessianC)) attr(data, "hessianC") <- hessianC
  if (!is.null(psi)) attr(data, "psi") <- psi
  if (!is.null(betab)) attr(data, "beta_b") <- betab
  if (!is.null(betaw)) attr(data, "beta_w") <- betaw
  if (!is.null(sbetab)) attr(data, "s_beta_b") <- sbetab
  if (!is.null(sbetaw)) attr(data, "s_beta_w") <- sbetaw
  if (!is.null(betabs)) attr(data, "beta_bs") <- betabs
  if (!is.null(betaws)) attr(data, "beta_ws") <- betaws
  if (!is.null(resamp)) attr(data, "resamp") <- resamp
  if (!is.null(erho)) attr(data, "erho") <- erho
  if (!is.null(esigma)) attr(data, "esigma") <- esigma
  if (!is.null(ebeta)) attr(data, "ebeta") <- ebeta
  if (!is.null(ealphab)) attr(data, "ealphab") <- ealphab
  if (!is.null(ealphaw)) attr(data, "ealphaw") <- ealphaw
  if (!is.null(numb)) attr(data, "numb") <- numb
  if (!is.null(Zb)) attr(data, "z_b") <- Zb
  if (!is.null(Zw)) attr(data, "z_w") <- Zw
  if (!is.null(truth)) attr(data, "truth") <- truth
  if (!is.null(precision)) attr(data, "precision") <- precision
  if (!is.null(id)) attr(data, "id") <- id
  if (!is.null(id)) attr(data, "Rfun") <- Rfun
  if (!is.null(id)) attr(data, "covs") <- covs

  data
}


validate_ei_tbl <- function(data) {
  if (!is.data.frame(data)) {
    cli::cli_abort("Not a data frame.")
  }
  if (!inherits(data, "ei_tbl")) {
    cli::cli_abort("Not a {.cls ei_tbl} object.")
  }

  # check for necessary attributes:
  stopifnot(!is.null(attr(data, "x")))
  stopifnot(!is.null(attr(data, "t")))
  stopifnot(!is.null(attr(data, "n")))

  data
}



reconstruct.ei_tbl <- function(data, old) {
  classes <- c("tbl_df", "tbl", "data.frame")

  if (inherits(data, "grouped_df")) {
    classes <- c("grouped_df", classes)
  }
  if (inherits(data, "rowwise_df")) {
    classes <- c("rowwise_df", classes)
  }
  if (inherits(data, "sf")) {
    classes <- c("sf", classes)
  }

  if (!missing(old)) {
    others <- setdiff(names(attributes(old)), names(attributes(data)))
    if (length(others) > 1) {
      for (i in seq_len(length(others))) {
        attr(data, others[i]) <- attr(old, others[i])
      }
    }
  }

  class(data) <- c("ei_tbl", classes)

  data
}

#' @rdname ei_tbl
#' @param x an object to be coerced
#' @export
as_ei_tbl <- function(x) {
  reconstruct.ei_tbl(x)
}

#' Convert to `ei_tbl` objects
#' @param ei.object list-based ei object to convert to tibble-based object
#' @rdname ei_tbl
#'
#' @export
#' @return ei_tbl object
#'
#' @concept tidy
#' @examples
#' data(sample_ei)
#' form <- t ~ x
#' dbuf <- ei(form, total = "n", data = sample_ei)
#' dbuf <- ei_as_ei_tbl(dbuf)
ei_as_ei_tbl <- function(ei.object) {
  ei <- tibble(
    x = ei.object$x, t = ei.object$t, n = ei.object$n,
    beta_b = ei.object$betab, beta_w = ei.object$betaw,
    s_beta_b = ei.object$sbetab, s_beta_w = ei.object$sbetaw,
    beta_bs = ei.object$betabs, beta_ws = ei.object$betaws,
    z_b = ei.object$Zb, z_w = ei.object$Zw,
    id = ei.object$id
  )

  ei <- new_ei_tbl(
    data = ei, x = "x", t = "t", n = "n",
    phi = ei.object$phi, hessian = ei.object$hessian,
    hessianC = ei.object$hessianC, psi = ei.object$psi,
    betab = "beta_b", betaw = "beta_w", sbetab = "s_beta_b",
    sbetaw = "s_beta_w", betabs = "beta_bs", betaws = "beta_ws",
    resamp = ei.object$resamp, erho = ei.object$erho,
    esigma = ei.object$esigma, ebeta = ei.object$ebeta,
    ealphab = ei.object$ealphab, ealphaw = ei.object$ealphaw,
    numb = ei.object$numb, Zb = "z_b", Zw = "z_w",
    truth = ei.object$truth, precision = ei.object$precision,
    id = ei.object$id
  )

  validate_ei_tbl(ei)
}


# dplyr internal stuff

# 'template' is the old df
#' @method dplyr_reconstruct ei_tbl
#' @export
dplyr_reconstruct.ei_tbl <- function(data, template) {
  reconstruct.ei_tbl(data, template)
}
