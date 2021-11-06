## Designed by Christopher T. Kenny
## Based on Cory McCartan's  betab,

ei_tbl <- function() {

}


# helpers ----
new_ei_tbl <- function(data, x = NULL, t = NULL, n = NULL, phi = NULL,
                       hessian = NULL, hessianC = NULL, psi = NULL,
                       betab = NULL, betaw = NULL, sbetab = NULL,
                       sbetaw = NULL, betabs = NULL, betaws = NULL, resamp = NULL,
                       erho = NULL, esigma = NULL, ebeta = NULL,
                       ealphab = NULL, ealphw = NULL, numb = NULL,
                       Zb = NULL, Zw = NULL, truth = NULL, precision = NULL,
                       id = NULL) {
  if (missing(data)) cli::cli_abort("`data` required for `new_ei_tbl`.")

  data <- reconstruct.ei_tbl(data)
  data <- add_ei_attr(data,  x, t, n, phi, hessian, hessianC, psi,
                      betab, betaw, sbetab, sbetaw, betabs, betaws, resamp,
                      erho, esigma, ebeta, ealphab, ealphw, numb,
                      Zb, Zw, truth, precision, id)

  data
}

#' attr control
#' @keywords internal
#' @noRd
add_ei_attr <- function(data, x, t, n, phi, hessian, hessianC, psi,
                        betab, betaw, sbetab, sbetaw, betabs, betaws, resamp,
                        erho, esigma, ebeta, ealphab, ealphw, numb,
                        Zb, Zw, truth, precision, id) {

  if (!is.null(x)) attr(data, "x") <- x
  if (!is.null(t)) attr(data, "t") <- t
  if (!is.null(n)) attr(data, "n") <- n
  if (!is.null(phi)) attr(data, "phi") <- phi
  if (!is.null(hessian)) attr(data, "hessian") <- hessian
  if (!is.null(hessianC)) attr(data, "hessianC") <- hessianC
  if (!is.null(psi)) attr(data, "psi") <- psi
  if (!is.null(betab)) attr(data, "betab") <- betab
  if (!is.null(betaw)) attr(data, "betaw") <- betaw
  if (!is.null(sbetab)) attr(data, "sbetab") <- sbetab
  if (!is.null(sbetaw)) attr(data, "sbetaw") <- sbetaw
  if (!is.null(betabs)) attr(data, "betabs") <- betabs
  if (!is.null(betaws)) attr(data, "betaws") <- betaws
  if (!is.null(resamp)) attr(data, "resamp") <- resamp
  if (!is.null(erho)) attr(data, "erho") <- erho
  if (!is.null(esigma)) attr(data, "esigma") <- esigma
  if (!is.null(ebeta)) attr(data, "ebeta") <- ebeta
  if (!is.null(ealphab)) attr(data, "ealphab") <- ealphab
  if (!is.null(ealphw)) attr(data, "ealphw") <- ealphw
  if (!is.null(numb)) attr(data, "numb") <- numb
  if (!is.null(Zb)) attr(data, "Zb") <- Zb
  if (!is.null(Zw)) attr(data, "Zw") <- Zw
  if (!is.null(truth)) attr(data, "truth") <- truth
  if (!is.null(precision)) attr(data, "precision") <- precision
  if (!is.null(id)) attr(data, "id") <- id

  data
}


validate_ei_tbl <- function(data) {
  if (!is.data.frame(data)) {
    stop("Not a data frame")
  }
  if (!inherits(data, "ei_tbl")) {
    stop("Not a `ei_tbl` object")
  }

  # check for necessary attributes:
  # ex: stopifnot(!is.null(attr(data, "ndists")))

  data
}



# 'template' is the old df
#' @method dplyr_reconstruct ei_tbl
#' @export
dplyr_reconstruct.ei_tbl <- function(data, template) {
  reconstruct_ei_tbl(data, template)
}

reconstruct.ei_tbl <- function(data, old) {
  classes <- c("tbl_df", "tbl", "data.frame")

  if (inherits(data, "grouped_df")) {
    classes <- c("grouped_df", classes)
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


