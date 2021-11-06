## Designed by Christopher T. Kenny
## Based on Cory McCartan's  betab,


new_ei_tbl <- function(data, x, t, n, phi, hessian, hessianC, psi,
                       betab, betaw, sbetab, sbetaw, betabs, betaws, resamp,
                       erho, esigma, ebta, ealphab, ealphw, numb,
                       Zb, Zw, truth, precision, id = NULL) {
  data <- reconstruct.ei_tbl(data)
  attr(data, "x") <- x
  attr(data, "t") <- t
  attr(data, "n") <- n

  #  phi
  # hessian
  # hessianC
  # psi
  # betab
  #
  #  betaw
  # sbetab
  # sbetaw
  # betabs
  # betaws
  # resamp
  #
  #  erho
  # esigma
  # ebta
  # ealphab
  # ealphw
  #
  #  numb
  # Zb
  # Zw
  # truth
  # precision
  id <- NULL

  data
}

validate_ei_tbl <- function(data) {
  if (!is.data.frame(data)) {
    stop("Not a data frame")
  }
  if (!inherits(data, "ei_tbl")) {
    stop("Not a `ei_tbl` object")
  }

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
