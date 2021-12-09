## Designed by Christopher T. Kenny
## Based on Cory McCartan's redist objects and r-spatial's sf tibbles

#' @export
dplyr::filter

#' @export
#' @importFrom dplyr rename
rename.ei_tbl <- function(.data, ...) {
  if (!requireNamespace("tidyselect", quietly = TRUE)) {
    cli::cli_abort("{.pkg tidyselect} required for tidy {.pkg ei} function. Run {.code install.packages('tidyselect')}.")
  }

  loc <- tidyselect::eval_rename(quote(c(...)), .data)

  x_name <- attr(.data, "x")
  x_loc <- match(x_name, names(.data))
  attr(.data, "x") <- ifelse(isTRUE(x_loc %in% loc), names(loc)[x_loc], x_name)

  t_name <- attr(.data, "t")
  t_loc <- match(t_name, names(.data))
  attr(.data, "t") <- ifelse(isTRUE(t_loc %in% loc), names(loc)[t_loc], t_name)

  n_name <- attr(.data, "n")
  n_loc <- match(n_name, names(.data))
  attr(.data, "n") <- ifelse(isTRUE(n_loc %in% loc), names(loc)[n_loc], n_name)

  betab_name <- attr(.data, "beta_b")
  betab_loc <- match(betab_name, names(.data))
  attr(.data, "betab") <- ifelse(isTRUE(betab_loc %in% loc), names(loc)[betab_loc], betab_name)

  betaw_name <- attr(.data, "beta_w")
  betaw_loc <- match(betaw_name, names(.data))
  attr(.data, "betaw") <- ifelse(isTRUE(betaw_loc %in% loc), names(loc)[betaw_loc], betaw_name)

  sbetab_name <- attr(.data, "s_beta_b")
  sbetab_loc <- match(sbetab_name, names(.data))
  attr(.data, "s_beta_b") <- ifelse(isTRUE(sbetab_loc %in% loc), names(loc)[sbetab_loc], sbetab_name)

  sbetaw_name <- attr(.data, "s_beta_w")
  sbetaw_loc <- match(sbetaw_name, names(.data))
  attr(.data, "s_beta_w") <- ifelse(isTRUE(sbetaw_loc %in% loc), names(loc)[sbetaw_loc], sbetaw_name)

  betabs_name <- attr(.data, "beta_bs")
  betabs_loc <- match(betabs_name, names(.data))
  attr(.data, "beta_bs") <- ifelse(isTRUE(betabs_loc %in% loc), names(loc)[betabs_loc], betabs_name)

  betaws_name <- attr(.data, "beta_ws")
  betaws_loc <- match(betaws_name, names(.data))
  attr(.data, "betaws") <- ifelse(isTRUE(betaws_loc %in% loc), names(loc)[betaws_loc], betaws_name)

  names(.data)[loc] <- names(loc)

  .data
}

#' @export
#' @importFrom dplyr mutate
mutate.ei_tbl <- function(.data, ...) {
  reconstruct.ei_tbl(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr transmute
transmute.ei_tbl <- function(.data, ...) {
  reconstruct.ei_tbl(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr filter
filter.ei_tbl <- function(.data, ..., .preserve = FALSE) {
  reconstruct.ei_tbl(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr arrange
arrange.ei_tbl <- function(.data, ...) {
  reconstruct.ei_tbl(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr distinct
distinct.ei_tbl <- function(.data, ..., .keep_all = FALSE) {
  reconstruct.ei_tbl(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr full_join
full_join.ei_tbl <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  reconstruct.ei_tbl(NextMethod(), x)
}

#' @export
#' @importFrom dplyr inner_join
inner_join.ei_tbl <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  reconstruct.ei_tbl(NextMethod(), x)
}

#' @export
#' @importFrom dplyr left_join
left_join.ei_tbl <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  reconstruct.ei_tbl(NextMethod(), x)
}

#' @export
#' @importFrom dplyr right_join
right_join.ei_tbl <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
  reconstruct.ei_tbl(NextMethod(), x)
}

#' @export
#' @importFrom dplyr slice
slice.ei_tbl <- function(.data, ...) {
  reconstruct.ei_tbl(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr group_by
group_by.ei_tbl <- function(.data, ..., .add = FALSE) {
  reconstruct.ei_tbl(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr ungroup
ungroup.ei_tbl <- function(x, ...) {
  reconstruct.ei_tbl(NextMethod(), x)
}

#' @export
#' @importFrom dplyr rowwise
rowwise.ei_tbl <- function(x, ...) {
  reconstruct.ei_tbl(NextMethod(), x)
}

rename_helper <- function(.data, attribute, loc) {
  a_name <- attr(.data, attribute)
  a_loc <- match(a_name, names(.data))
  attr(.data, attribute) <- ifelse(isTRUE(a_loc %in% loc), names(loc)[a_loc], a_name)
  .data
}
