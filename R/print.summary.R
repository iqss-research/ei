#' @noRd
#' @export
print.summary <- function(x, ...) {
  # Show the output
  cli::cli_h1("Summary")
  for (name in names(x)) {
    cli::cli_h2(name)
    if (length(x[[name]]) == 1) {
      cli::cli_text(x[[name]])
    } else {
      print(x[[name]])
    }
    cli::cli_text("")
  }
}

# prints summary function
# print.summary <- function(x, ...) {
#   sum.object <- x
#   dec <- sum.object$precision
#   sum.object <- sum.object[1:length(sum.object) - 1]
#   top <- list()
#   num <- 1
#   for (key in names(sum.object)) {
#     val <- sum.object[[key]]
#     if ((is.character(val) || is.numeric(val)) && length(val) < 2) {
#       top[[num]] <- sum.object[[key]]
#       names(top[[num]]) <- key
#       num <- num + 1
#     }
#   }
#   if ("Resamp" %in% names(sum.object)) {
#     cli::cat_print(cat(
#       names(top[[1]]), "=", top[[1]], ",", names(top[[2]]),
#       "=", top[[2]], ",", names(top[[3]]), "=", top[[3]], ",",
#       names(top[[4]]), "=", top[[4]], ",", names(top[[5]]), "=", top[[5]]
#     ))
#   }

#   if (!("Resamp" %in% names(sum.object))) {
#     cli::cat_print(cat(
#       names(top[[1]]), "=", top[[1]], ",", names(top[[2]]),
#       "=", top[[2]], ",", names(top[[3]]), "=", top[[3]], ",",
#       names(top[[4]]), "=", top[[4]]
#     ))
#   }
#   for (key in names(sum.object)) {
#     val <- sum.object[[key]]
#     if (!((is.character(val) || is.numeric(val)) && length(val) < 2)) {
#       cli::cat_print(key)
#       print(floor(val * 10^dec) / (10^dec))
#     }
#   }
# }
