#' Visualizing EI
#'
#' @param ei.object The output of \code{ei()}
#' @param options The list of options
#' @export
plot_truth <- function(ei.object, options = list()) {
  options <- plot_truth_options(options)

  p <- plot_truth_base(ei.object, options)

  return(p)
}

plot_truth_options <- function(options) {
  return(options)
}

aggs <- function(ei.object) {
  if (!("betabs" %in% names(ei.object))) {
    stop("This eiread function requires an ei.sim object.")
  }
  ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betab <- ei.object$betabs[ok, ]
  betaw <- ei.object$betaws[ok, ]
  omx <- 1 - x
  Nb <- n * x
  Nw <- n * omx
  Bbgg <- vector(mode = "numeric", length = dim(betab)[2])
  for (i in 1:dim(betab)[2]) {
    Bbgg[i] <- weighted.mean(betab[, i], Nb)
  }
  Bwgg <- vector(mode = "numeric", length = dim(betaw)[2])
  for (i in 1:dim(betaw)[2]) {
    Bwgg[i] <- weighted.mean(betaw[, i], Nw)
  }
  return(cbind(Bbgg, Bwgg) %>% tibble::as_tibble())
}


#' @import ggplot2
#' @importFrom rlang .data
plot_truth_base <- function(ei.object, options) {
  n <- ei.object$n
  x <- ei.object$x
  omx <- 1 - x
  truebb <- ei.object$truth[, 1]
  truebw <- ei.object$truth[, 2]
  betabs <- ei.object$betabs
  betaws <- ei.object$betaws
  betab <- ei.object$betab
  betaw <- ei.object$betaw
  truthbb <- sum(truebb * n) / sum(n)
  truthbw <- sum(truebw * n) / sum(n)
  minn <- min(x * n)
  maxn <- max(x * n)

  ag <- aggs(ei.object)
  res <- tibble::tibble(
    x = x,
    betab = betab,
    betaw = betaw,
    truebb = truebb,
    truebw = truebw,
    radius = (n * x - minn + 1) / (1 + maxn - minn)
  )


  # Density of Bb Posterior and Truth
  p1 <- ggplot(ag, aes(x = .data$Bbgg)) +
    geom_density() +
    geom_vline(xintercept = truthbb, color = "red", size = 1.1) +
    labs(x = latex2exp::TeX("$\\beta_B$"), y = latex2exp::TeX("Density")) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  # Density of Bw Posterior and Truth
  p2 <- ggplot(ag, aes(x = .data$Bwgg)) +
    geom_density() +
    geom_vline(xintercept = truthbw, color = "red", size = 1.1) +
    labs(x = latex2exp::TeX("$\\beta_W$"), y = latex2exp::TeX("Density")) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  # True betab
  ci80b <- calc_CI(ei.object, alpha = 1 - 0.8)$betabcd
  low <- mean(abs(ci80b[, 1] - betab))
  high <- mean(abs(ci80b[, 2] - betab))

  p3 <- ggplot(res, aes(x = betab, y = truebb)) +
    geom_point(shape = 21, aes(size = radius), show.legend = FALSE) +
    geom_abline(intercept = 0, slope = 1) +
    geom_abline(intercept = -low, slope = 1, linetype = "dashed") +
    geom_abline(intercept = high, slope = 1, linetype = "dashed") +
    labs(x = latex2exp::TeX("Estimated $\\beta_B$"), y = latex2exp::TeX("True $\\beta_B$")) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  # True betaw
  ci80w <- calc_CI(ei.object, alpha = 1 - 0.8)$betawcd
  low <- mean(abs(ci80w[, 1] - betaw))
  high <- mean(abs(ci80w[, 2] - betaw))

  p4 <- ggplot(res, aes(x = betaw, y = truebw)) +
    geom_point(shape = 21, aes(size = radius), show.legend = FALSE) +
    geom_abline(intercept = 0, slope = 1) +
    geom_abline(intercept = -low, slope = 1, linetype = "dashed") +
    geom_abline(intercept = high, slope = 1, linetype = "dashed") +
    labs(x = latex2exp::TeX("Estimated $\\beta_W$"), y = latex2exp::TeX("True $\\beta_W$")) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  p <- (p1 + p2) / (p3 + p4)
  return(p)
}
