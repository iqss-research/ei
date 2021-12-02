#' Visualizing EI
#'
#' @param ei.object The output of \code{ei()}
#' @param options The list of options
#' @export
plot_bound <- function(ei.object, options = list()) {

  if ("hessian" %in% names(ei.object)) {
    # 2x2 case
    options <- plot_bound_options(options)
    p <- plot_bound_base(ei.object, options)
  } else {
    p <- plot_bound_baseRxC(ei.object, options)
  }

  return(p)
}

plot_bound_options <- function(options) {
  if (!"parameter" %in% names(options)) {
    stop("Please specify `parameter` in the option")
  }

  if (!options$parameter %in% c("betab", "betaw")) {
    stop("`options$parameter` takes `betab` or `betaw`")
  }

  return(options)
}


#' @import ggplot2
#' @importFrom rlang .data
plot_bound_base <- function(ei.object, options) {
  if (options$parameter == "betab") {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    truebb <- ei.object$truth[, 1]
    bounds <- bounds1(x, t, n)

    res <- tibble::tibble(
      x = x,
      lower = bounds[, 1],
      upper = bounds[, 2],
      true = truebb
    )

    title <- "Aggregation Bias for $\\beta_B$"
    ylab <- "True $\\beta_B$"
  } else if (options$parameter == "betaw") {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    truebw <- ei.object$truth[, 2]
    bounds <- bounds1(x, t, n)

    res <- tibble::tibble(
      x = x,
      lower = bounds[, 3],
      upper = bounds[, 4],
      true = truebw
    )

    title <- "Aggregation Bias for $\\beta_W$"
    ylab <- "True $\\beta_W$"
  } else {
    stop("Invalid argument in `options$parameter`")
  }

  fit <- lm(true ~ x, data = res)
  dat <- tibble::tibble(
    x = seq(0, 1, 0.1)
  )
  dat$y <- stats::predict(fit, newdata = data.frame(x = dat$x))
  dat$lower <- 0
  dat$upper <- 0

  p <- ggplot(res, aes(x = .data$x, y = .data$true, ymin = .data$lower, ymax = .data$upper)) +
    geom_point() +
    geom_linerange() +
    geom_line(
      data = dat, aes(x = .data$x, y = .data$y),
      linetype = "dashed"
    ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = latex2exp::TeX("$X$"), y = latex2exp::TeX(ylab)) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  return(p)
}


#
# RxC case
#

#' @import ggplot2
#' @importFrom rlang .data
plot_bound_baseRxC <- function(dbuf, options) {
  formatted <- format_bound_RxC(dbuf)
  res <- formatted$res
  names <- formatted$names

  xlab <- paste("Percent", names[1], "Vote Democrat")
  ylab <- paste("Percent", names[2], "Vote Democrat")

  p <- ggplot(res, aes(x = .data$x, y = .data$y)) +
    geom_polygon(
      aes(group = .data$id),
      fill = rgb(res$red, 0, res$blue),
      alpha = res$alpha,
      show.legend = FALSE
    ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = latex2exp::TeX(xlab), y = latex2exp::TeX(ylab)) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  return(p)
}


format_bound_RxC <- function(dbuf) {
  form <- dbuf$formula
  total <- dbuf$total
  data <- dbuf$data
  n <- nrow(data)
  covariate <- NA

  # Make the bounds
  rows <- c(all.names(form)[6:(length(all.names(form)))])
  names <- rows
  cols <- c(all.names(form)[3])
  if (sum(data[, rows][, 1] < 1.1) == length(data[, rows][, 1])) {
    data <- round(data * data[, total])
  }
  bnds <- eiPack_bounds(form, data = data, rows = rows, column = cols, threshold = 0)

  # Totals
  dv <- data[, all.names(form)[3]]

  # Assign other category
  bndsoth <- bnds$bounds[[3]]
  oth <- data[, all.names(form)[length(all.names(form))]]

  # Assign x-axis category
  bndsx <- bnds$bounds[[1]]
  xcat <- data[, all.names(form)[6]]

  # Assign y-axis category
  bndsy <- bnds$bounds[[2]]
  ycat <- data[, all.names(form)[7]]

  # Minimums & Maximums
  minx <- bndsx[, 1]
  miny <- bndsy[, 1]
  minoth <- bndsoth[, 1]
  maxx <- bndsx[, 2]
  maxy <- bndsy[, 2]
  maxoth <- bndsoth[, 2]

  #####
  # Starting point when x is at minimum
  ##
  # Holding x at its minimum, what are the bounds on y?

  # When x is at its minimum, the new dv and total are:
  newdv <- dv - (minx * xcat)
  newtot <- oth + ycat
  t <- newdv / newtot
  y <- ycat / newtot

  # The new bounds on the y category are:
  lby <- cbind(miny, (t - maxoth * oth / newtot) / (y))
  lby[, 2] <- ifelse(y == 0, 0, lby[, 2])
  lowy <- apply(lby, 1, max)
  hby <- cbind((t - minoth * oth / newtot) / y, maxy)
  highy <- apply(hby, 1, min)

  ####
  # Starting point when x is at maximum
  ##
  # Holding x at its maximum, what are the bounds on y?
  # The new bounds on x are:
  newtot <- oth + xcat
  newdv <- dv - (miny * ycat)
  x <- xcat / newtot
  t <- newdv / newtot
  lbx <- cbind(minx, (t - maxoth * oth / newtot) / x)
  lbx[, 2] <- ifelse(x == 0, 0, lbx[, 2])
  lowx <- apply(lbx, 1, max)
  hbx <- cbind((t - minoth * oth / newtot) / x, maxx)
  highx <- apply(hbx, 1, min)

  ###
  # Graph starting points
  ####

  # High starting points
  hstr <- cbind(minx, highy)
  # High ending points
  hend <- cbind(highx, miny)

  # Low starting points
  lstr <- cbind(minx, lowy)
  lend <- cbind(lowx, miny)

  # Colors for covariates
  if (!is.na(covariate)) {
    redg <- data[, covariate] / (oth - data[, covariate]) / max(data[, covariate] / (oth - data[, covariate]))
    blug <- 1 - redg
  }
  if (is.na(covariate)) {
    redg <- rep(.5, length(minx))
    blug <- rep(.5, length(minx))
  }

  # Graph labels
  xl <- paste("Percent", names[1], "Vote Democrat")
  yl <- paste("Percent", names[2], "Vote Democrat")
  mn <- paste("Tomography Plot in a 2x3 Table (", names[3], " Other Category)", sep = "")

  # Only non-NA starting pts are OK
  ok <- !is.na(hstr[, 2]) & !is.na(hend[, 1])

  # All different types of polygons
  exp1 <- hstr[ok, 2] >= maxy[ok]
  exp2 <- hend[ok, 1] >= maxx[ok]
  exp3 <- lstr[ok, 2] <= miny[ok]
  exp4 <- lend[ok, 1] <= minx[ok]

  # Subsets all the variables
  hstr <- hstr[ok, ]
  hend <- hend[ok, ]
  lstr <- lstr[ok, ]
  lend <- lend[ok, ]
  dv <- dv[ok]
  ycat <- ycat[ok]
  oth <- oth[ok]
  minoth <- minoth[ok]
  xcat <- xcat[ok]
  maxy <- maxy[ok]
  maxx <- maxx[ok]
  redg <- redg[ok]
  blug <- blug[ok]

  res <- list()

  for (i in 1:dim(hstr)[1]) {
    # 4 corner polygon
    if ((exp1[i] + exp2[i] + exp3[i] + exp4[i]) == 0) {
      xaxs <- c(hstr[i, 1], lstr[i, 1], lend[i, 1], hend[i, 1])
      yaxs <- c(hstr[i, 2], lstr[i, 2], lend[i, 2], hend[i, 2])
      side1 <- c(xaxs, xaxs[1])
      side2 <- c(yaxs, yaxs[1])
      c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 1)])
      area <- abs(c1 - c2) / 2
      if (area > 0 & !is.nan(area)) {
        alpha <- min(.5 / (area * (n)), 1)
      }
      if (area == 0 | is.nan(area)) {
        alpha <- .05
      }
      border <- alpha
    }

    # Create more corners
    if ((exp1[i] == 1) & (exp2[i]) == 0) {
      cut <- (dv[i] - (oth[i]) * minoth[i]) / xcat[i] - maxy[i] * ycat[i] / xcat[i]
      kink1x <- c(cut)
      kink1y <- c(maxy[i])
    }
    if ((exp2[i] == 1) & (exp1[i]) == 0) {
      cut <- (dv[i] - (oth[i]) * minoth[i]) / ycat[i] - maxx[i] * xcat[i] / ycat[i]
      kink1x <- c(maxx[i])
      kink1y <- c(cut)
    }
    if ((exp2[i] == 1 & exp1[i] == 1)) {
      cut <- (dv[i] - (oth[i]) * minoth[i]) / ycat[i] - maxx[i] * xcat[i] / ycat[i]
      cut2 <- (dv[i] - (oth[i]) * minoth[i]) / xcat[i] - maxy[i] * ycat[i] / xcat[i]
      kink1x <- c(maxx[i], cut2)
      kink1y <- c(cut, maxy[i])
    }
    if ((exp3[i] == 1) & (exp4[i]) == 0) {
      cut <- (dv[i] - (oth[i]) * maxoth[i]) / xcat[i] - miny[i] * ycat[i] / xcat[i]
      kink2x <- c(cut)
      kink2y <- c(miny[i])
    }
    if ((exp4[i] == 1) & (exp3[i]) == 0) {
      cut <- (dv[i] - (oth[i]) * maxoth[i]) / ycat[i] - minx[i] * xcat[i] / ycat[i]
      kink2x <- c(minx[i])
      kink2y <- c(cut)
    }
    if ((exp3[i] == 1 & exp4[i] == 1)) {
      cut <- (dv[i] - (oth[i]) * maxoth[i]) / ycat[i] - minx[i] * xcat[i] / ycat[i]
      cut2 <- (dv[i] - (oth[i]) * maxoth[i]) / xcat[i] - miny[i] * ycat[i] / xcat[i]
      kink2x <- c(minx[i], cut2)
      kink2y <- c(cut, miny[i])
    }

    # Plot 5-sided polygon
    if ((exp3[i] + exp4[i]) == 0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i]) != 0) {
      xaxs <- c(hstr[i, 1], lstr[i, 1], lend[i, 1], hend[i, 1], kink1x)
      xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
      yaxs <- c(hstr[i, 2], lstr[i, 2], lend[i, 2], hend[i, 2], kink1y)
      yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
      side1 <- c(xaxs, xaxs[1])
      side2 <- c(yaxs, yaxs[1])
      c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 1)])
      area <- abs(c1 - c2) / 2
      if (area > 0 & !is.nan(area)) {
        alpha <- min(.5 / (area * (n)), 1)
      }
      if (area == 0 | is.nan(area)) {
        alpha <- .05
      }
      border <- alpha
    }

    # Another 5-sided polygon
    if ((exp1[i] + exp2[i]) == 0 & (exp1[i] + exp2[i] + exp3[i] + exp4[i]) != 0) {
      xaxs <- c(hstr[i, 1], lstr[i, 1], kink2x, lend[i, 1], hend[i, 1])
      xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
      yaxs <- c(hstr[i, 2], lstr[i, 2], kink2y, lend[i, 2], hend[i, 2])
      yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
      side1 <- c(xaxs, xaxs[1])
      side2 <- c(yaxs, yaxs[1])
      c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 1)])
      area <- abs(c1 - c2) / 2
      if (area > 0 & !is.nan(area)) {
        alpha <- min(.5 / (area * (n)), 1)
      }
      if (area == 0 | is.nan(area)) {
        alpha <- .05
      }
      border <- alpha
    }

    # Plot 6-sided polygons
    if ((exp1[i] + exp2[i]) != 0 & (exp3[i] + exp4[i]) != 0) {
      xaxs <- c(hstr[i, 1], lstr[i, 1], kink2x, lend[i, 1], hend[i, 1], kink1x)
      xaxs <- ifelse(is.nan(xaxs), 0, xaxs)
      yaxs <- c(hstr[i, 2], lstr[i, 2], kink2y, lend[i, 2], hend[i, 2], kink1y)
      yaxs <- ifelse(is.nan(yaxs), 0, yaxs)
      side1 <- c(xaxs, xaxs[1])
      side2 <- c(yaxs, yaxs[1])
      c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 1)])
      area <- abs(c1 - c2) / 2
      if (area > 0 & !is.nan(area)) {
        alpha <- min(.5 / (area * (n)), 1)
      }
      if (area == 0 | is.nan(area)) {
        alpha <- .05
      }
      border <- alpha
    }

    res[[i]] <- tibble::tibble(
      x = xaxs,
      y = yaxs,
      id = i,
      alpha = alpha,
      red = redg[i],
      blue = blug[i]
    )
  }

  return(list(res = dplyr::bind_rows(res), name = names))
}
