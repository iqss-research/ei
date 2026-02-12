#' Visualizing EI (tomography plot for RxC)
#'
#' A tomography plot for an estimated Ecological Inference model in RxC data.
#' This function supports the 2x3 case.
#'
#' @param formula A formula of the form \code{cbind(col1,
#' col2,...)~cbind(row1,row2,...)}
#' @param data data that contains the data that corresponds to the formula
#' @param total `total' is the name of the variable in the dataset that contains the number of individuals in each unit
#'
#' @return a ggplot object
#' @concept visualization
#' @examples
#' data(RxCdata)
#' formula <- cbind(turnout, noturnout) ~ cbind(white, black, hisp)
#' plot_tomogRxC(formula, RxCdata)
#' @export
plot_tomogRxC <- function(formula, data, total = NULL) {
  p <- plot_tomogRxC_base(formula, data, total)

  return(p)
}

#' @import magrittr
#' @import ggplot2
#' @import tibble
#' @importFrom rlang .data
#' @import dplyr
plot_tomogRxC_base <- function(form, data, total) {
  refine <- 500
  formatted <- format_tomog_RxC(form, data, total, refine)

  heatmap <- formatted$contour$z %>% tibble::as_tibble(.name_repair = "minimal")
  colnames(heatmap) <- as.character(1:ncol(heatmap))
  heatmap$x <- as.character(1:nrow(heatmap))
  heatmap %>%
    tidyr::pivot_longer(-"x") %>%
    dplyr::rename(y = "name") %>%
    dplyr::relocate("x", "y") %>%
    dplyr::mutate(
      x = as.numeric(.data$x) / refine,
      y = as.numeric(.data$y) / refine,
    ) -> heatmap

  color_heat <- sort(heat.colors(refine), decreasing = TRUE)

  p <- ggplot(formatted$res, aes(x = .data$x, y = .data$y)) +
    geom_tile(
      data = heatmap,
      aes(x = .data$x, y = .data$y, fill = .data$value),
      show.legend = FALSE
    ) +
    geom_polygon(
      aes(group = .data$id),
      alpha = 0.3,
      colour = "black",
      fill = NA,
      linewidth = 0.12,
      show.legend = FALSE
    ) +
    scale_fill_gradient(
      low = color_heat[1],
      high = color_heat[refine]
    ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      x = latex2exp::TeX(formatted$names$xlab),
      y = latex2exp::TeX(formatted$names$ylab)
    ) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  attr(p, "ei_data") <- formatted

  return(p)
}

format_tomog_RxC <- function(form, data, total, refine) {
  # This function supports the 2x3 case
  #   https://github.com/cran/ei/blob/master/R/tomogRxc3d.R#L1

  # Checking the formula
  pos_cbind <- which(all.names(form) == "cbind")
  if (!all(pos_cbind == c(2, 5))) {
    cli::cli_abort("`plot_tomogRxC()` only supports the 2x3 case.")
  }

  noinfocount <- 0
  dvname <- terms.formula(form)[[2]]
  covariate <- NA

  # Make the bounds
  rows <- c(all.names(form)[6:(length(all.names(form)))])
  names <- rows
  cols <- c(all.names(form)[3])
  suppressWarnings({
    bnds <- eiPack_bounds(form, data = data, rows = rows, column = cols, threshold = 0)
  })

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

  #####
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

  #
  # Graph starting points
  #
  # High starting points
  hstr <- cbind(minx, highy)
  # High ending points
  hend <- cbind(highx, miny)

  # Low starting points
  lstr <- cbind(minx, lowy)
  lend <- cbind(lowx, miny)

  xl <- paste("Percent", names[1], dvname[2])
  yl <- paste("Percent", names[2], dvname[2])
  mn <- paste("Tomography Plot in a 2x3 Table (", names[3], " Other Category)*", sep = "")
  plot_labels <- list(xlab = xl, ylab = yl, title = mn)

  #
  # Prepare plotting Data
  #
  ok <- !is.na(hstr[, 2]) & !is.na(hend[, 1])
  exp1 <- hstr[ok, 2] >= maxy[ok]
  exp2 <- hend[ok, 1] >= maxx[ok]
  exp3 <- lstr[ok, 2] <= miny[ok]
  exp4 <- lend[ok, 1] <= minx[ok]
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

  contourx <- seq(0, 1, by = 1 / refine)
  contoury <- seq(0, 1, by = 1 / refine)
  contourz <- matrix(0, nrow = length(contourx), ncol = length(contoury))

  for (i in 1:dim(hstr)[1]) {
    if ((exp1[i] + exp2[i] + exp3[i] + exp4[i]) == 0) {
      xaxs <- c(hstr[i, 1], lstr[i, 1], lend[i, 1], hend[i, 1])
      yaxs <- c(hstr[i, 2], lstr[i, 2], lend[i, 2], hend[i, 2])
      side1 <- c(xaxs, xaxs[1])
      side2 <- c(yaxs, yaxs[1])
      c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 1)])
      area <- abs(c1 - c2) / 2
      if (area == 1) {
        noinfocount <- noinfocount + 1
      }
      for (j in 1:length(as.vector(contourx))) {
        contourz[j, ] <- contourz[j, ] + ifelse(point.in.polygon(rep(contourx[j], length(contourx)), contoury, xaxs, yaxs) == 1, 1 + 1 / (.1 + area), 0)
      }
    }
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
      if (area == 1) {
        noinfocount <- noinfocount + 1
      }
      for (j in 1:length(as.vector(contourx))) {
        contourz[j, ] <- contourz[j, ] + ifelse(point.in.polygon(rep(contourx[j], length(contourx)), contoury, xaxs, yaxs) == 1, 1 + 1 / (.1 + area), 0)
      }
    }
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
      if (area == 1) {
        noinfocount <- noinfocount + 1
      }
      for (j in 1:length(as.vector(contourx))) {
        contourz[j, ] <- contourz[j, ] + ifelse(point.in.polygon(rep(contourx[j], length(contourx)), contoury, xaxs, yaxs) == 1, 1 + 1 / (.1 + area), 0)
      }
    }
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
      if (area == 1) {
        noinfocount <- noinfocount + 1
      }
      for (j in 1:length(as.vector(contourx))) {
        contourz[j, ] <- contourz[j, ] + ifelse(point.in.polygon(rep(contourx[j], length(contourx)), contoury, xaxs, yaxs) == 1, 1 + 1 / (.1 + area), 0)
      }
    }
  }


  res_polygon <- list()
  for (i in 1:dim(hstr)[1]) {
    if ((exp1[i] + exp2[i] + exp3[i] + exp4[i]) == 0) {
      xaxs <- c(hstr[i, 1], lstr[i, 1], lend[i, 1], hend[i, 1])
      yaxs <- c(hstr[i, 2], lstr[i, 2], lend[i, 2], hend[i, 2])
      side1 <- c(xaxs, xaxs[1])
      side2 <- c(yaxs, yaxs[1])
      c1 <- sum(side1[1:(length(side1) - 1)] * side2[2:length(side2)])
      c2 <- sum(side1[2:(length(side1))] * side2[1:(length(side2) - 1)])
      area <- abs(c1 - c2) / 2
    }
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
    }
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
    }
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
    }

    res_polygon[[i]] <- tibble::tibble(
      x = xaxs,
      y = yaxs,
      id = i,
      colour = heat.colors(refine)[i]
    )
  }

  if (noinfocount > 0) {
    message(paste("There are", noinfocount, "tomography polygons with no information"))
  }

  res <- list(
    res = dplyr::bind_rows(res_polygon), names = plot_labels,
    contour = tibble::tibble(
      x = contourx[2:refine], y = contoury[2:refine], z = contourz[2:refine, 2:refine]
    )
  )

  return(res)
}
