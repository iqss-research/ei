plot_tomog <- function(ei.object, title = "Tomography Plot with the Data", lci = TRUE) {
  plot_tomogd(ei.object$x, ei.object$t, ei.object$n, title, lci)
}

plot_tomogd <- function(x, t, n, title, lci = TRUE) {
  # Take out the bounds
  bounds <- bounds1(x, t, n)

  tb <- tibble::tibble(
    b_bounds = cbind(bounds[, 1], bounds[, 2]),
    w_bounds = cbind(bounds[, 4], bounds[, 3]),
    length = sqrt(abs(b_bounds[, 1] - b_bounds[, 2])^2 +
      abs(w_bounds[, 1] - w_bounds[, 2])^2),
    scale = ((length - min(length)) / (max(length) - min(length))) * 100
  )

  # Plot
  p <- tb %>%
    ggplot2::ggplot() +
    ggplot2::guides(color = "none") +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = latex2exp::TeX("$\\beta_B$"), y = latex2exp::TeX("$\\beta_W$")) +
    theme_ei()

  if (lci) {
    p <- p +
      ggplot2::geom_segment(aes(
        x = b_bounds[, 1], y = w_bounds[, 1],
        xend = b_bounds[, 2], yend = w_bounds[, 2],
        color = hcl(h = 30, c = 100, l = scale)
      )) +
      ggplot2::scale_color_manual(values = hcl(h = 30, c = 100, l = scale))
  } else {
    p <- p +
      ggplot2::geom_segment(aes(
        x = b_bounds[, 1], y = w_bounds[, 1],
        xend = b_bounds[, 2], yend = w_bounds[, 2]
      ),
      color = "yellow"
      )
  }

  return(p)
}


plot_tomogl <- function(ei.object, lci = TRUE) {
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  Zb <- ei.object$Zb
  Zw <- ei.object$Zw
  phi <- ei.object$phi
  p <- plot_tomogd(x, t, n, "Tomography Plot with ML Contours", lci = lci)
  numb <- dim(Zb)[2]
  numw <- dim(Zw)[2]
  Bb0 <- phi[1]
  Bw0 <- phi[2]
  sb0 <- phi[3]
  sw0 <- phi[4]
  rho0 <- phi[5]
  Bb0v <- phi[6:(5 + numb)]
  Bw0v <- phi[(6 + numb):length(phi)]
  vars <- repar(Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, Zw)
  bb <- vars[1:length(x)]
  bw <- vars[(length(x) + 1):(2 * length(x))]
  sb <- vars[2 * length(x) + 1]
  sw <- vars[2 * length(x) + 2]
  rho <- vars[2 * length(x) + 3]
  .tomog3 <- function(bb, bw, sb, sw, rho) {
    lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
      scale = c(sb, sw),
      centre = c(mean(bb), mean(bw)), level = .914
    ), col = "blue", lwd = 4)
    lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
      scale = c(sb, sw),
      centre = c(mean(bb), mean(bw)), level = .35
    ), col = "red", lwd = 4)
    points(mean(bb), mean(bw), col = "pink", pch = 15)
  }

  .tomog3(bb, bw, sb, sw, rho)
}

.tomog3 <- function(bb, bw, sb, sw, rho) {
  lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
    scale = c(sb, sw), centre = c(mean(bb), mean(bw)), level = .914
  ),
  col = "blue", lwd = 4
  )
  lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
    scale = c(sb, sw),
    centre = c(mean(bb), mean(bw)), level = .35
  ), col = "red", lwd = 4)
  points(mean(bb), mean(bw), col = "pink", pch = 15)
}


#
# For plot_tomog80CI()
#

tomog80CI <- function(ei.object) {
  if (!("betabs" %in% names(ei.object))) {
    message("Error: This plot function requires an ei.sim object.")
  }
  if ("betabs" %in% names(ei.object)) {
    # Only consider precincts that are heterogeneous
    ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
    x <- ei.object$x[ok]
    t <- ei.object$t[ok]
    n <- ei.object$n[ok]
    betabs <- ei.object$betabs[ok, ]
    betaws <- ei.object$betaws[ok, ]
    # .tomogd(x,t,n,"Tomography Plot with 80% CIs",lci=F)
    # Create confidence intervales
    betabcd <- apply(betabs, 1, function(x) quantile(x, probs = c(.1, .9)))
    betawcd <- apply(betaws, 1, function(x) quantile(x, probs = c(.1, .9)))
    n <- dim(betabcd)[2]
    # for(i in 1:n){
    # lines(betabcd[,i], sort(betawcd[,i],decreasing=T), col="red",
    # lwd=3)
    # }
    return(list(x = x, t = t, n = n, betabcd = betabcd, betawcd = betawcd))
  }
}

tomogd <- function(x, t, n, title, lci = T) {
  bounds <- bounds1(x, t, n)
  bbounds <- cbind(bounds[, 1], bounds[, 2])
  wbounds <- cbind(bounds[, 4], bounds[, 3])
  n <- dim(bounds)[1]
  return(bounds)
}

bounds <- function(x, t, n) {
  # set basic values
  homindx <- NULL
  tx <- NULL
  tomx <- NULL
  LbetaB <- NULL
  UbetaB <- NULL
  LbetaW <- NULL
  UbetaW <- NULL
  omx <- 1 - x
  Nb <- x * n
  Nw <- omx * n
  p <- length(x)
  homoindx <- ifelse(x == 0, 1, 0)
  homoindx <- ifelse(x == 1, 2, homoindx)


  # Heterogenous precincts
  tx <- as.matrix(t / x)
  tomx <- as.matrix(t / omx)
  tomxx <- as.matrix(tx - (omx / x))
  txx <- as.matrix(tomx - x / (1 - x))
  LbetaB <- apply(tomxx, 1, function(x) max(0, x))
  UbetaB <- apply(tx, 1, function(x) min(x, 1))
  LbetaW <- apply(txx, 1, function(x) max(0, x))
  UbetaW <- apply(tomx, 1, function(x) min(x, 1))

  # Homogenously black
  bl <- homoindx == 2
  LbetaB[bl] <- t[bl]
  UbetaB[bl] <- t[bl]
  LbetaW[bl] <- NA
  UbetaW[bl] <- NA


  # Homogenously white
  wh <- homoindx == 1
  LbetaB[wh] <- NA
  UbetaB[wh] <- NA
  LbetaW[wh] <- t[wh]
  UbetaW[wh] <- t[wh]

  return(cbind(LbetaB, UbetaB, LbetaW, UbetaW))
}
