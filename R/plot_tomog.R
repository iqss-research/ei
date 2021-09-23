#' @export
plot_tomog <- function() {

}


#' @export
plot_tomogl <- function() {

}

#' @import magrittr
#' @importFrom rlang .data
#' @export
plot_tomog80CI <- function(dbuf) {
  # tomog80CI
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

  res_tomog80CI <- tomog80CI(dbuf) # to replicate vignette
  bounds <- bounds(dbuf$x, dbuf$t, dbuf$n)

  bbounds <- cbind(bounds[, 1], bounds[, 2]) %>%
    as_tibble(.name_repair = ~ c("b_start", "b_end"))
  wbounds <- cbind(bounds[, 4], bounds[, 3]) %>%
    as_tibble(.name_repair = ~ c("w_start", "w_end"))
  bind_cols(
    bbounds,
    wbounds
  ) -> tomo_res_bounds

  # Visualize attempt end
  tomo_res_bounds %>%
    mutate(length = sqrt((b_end - b_start)^2 + (w_end - w_start)^2)) %>%
    mutate(inv_length = 1 / length) -> tomo_res_bounds
  # Visualize attempt end


  b <- res_tomog80CI$betabcd %>% t()
  w <- res_tomog80CI$betawcd
  w <- apply(w, 2, sort, decreasing = TRUE) %>% t()

  bind_cols(
    b %>%
      as_tibble(.name_repair = ~ c("b_start", "b_end")),
    w %>%
      as_tibble(.name_repair = ~ c("w_start", "w_end"))
  ) -> tomo_res_CI

  # Visualize attempt end
  tomo_res_CI %>%
    mutate(length = sqrt((b_end - b_start)^2 + (w_end - w_start)^2)) %>%
    mutate(inv_length = 1 / length) -> tomo_res_CI
  # Visualize attempt end

  ggplot() +
    geom_segment(
      data = tomo_res_bounds,
      aes(x = b_start, y = w_start, xend = b_end, yend = w_end, alpha = inv_length),
      show.legend = FALSE
    ) +
    # geom_segment(
    #   data = tomo_res_CI,
    #   aes(x = b_start, y = w_start, xend = b_end, yend = w_end),
    #   color = "#F8766D", show.legend = FALSE
    # ) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("betaB") +
    ylab("betaW") +
    theme_tomog() -> p
    return(p)
}


#' @export
plot_tomog95CI <- function() {

}


#' @export
plot_tomogE <- function() {

}


#' @export
plot_tomogP2 <- function() {

}
