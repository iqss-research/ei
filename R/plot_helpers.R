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
