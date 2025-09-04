#' log-likelihood for a dyad independence STERGM
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param theta.p a numeric vector, the formation model parameters
#' @param theta.m a numeric vector, the persistence model parameters
#'
#' @return a number, the log-likelihood value
#'
#' @export
di.stergm.loglik <- function(net.df, theta.p, theta.m) {

  #theta.p <- c(0, 0); theta.m <- c(0, 0)

  # ensure data is sorted by time
  net.df <- net.df[order(net.df$time, net.df$head, net.df$tail), ]
  tau <- max(net.df$time)

  total <- 0

  for (k in 1:(tau-1)) {

    # subset data for time k and k+1
    df.k   <- subset(net.df, time == k)
    df.k1  <- subset(net.df, time == k+1)

    # compute linear predictors
    X <- as.matrix(df.k1[, grep("^X", names(df.k1))])
    eta.p  <- X %*% theta.p
    eta.m <- X %*% theta.m

    # edges at time k and k+1
    a.k  <- df.k$edge
    a.k1 <- df.k1$edge

    # Contribution for this time step
    term <- (1 - a.k) * (a.k1 * eta.p - log1p(exp(eta.p))) +
      (    a.k) * (a.k1 * eta.m - log1p(exp(eta.m)))

    total <- total + sum(term)
  }

  return(total)
}

#' score function for a dyad independence STERGM
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param theta.p a numeric vector, the formation model parameters
#' @param theta.m a numeric vector, the persistence model parameters
#'
#' @return a numeric vector, the score value
#'
#' @export
di.stergm.score <- function(net.df, theta.p, theta.m) {

  #theta.p <- c(0, 0); theta.m <- c(0, 0)

  # ensure data is sorted by time
  net.df <- net.df[order(net.df$time, net.df$head, net.df$tail), ]
  tau <- max(net.df$time)

  # Initialize score vectors
  l.p <- length(theta.p)
  l.m <- length(theta.m)
  score.p <- rep(0, l.p)
  score.m <- rep(0, l.m)

  for (k in 1:(tau-1)) {

    # subset data for time k and k+1
    df.k   <- subset(net.df, time == k)
    df.k1  <- subset(net.df, time == k+1)

    # covariate matrix
    X <- as.matrix(df.k1[, grep("^X", names(df.k))])

    # linear predictors
    eta.p <- X %*% theta.p
    eta.m <- X %*% theta.m

    # edges at time k and k+1
    a.k  <- df.k$edge
    a.k1 <- df.k1$edge

    # contributions to score
    score.p <- score.p + t(X) %*% ((1 - a.k) * (a.k1 - plogis(eta.p)))
    score.m <- score.m + t(X) %*% (    a.k  * (a.k1 - plogis(eta.m)))
  }

  # return combined score vector
  return(c(score.p, score.m))
}

#' Hessian matrix for a dyad independence STERGM
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param theta.p a numeric vector, the formation model parameters
#' @param theta.m a numeric vector, the persistence model parameters
#'
#' @return a numeric matrix, the Hessian value
#'
#' @export
di.stergm.hessian <- function(net.df, theta.p, theta.m) {

  #theta.p <- c(0, 0); theta.m <- c(0, 0)

  # Ensure data is sorted by time
  net.df <- net.df[order(net.df$time, net.df$head, net.df$tail), ]
  tau <- max(net.df$time)

  # Number of covariates
  dim.x <- length(theta.p)

  # Initialize Hessians
  H.p <- matrix(0, nrow = dim.x, ncol = dim.x)
  H.m <- matrix(0, nrow = dim.x, ncol = dim.x)

  for (k in 1:(tau-1)) {

    # subset data for time k and k+1
    df.k   <- subset(net.df, time == k)
    df.k1  <- subset(net.df, time == k+1)

    # covariate matrix
    X <- as.matrix(df.k[, grep("^X", names(df.k))])

    # linear predictors
    eta.p <- X %*% theta.p
    eta.m <- X %*% theta.m

    # edges at time k
    a.k <- df.k$edge

    # plogis values
    mu.p <- plogis(eta.p)
    mu.m <- plogis(eta.m)

    # weights for Hessian
    w.p <- (1 - a.k) * mu.p * (1 - mu.p)
    w.m <- a.k       * mu.m * (1 - mu.m)

    # contribution to Hessian: t(X * w) %*% X
    H.p <- H.p + t(sweep(X, 1, w.p, "*")) %*% X
    H.m <- H.m + t(sweep(X, 1, w.m, "*")) %*% X
  }

  # combine into block-diagonal matrix
  H <- matrix(0, nrow = 2*dim.x, ncol = 2*dim.x)
  H[1:dim.x, 1:dim.x] <- H.p
  H[(dim.x+1):(2*dim.x), (dim.x+1):(2*dim.x)] <- H.m

  return(-1 * H)
}

#' fit dyad independence STERGM
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param theta.p.start a numeric vector, starting value for formation model parameters
#' @param theta.m.start a numeric vector, starting value for the persistence model parameters
#' @param offset.p a numeric vector, the indices of offset components of theta.p
#' @param offset.m a numeric vector, the indices of offset components of theta.m
#'
#'
#' @return a numeric matrix, the Hessian value
#'
#' @importFrom rootSolve multiroot
#'
#' @export
di.stergm.fit <- function(net.df, theta.p.start, theta.m.start,
                          offset.p = NULL, offset.m = NULL) {
  
  # Indices of free parameters
  free.p.idx <- setdiff(seq_along(theta.p.start), offset.p)
  free.m.idx <- setdiff(seq_along(theta.m.start), offset.m)
  
  # Starting values for free parameters
  theta.free.start <- c(theta.p.start[free.p.idx], theta.m.start[free.m.idx])
  
  # Function to reconstruct full theta vector
  reconstruct_theta <- function(theta.free) {
    theta <- c(theta.p.start, theta.m.start)
    theta[free.p.idx] <- theta.free[seq_along(free.p.idx)]
    theta[length(theta.p.start) + free.m.idx] <- theta.free[(length(free.p.idx) + 1):length(theta.free)]
    theta
  }
  
  # Objective using only free parameters
  f.free <- function(theta.free) {
    theta <- reconstruct_theta(theta.free)
    score.full <- di.stergm.score(
      net.df = net.df,
      theta.p = head(theta, length(theta.p.start)),
      theta.m = tail(theta, length(theta.m.start))
    )
    c(score.full[free.p.idx], score.full[length(theta.p.start) + free.m.idx])
  }
  
  # Hessian using only free parameters
  jac.free <- function(theta.free) {
    theta <- reconstruct_theta(theta.free)
    hess.full <- di.stergm.hessian(
      net.df = net.df,
      theta.p = head(theta, length(theta.p.start)),
      theta.m = tail(theta, length(theta.m.start))
    )
    idx <- c(free.p.idx, length(theta.p.start) + free.m.idx)
    hess.full[idx, idx, drop = FALSE]
  }
  
  # Run root solver
  root.res <- multiroot(
    f = f.free,
    jacfunc = jac.free,
    start = theta.free.start
  )
  
  # Reconstruct full theta including offsets
  theta.hat <- reconstruct_theta(root.res$root)
  
  # Hessian matrix for covariance estimate
  hess <- jac.free(root.res$root)
  
  # Return list
  return(list(
    theta.hat = theta.hat,
    hessian = hess,
    offset.p = offset.p,
    offset.m = offset.m,
    root.res = root.res
  ))
}





di.stergm.fit.nooffset <- function(net.df, theta.p.start, theta.m.start) {

   multiroot(
    f = function(theta) {
      di.stergm.score(net.df = net.df,
                      theta.p = head(theta, length(theta.p.start)),
                      theta.m = tail(theta, length(theta.m.start))) },
    jacfunc = function(theta) {
      di.stergm.hessian(net.df = net.df,
                        theta.p = head(theta, length(theta.p.start)),
                        theta.m = tail(theta, length(theta.m.start))) },
    start = c(theta.p.start, theta.m.start))
}
