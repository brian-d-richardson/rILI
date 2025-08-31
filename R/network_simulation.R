#' initialize a network data frame
#'
#' A one-sentence description of what the function does.
#'
#' @param n a positive integer, the sample size
#' @param tau a positive integer, the number of time points
#'
#' @return a data frame
#'
#' @export
initialize.net.df <- function(n, tau) {

  # generate all unique dyads (tail < head)
  dyads <- t(combn(1:n, 2))  # each row: c(tail, head)

  # Repeat dyads for all time points
  net.df <- data.frame(
    time = rep(0:tau, each = n * (n-1) / 2),
    head = rep(dyads[,2], times = tau + 1),
    tail = rep(dyads[,1], times = tau + 1),
    edge = 0
  )

  return(net.df)

}

#' add edge-level covariates to a network data frame
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#'
#' @return a data frame
#'
#' @export
add.edgecov <- function(net.df) {

  net.df$X1 <- 1                               # intercept
  net.df$X2 <- net.df$Zhead + net.df$Ztail     # number treated
  net.df$X3 <- net.df$Yhead + net.df$Ytail     # number infected
  net.df$X4 <- (net.df$Zhead * net.df$Yhead) + # number treated and infected
    (net.df$Ztail * net.df$Ytail)

  return(net.df)
}


#' get infection probability given number of infected neighbors
#'
#' A one-sentence description of what the function does.
#'
#' @param s a numeric vector, counts of infected neighbors
#' @param H0Y a logical indicator for whether H_0^Y is true
#' @param EY a number in [0,1], the (constant) probability of infection if `H0Y = TRUE`
#' @param g.shift a number, the shift parameter in `g` if `H0Y = FALSE`
#' @param g.scale a number, the scale parameter in `g` if `H0Y = FALSE`
#'
#' @return a numeric vector, the infection probabilities
#'
#' @export
g <- function(s, H0Y, EY = NULL, g.shift = 8.5, g.scale = 0.8) {

  # under H0^Y, infection probs don't depend on # infected neighbors
  if (H0Y) {

    return(rep(EY, length(s)))

    # otherwise, infection probs increase with # infected neighbors
  } else {

    ss <- plogis(s - g.shift) * (g.scale) + (1 - g.scale)

    return(ss)
  }
}

#' count the number of infected neighbors for each individual at time k
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param k a positive integer, the time point at which to count infected neighbors
#'
#' @return a numeric vector, the counts of infected neighbors
#'
#' @export
count.infected.neighbors <- function(net.df, k) {

  # restrict to time k
  df.k <- net.df[net.df$time == k, ]
  n <- max(df.k$head, df.k$tail)

  # build infection status per node from dyad-level indicators
  Y.k <- numeric(n)
  infected.idx <- unique(c(
    df.k$head[df.k$Yhead == 1],
    df.k$tail[df.k$Ytail == 1]))
  Y.k[infected.idx] <- 1

  # sparse adjacency (undirected edges)
  ed <- which(df.k$edge == 1L)
  i <- df.k$head[ed]; j <- df.k$tail[ed]
  A.k <- Matrix::sparseMatrix(i = c(i, j), j = c(j, i), x = 1L, dims = c(n, n))

  # infected-neighbor counts
  s.k <- as.vector(A.k %*% Y.k)
  return(s.k)

}


#' simulate data according to a dyad independence STERGM with prespecified covariates
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df an initialized network data frame
#' @param theta.p a numeric vector, the formation model parameters
#' @param theta.m a numeric vector, the persistence model parameters
#' @param print.progress a logical indicator for whether to print progress, default is `TRUE`
#'
#' @return a network data frame
#'
#' @export
simulate.fixedcov <- function(net.df, theta.p, theta.m, print.progress = TRUE) {

  # Ensure data is sorted by time
  net.df <- net.df[order(net.df$time, net.df$head, net.df$tail), ]

  # covariate columns
  Xcols <- grep("^X", names(net.df))

  tau <- max(net.df$time)
  edges <- numeric(nrow(net.df))

  # ----- Time 0: Erdos-Renyi -----
  idx.0 <- which(net.df$time == 1)
  A.k <- rbinom(length(idx.0), 1, 0.5)
  edges[idx.0] <- A.k
  A.k1 <- A.k

  if (print.progress) {
    message(paste0("time 0: total edges = ", sum(A.k)))
  }

  # ----- Times 1..tau: STERGM -----
  for (k in 0:(tau-1)) {

    # rows for time k+1
    idx.k1 <- which(net.df$time == k+1)
    df.k1  <- net.df[idx.k1, ]
    X.k1   <- as.matrix(df.k1[, Xcols])

    # linear predictors
    eta.p <- X.k1 %*% theta.p  # formation logits
    eta.m <- X.k1 %*% theta.m  # persistence logits

    # edge probabilities conditional on previous edges
    probs <- ifelse(A.k == 0, plogis(eta.p), plogis(eta.m))

    # draw edges
    A.k1 <- rbinom(length(probs), 1, probs)
    edges[idx.k1] <- A.k1

    # advance
    A.k <- A.k1

    if (print.progress) {
      message(paste0("time ", k+1, ": total edges = ", sum(A.k1)))
    }
  }

  net.df$edge <- edges
  return(net.df)
}

#' simulate data according to a dyad independence STERGM with new randomizations and infection statuses
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df an initialized network data frame
#' @param theta.p a numeric vector, the formation model parameters
#' @param theta.m a numeric vector, the persistence model parameters
#' @param randomize a function, the randomization mechanism
#' @param arguments a list, arguments for the randomization mechanism
#' @param H0Y a logical indicator for whether H_0^Y is true
#' @param EY a number in [0,1], the (constant) probability of infection if `H0Y = TRUE`
#' @param g.shift a number, the shift parameter in `g` if `H0Y = FALSE`
#' @param g.scale a number, the scale parameter in `g` if `H0Y = FALSE`
#' @param print.progress a logical indicator for whether to print progress, default is `TRUE`
#'
#' @return a list with randomization assignments `Z` and a network data frame `net.df`
#'
#' @export
simulate.newcov <- function(n, tau,
                          theta.p, theta.m,
                          randomize, arguments,
                          H0Y,
                          EY = NULL,
                          g.shift = NULL,
                          g.scale = NULL,
                          print.progress = TRUE) {

  net.df <- initialize.net.df(n = n, tau = tau)

  # unit-level covariates
  Z <- randomize(arguments = arguments)
  net.df$Zhead <- Z[net.df$head]
  net.df$Ztail <- Z[net.df$tail]

  # initialize edge and infection columns
  net.df$Yhead <- 0
  net.df$Ytail <- 0
  net.df$X4 <- net.df$X3 <- net.df$X2 <- net.df$X1 <- 0
  Xcols <- paste0("X", 1:4)

  # initial time 0 edges (Erdos-Renyi)
  idx.0 <- which(net.df$time == 0)
  A.k <- rbinom(length(idx.0), 1, 0.5)
  net.df$edge[idx.0] <- A.k

  # initial infections at time 0
  Y.k <- rbinom(n, 1, EY)
  net.df$Yhead[idx.0] <- Y.k[net.df$head[idx.0]]
  net.df$Ytail[idx.0] <- Y.k[net.df$tail[idx.0]]

  # compute edge-level covariates at time 0
  net.df[idx.0, ] <- add.edgecov(net.df[idx.0, ])

  # iterate over time
  for (k in 0:(tau-1)) {

    # ----- Step 1: Identify rows for time k+1 -----
    idx.k1 <- which(net.df$time == k + 1)

    # ----- Step 2: Compute infection statuses Y^{k+1} -----
    s.k <- count.infected.neighbors(net.df, k = k)     # count infected nbrs
    g.k <- g(s = s.k, H0Y = H0Y,                       # transformation to [0,1]
             EY = EY,
             g.shift = g.shift, g.scale = g.scale)
    Y.k1 <- rbinom(n = n, size = 1, prob = g.k)        # draw new infections
    net.df$Yhead[idx.k1] <- Y.k1[net.df$head[idx.k1]]  # update head infections
    net.df$Ytail[idx.k1] <- Y.k1[net.df$tail[idx.k1]]  # update tail infections

    # ----- Step 3: Update edge-level covariates X1-X4 at time k+1 -----
    net.df[idx.k1, ] <- add.edgecov(net.df[idx.k1, ])
    X.k1 <- as.matrix(net.df[idx.k1, Xcols])

    # ----- Step 4: Compute linear predictors -----
    eta.p <- X.k1 %*% theta.p   # formation logits for edges not present at k
    eta.m <- X.k1 %*% theta.m   # persistence logits for edges present at k

    # ----- Step 5: Compute probabilities for edges A^{k+1} -----
    prob.k1 <- ifelse(A.k == 0, plogis(eta.p), plogis(eta.m))

    # ----- Step 6: Draw edges A^{k+1} -----
    A.k1 <- rbinom(length(prob.k1), 1, prob.k1)
    net.df$edge[idx.k1] <- A.k1

    # ----- Step 7: Advance time -----
    A.k <- A.k1

    if (print.progress) {
      message(paste0("k = ", k + 1,
                     "; total infections = ", sum(Y.k1),
                     "; total edges = ", sum(A.k1)))
    }
  }

  return(list(Z = Z,
              net.df = net.df))
}



