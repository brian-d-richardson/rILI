#' compute a p-value given observed test statistic and null distribution
#'
#' A one-sentence description of what the function does.
#'
#' @param T.obs a numeric vector, observed values of test statistics
#' @param T.distr a numeric matrix, the null distributions of test statistics
#'
#' @return a numeric vector, p-values
#'
#' @export
get.rho <- function(T.obs, T.distr) {
  colMeans(sweep(T.distr, 2, T.obs, `<=`))
}

#' simulate the null distribution of test statistics under the sharp null
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param randomize a function, the randomization mechanism
#' @param arguments a list, arguments for the randomization mechanism
#' @param B a positive integer, the number of MC replicates
#' @param pb a logical indicator for whether to print a progress bar, default is `TRUE`
#'
#' @importFrom pbapply pbvapply
#'
#' @return a numeric vector, p-values
#'
#' @export
get.T.sharp <- function(net.df, randomize, arguments, B,
                        pb = T) {

  #randomize <- randomize.bernoulli; arguments <- list(n = n, p = 0.5)

  # select whether to show progress bar
  if (pb) {
    my.vapply <- pbvapply
  } else {
    my.vapply <- vapply
  }

  T.null <- my.vapply(
    X = 1:B,
    FUN.VALUE = numeric(8),
    FUN = function(jj) {

      # re-randomization
      Z.new <- randomize(arguments = arguments)

      # compute test stats
      get.T.all(net.df = net.df, Z = Z.new)
    })

  return(t(T.null))
}

#' simulate the null distribution of test statistics under the nonsharp null
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param randomize a function, the randomization mechanism
#' @param arguments a list, arguments for the randomization mechanism
#' @param B a positive integer, the number of MC replicates
#' @param theta.p a numeric vector, the formation model parameters
#' @param theta.m a numeric vector, the persistence model parameters
#' @param pb a logical indicator for whether to print a progress bar, default is `TRUE`
#'
#' @importFrom pbapply pbvapply
#'
#' @return a numeric vector, p-values
#'
#' @export
get.T.nonsharp <- function(net.df, randomize, arguments,
                           theta.p, theta.m, B,
                           pb = T) {

  #randomize <- randomize.bernoulli; arguments <- list(n = n, p = 0.5)
  #theta <- theta.0

  # select whether to show progress bar
  if (pb) {
    my.vapply <- pbvapply
  } else {
    my.vapply <- vapply
  }
  
  # get the mean number of edges to use for baseline Erdos-Renyi network
  #mean.edges <- mean(net.df$edge)

  # get the null distribution of T stats
  T.null <- my.vapply(
    X = 1:B,
    FUN.VALUE = numeric(8),
    FUN = function(jj) {

      # re-randomization
      Z.new <- randomize(arguments = arguments)

      # simulate new networks
      net.df.new <- net.df
      net.df.new$Zhead <- Z.new[net.df.new$head]
      net.df.new$Ztail <- Z.new[net.df.new$tail]
      net.df.new <- add.edgecov(net.df.new)
      net.df.new <- simulate.fixedcov(
        net.df = net.df.new,
        theta.p = theta.p,
        theta.m = theta.m,
        #edge0.prob = mean.edges,
        print.progress = F)

      # compute test stats
      get.T.all(net.df = net.df.new, Z = Z.new)

    })

  return(t(T.null))
}

#' get bootstrap adjusted p-value testing the nonsharp null
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param randomize a function, the randomization mechanism
#' @param arguments a list, arguments for the randomization mechanism
#' @param B a positive integer, the number of MC replicates
#' @param R a positive integer, the number of bootstrap replicates
#' @param theta.p a numeric vector, the formation model parameters
#' @param theta.m a numeric vector, the persistence model parameters
#' @param rho.plug a numeric vector, the plug-in p-values testing the nonsharp null
#' @param pb a logical indicator for whether to print a progress bar, default is `TRUE`
#'
#' @importFrom pbapply pbvapply
#'
#' @return a numeric vector, p-values
#'
#' @export
get.rho.adj <- function(net.df, Z, randomize, arguments, B, R,
                        theta.p, theta.m, rho.plug,
                        pb = F) {

  # theta.p.hat = theta.hat[1:4]; theta.m.hat = theta.hat[5:8]

  # select whether to show progress bar
  if (pb) {
    my.vapply <- pbvapply
  } else {
    my.vapply <- vapply
  }
  
  # get the mean number of edges to use for baseline Erdos-Renyi network
  mean.edges <- mean(net.df$edge)

  rho.star <- my.vapply(
    X = 1:R,
    FUN.VALUE = numeric(8),
    FUN = function(rr) {

      # re-ranomization
      Z.star <- randomize(arguments = arguments)

      # generate A*_r from theta.hat and Z*_r
      net.df.star <- net.df
      net.df.star$Zhead <- Z.star[net.df.star$head]
      net.df.star$Ztail <- Z.star[net.df.star$tail]
      net.df.star <- add.edgecov(net.df.star)
      net.df.star <- simulate.fixedcov(
        net.df = net.df.star,
        theta.p = theta.p,
        theta.m = theta.m,
        edge0.prob = mean.edges,
        print.progress = F)

      # calculate T*_r
      T.star <- get.T.all(net.df.star, Z = Z.star)

      # estimate theta using A*_r
      stergm.fit.star <- di.stergm.fit(
        net.df = net.df.star,
        theta.p.start = theta.p,
        theta.m.start = theta.m)

      theta.hat.star <- stergm.fit.star$root

      # compute distribution of T-stats under the nonsharp null, theta*_r
      T.nonsharp.star <- get.T.nonsharp(
        net.df = net.df.star,
        randomize = randomize,
        arguments = arguments,
        theta.p = head(theta.hat.star, length(theta.p)),
        theta.m = head(theta.hat.star, length(theta.m)),
        B = B,
        pb = F)

      # calculate rho*
      get.rho(T.obs = T.star, T.distr = T.nonsharp.star)

    })

  # calculate adjusted p-value (compare rho.plug to bootstrap distribution)
  rho.nonsharp.adj <- get.rho(T.obs = rho.plug, T.distr = t(rho.star))

  return(rho.nonsharp.adj)
}



