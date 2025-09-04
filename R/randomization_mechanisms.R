#' Bernoulli randomization
#'
#' A one-sentence description of what the function does.
#'
#' @param arguments a list with `n` (sample size) and `p` (probability of treatment)
#'
#' @return a numeric vector, the randomization assignments
#'
#' @export
randomize.bernoulli <- function(arguments) {
  
  # unpack the arguments
  n <- arguments$n
  p <- arguments$p
  
  if (is.null(n)) stop("Need to provide arguments$n for Bernoulli randomization")
  if (is.null(p)) stop("Need to provide arguments$p for Bernoulli randomization")
  
  # randomize Z according to Bernoulli(p)
  Z <- rbinom(n, size = 1, prob = p)
  
  return(Z)
}


#' cluster randomization
#'
#' A one-sentence description of what the function does.
#'
#' @param arguments a list with `clusters` (cluster assignments) and `halls`
#'
#' @return a numeric vector, the randomization assignments
#'
#' @export
randomize.cluster <- function(arguments) {
  
  # extract cluster and halls from arguments
  clusters <- arguments$clusters
  halls <- arguments$halls
  
  # initialize Z vector
  Z <- numeric(length(clusters))
  
  # loop through halls
  for (h in unique(halls)) {
    
    # select clusters in hall h
    clust_in_hall <- unique(clusters[halls == h])
    
    # sample half of clusters to get Z = 1
    nh <- length(clust_in_hall)
    nh1 <- floor(nh / 2)
    if (nh %% 2 == 1 && sample(c(TRUE, FALSE), 1)) nh1 <- nh1 + 1
    z1 <- sample(clust_in_hall, nh1)
    Z[clusters %in% z1 & halls == h] <- 1
  }
  
  return(Z)
  
}
