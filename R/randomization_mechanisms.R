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
  Z.new <- rbinom(n, size = 1, prob = p)
  
  return(Z.new)
}
