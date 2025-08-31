#' get test statistics using 12 different definitions of transmission event
#'
#' A one-sentence description of what the function does.
#'
#' @param net.df a network data frame
#' @param theta.p a binary numeric vector, randomization assignments
#'
#' @return a numeric vector, the values of the 12 test statistics
#'
#' @export
get.T.all <- function(net.df, Z) {
  
  tau <- max(net.df$time)
  idx <- split(seq_len(nrow(net.df)), net.df$time)
  nT <- 12
  num_sums <- numeric(nT)
  den_sums <- numeric(nT)
  
  # add supplied randomization to network df
  net.df$Zhead <- Z[net.df$head]
  net.df$Ztail <- Z[net.df$tail]
  
  for (k in 2:tau) {
    
    df.k   <- net.df[idx[[k]], ]
    df.km1 <- net.df[idx[[k-1]], ]
    
    # dyad variables at times k and k-1
    YH.k <- df.k$Yhead
    YT.k <- df.k$Ytail
    YH.km1 <- df.km1$Yhead
    YT.km1 <- df.km1$Ytail
    
    A.k  <- df.k$edge
    A.km1 <- df.km1$edge
    A.eith <- as.numeric((A.k + A.km1) > 0)
    A.both <- as.numeric((A.k + A.km1) == 2)
    ZH <- df.k$Zhead
    ZT <- df.k$Ztail
    
    # ----- T11: contact at k-1; from infected to infected -----
    num11 <- ZH * YH.km1 * A.km1 * YT.k +
      ZT * YT.km1 * A.km1 * YH.k
    den11 <- YH.km1 * A.km1 * YT.k +
      YT.km1 * A.km1 * YH.k
    num_sums[1] <- num_sums[1] + sum(num11)
    den_sums[1] <- den_sums[1] + sum(den11)
    
    # ----- T12 contact at k-1; from infected ----
    num12 <- ZH * YH.km1 * A.km1 +
      ZT * YT.km1 * A.km1
    den12 <- YH.km1 * A.km1 +
      YT.km1 * A.km1
    num_sums[2] <- num_sums[2] + sum(num12)
    den_sums[2] <- den_sums[2] + sum(den12)
    
    # ----- T13 contact at k-1; to infected -----
    num13 <- ZH * A.km1 * YT.k +
      ZT * A.km1 * YH.k
    den13 <- A.km1 * YT.k +
      A.km1 * YH.k
    num_sums[3] <- num_sums[3] + sum(num13)
    den_sums[3] <- den_sums[3] + sum(den13)
    
    # ----- T21: contact at k; from infected to infected -----
    num21 <- ZH * YH.km1 * A.k * YT.k +
      ZT * YT.km1 * A.k * YH.k
    den21 <- YH.km1 * A.k * YT.k +
      YT.km1 * A.k * YH.k
    num_sums[4] <- num_sums[4] + sum(num21)
    den_sums[4] <- den_sums[4] + sum(den21)
    
    # ----- T22 contact at k; from infected ----
    num22 <- ZH * YH.km1 * A.k +
      ZT * YT.km1 * A.k
    den22 <- YH.km1 * A.k +
      YT.km1 * A.k
    num_sums[5] <- num_sums[5] + sum(num22)
    den_sums[5] <- den_sums[5] + sum(den22)
    
    # ----- T23 contact at k; to infected -----
    num23 <- ZH * A.k * YT.k +
      ZT * A.k * YH.k
    den23 <- A.k * YT.k +
      A.k * YH.k
    num_sums[6] <- num_sums[6] + sum(num23)
    den_sums[6] <- den_sums[6] + sum(den23)
    
    # ----- T31: contact at k or k-1; from infected to infected -----
    num31 <- ZH * YH.km1 * A.eith * YT.k +
      ZT * YT.km1 * A.eith * YH.k
    den31 <- YH.km1 * A.eith * YT.k +
      YT.km1 * A.eith * YH.k
    num_sums[7] <- num_sums[7] + sum(num31)
    den_sums[7] <- den_sums[7] + sum(den31)
    
    # ----- T32 contact at k or k-1; from infected ----
    num32 <- ZH * YH.km1 * A.eith +
      ZT * YT.km1 * A.eith
    den32 <- YH.km1 * A.eith +
      YT.km1 * A.eith
    num_sums[8] <- num_sums[8] + sum(num32)
    den_sums[8] <- den_sums[8] + sum(den32)
    
    # ----- T33 contact at k or k-1; to infected -----
    num33 <- ZH * A.eith * YT.k +
      ZT * A.eith * YH.k
    den33 <- A.eith * YT.k +
      A.eith * YH.k
    num_sums[9] <- num_sums[9] + sum(num33)
    den_sums[9] <- den_sums[9] + sum(den33)
    
    # ----- T41: contact at k and k-1; from infected to infected -----
    num41 <- ZH * YH.km1 * A.eith * YT.k +
      ZT * YT.km1 * A.eith * YH.k
    den41 <- YH.km1 * A.eith * YT.k +
      YT.km1 * A.eith * YH.k
    num_sums[10] <- num_sums[10] + sum(num41)
    den_sums[10] <- den_sums[10] + sum(den41)
    
    # ----- T42 contact at k and k-1; from infected ----
    num42 <- ZH * YH.km1 * A.eith +
      ZT * YT.km1 * A.eith
    den42 <- YH.km1 * A.eith +
      YT.km1 * A.eith
    num_sums[11] <- num_sums[11] + sum(num42)
    den_sums[11] <- den_sums[11] + sum(den42)
    
    # ----- T43 contact at k and k-1; to infected -----
    num43 <- ZH * A.eith * YT.k +
      ZT * A.eith * YH.k
    den43 <- A.eith * YT.k +
      A.eith * YH.k
    num_sums[12] <- num_sums[12] + sum(num43)
    den_sums[12] <- den_sums[12] + sum(den43)
    
  }
  
  res <- num_sums / den_sums
  names(res) <- c(
    "T11", "T12", "T13",
    "T21", "T22", "T23",
    "T31", "T32", "T33",
    "T41", "T42", "T43")
  return(res)
}