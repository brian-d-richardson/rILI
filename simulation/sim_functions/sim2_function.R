###############################################################################
###############################################################################

# rILI Simulation 2 Script

# Brian Richardson

# 2025-09-02

# The goal of this simulation is to test the empirical performance of the
# proposed randomization-based inference procedure

###############################################################################
###############################################################################

sim2.function <- function(n, tau, EY, B, R, H0A, H0Y, seed,
                          cluster.size, hall.size,
                          compute.rho.adj = T,
                          print.progress = F) {
  
  # choose STERGM parameters theta depending on H0A
  if (H0A) {
    
    # A does not depend on Z (even through Y)
    theta.p <- c(-0.5, 0, 0, 0)  # formation parameters
    theta.m <- c(-0.5, 0, 0, 0)  # persistence parameters
    
  } else {
    
    # edges less likely when touching treated infected person(s)
    theta.p <- c(-0.2, 0, 0, -2)  # formation parameters
    theta.m <- c(-0.2, 0, 0, -2)  # persistence parameters
  }
  
  # combine STERGM parameters
  theta.0 <- c(theta.p, theta.m)    
  
  # set seed
  set.seed(seed)
  
  # assign clusters and halls
  clusters <- rep(seq(1, ceiling(n / cluster.size)),
                  each = cluster.size)
  halls <- rep(seq(1, ceiling(n / hall.size)),
               each = hall.size)
  
  # arguments for cluster randomization
  randomize <- randomize.cluster
  arguments <- list(clusters = clusters, halls = halls)
  
  # generate networks -------------------------------------------------------
  
  sim.dat <- simulate.newcov(
    n = n,
    tau = tau,
    theta.p = theta.p,
    theta.m = theta.m,
    randomize = randomize,
    arguments = arguments,
    H0Y = H0Y,
    EY = EY,
    print.progress = print.progress)
  
  # keep only time starting at 1
  net.df <- sim.dat$net.df[sim.dat$net.df$time > 0,]
  Z <- sim.dat$Z
  
  # sharp null p-value ------------------------------------------------------
  
  # compute observed T-stats (1-12)
  T.obs <- get.T.all(net.df = net.df, Z = Z)

  # compute distribution of T-stats under the sharp null
  T.sharp <- get.T.sharp(
    net.df = net.df,
    randomize = randomize,
    arguments = arguments,
    B = B,
    pb = print.progress)
  
  # sharp null p-value
  rho.sharp <- get.rho(T.obs = T.obs, T.distr = T.sharp)

  # oracle nonsharp null p-value --------------------------------------------
  
  # compute distribution of T-stats under the nonsharp null, true theta
  T.nonsharp.0 <- get.T.nonsharp(
    net.df = net.df,
    randomize = randomize,
    arguments = arguments,
    theta.p = theta.p,
    theta.m = theta.m,
    B = B,
    pb = print.progress)
  
  # oracle nonsharp null p-value
  rho.oracle <- get.rho(T.obs = T.obs, T.distr = T.nonsharp.0)

  # estimate STERGM parameters ----------------------------------------------
  
  # estimate theta using root solving
  stergm.fit <- di.stergm.fit(
    net.df = net.df,
    theta.p = numeric(4),
    theta.m = numeric(4))
  
  theta.hat <- stergm.fit$root
  
  # compare true and estimated thetas
  #plot(theta.0, theta.hat); abline(a = 0, b = 1)
  
  # plug-in nonsharp null p-value -------------------------------------------
  
  # compute distribution of T-stats under the nonsharp null, estimated theta
  T.nonsharp.est <- get.T.nonsharp(
    net.df = net.df,
    randomize = randomize,
    arguments = arguments,
    theta.p = theta.hat[1:4],
    theta.m = theta.hat[5:8],
    B = B,
    pb = print.progress)
  
  # plug-in nonsharp null p-value
  rho.plug <- get.rho(T.obs = T.obs, T.distr = T.nonsharp.est)

  # adjusted nonsharp null p-value ------------------------------------------
  
  # adjusted p-value
  if (compute.rho.adj) {
    rho.adj <- get.rho.adj(
      net.df = net.df,
      Z = Z,
      randomize = randomize,
      arguments = arguments,
      B = B,
      R = R,
      theta.p = theta.hat[1:4],
      theta.m = theta.hat[5:8],
      rho.plug = rho.plug,
      pb = print.progress)
  } else {
    rho.adj <- rep(NA, length(rho.oracle))
  }
  
  # combine output ----------------------------------------------------------

  # combine p-values
  rho.all <- data.frame(
    tstat = names(rho.sharp),
    rho.sharp = rho.sharp,
    rho.oracle = rho.oracle,
    rho.plug = rho.plug,
    rho.adj = rho.adj) %>% 
    pivot_longer(cols = starts_with("rho"),
                 names_to = "ptype",
                 values_to = "value") %>% 
    mutate(n = n,
           tau = tau,
           EY = EY,
           B = B,
           R = R,
           H0A = H0A,
           H0Y = H0Y,
           seed = seed)
  
  return(rho.all)
}

