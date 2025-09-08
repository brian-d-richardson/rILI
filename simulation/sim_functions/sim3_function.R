###############################################################################
###############################################################################

# rILI Simulation 3 Script

# Brian Richardson

# 2025-09-03

# The goal of this simulation is to test the empirical performance of the
# proposed randomization-based inference procedure

###############################################################################
###############################################################################

sim3.function <- function(n, tau, EY, B, R, H0A, H0Y, seed,
                          cluster.size, hall.size,
                          compute.rho.adj = T,
                          print.progress = F) {
  
  # choose STERGM parameters theta depending on H0A
  if (H0A) {
    
    # A does not depend on Z (even through Y)
    theta.p <- c(-0.5, 0, 0, 0, 1E6)  # formation parameters
    theta.m <- c(-0.5, 0, 0, 0, 1E6)  # persistence parameters
    
  } else {
    
    # edges less likely when touching treated infected person(s)
    theta.p <- c(-0.2, 0, 0, -1, 1E6)  # formation parameters
    theta.m <- c(-0.2, 0, 0, -1, 1E6)  # persistence parameters
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
  
  # choose roommate indices
  roommates <- sample(1:(n * (n-1) / 2), size = 5, replace = F)
  
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
    fixed.edges = roommates,
    print.progress = print.progress)
  
  # keep only time starting at 1
  net.df <- sim.dat$net.df[sim.dat$net.df$time > 0,]
  Z <- sim.dat$Z
  
  # sharp null p-value ------------------------------------------------------
  
  # compute observed T-stats
  T.obs <- get.T.all(net.df = net.df, Z = Z)

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
    theta.p.start = c(0, 0, 0, 0, 1E6),
    theta.m.start = c(0, 0, 0, 0, 1E6),
    offset.p = 5,
    offset.m = 5)
  
  theta.hat <- stergm.fit$theta.hat
  
  # compare true and estimated thetas
  #plot(theta.0[-c(5, 10)], theta.hat[-c(5, 10)]); abline(a = 0, b = 1)
  
  # plug-in nonsharp null p-value -------------------------------------------
  
  # compute distribution of T-stats under the nonsharp null, estimated theta
  T.nonsharp.est <- get.T.nonsharp(
    net.df = net.df,
    randomize = randomize,
    arguments = arguments,
    theta.p = theta.hat[1:5],
    theta.m = theta.hat[6:10],
    B = B,
    pb = print.progress)
  
  # plug-in nonsharp null p-value
  rho.plug <- get.rho(T.obs = T.obs, T.distr = T.nonsharp.est)

  # combine output ----------------------------------------------------------

  # combine p-values
  rho.res <- data.frame(
    tstat = names(rho.sharp),
    rho.sharp = rho.sharp,
    rho.oracle = rho.oracle,
    rho.plug = rho.plug) %>% 
    pivot_longer(cols = starts_with("rho"),
                 names_to = "ptype",
                 values_to = "value") %>% 
    mutate(n = n,
           tau = tau,
           EY = EY,
           B = B,
           H0A = H0A,
           H0Y = H0Y,
           seed = seed)
  
  # format STERGM estimation results
  theta.hat.free <- theta.hat[-c(5, 10)]
  theta.hat.sd <- sqrt(diag(-solve(stergm.fit$hessian)))
  
  theta.res <- data.frame(
    param = c(paste0("p", 1:4), paste0("m", 1:4)),
    est = theta.hat.free,
    true = c(theta.p[1:4], theta.m[1:4]),
    sd = theta.hat.sd) %>% 
    mutate(n = n,
           tau = tau,
           EY = EY,
           B = B,
           H0A = H0A,
           H0Y = H0Y,
           seed = seed)
  
  return(list(
    rho.res = rho.res,
    theta.res = theta.res))
}

