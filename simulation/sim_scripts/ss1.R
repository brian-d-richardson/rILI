###############################################################################
###############################################################################

# rILI Simulation 1 Script

# Brian Richardson

# 2025-08-31

# The goal of this simulation is to test the empirical performance of the
# proposed randomization-based inference procedure

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

# clear workspace
rm(list = ls())

# load necessary packages
library(devtools)
library(dplyr)
library(tidyr)
library(pbapply)
library(numDeriv)
library(rootSolve)
library(Matrix)

# indicator for whether this R script is being run on the cluster
cluster.id <- as.numeric(commandArgs(TRUE))
on.cluster <- length(cluster.id) > 0
if (!on.cluster) {
  setwd("C:/Users/brich/OneDrive - University of North Carolina at Chapel Hill/Desktop/CIRL/RINI/rILI/simulation/sim_scripts")
  cluster.id <- 0
}
setwd(dirname(dirname(getwd())))

# load rILI package
load_all()

# load simulation function
source("simulation/sim_functions/sim1_function.R")

# simulation parameters ---------------------------------------------------

## baseline seed (specific to cluster)
base.seed <- 10^6 * as.integer(cluster.id)

## small sample setting
setting1 <- list(
  n = 20,        # sample size
  tau = 5,       # number of time points
  EY = 0.5,      # marginal infection probability under H0Y
  B = 200,       # inner bootstrap replicates
  R = 100,       # outer bootstrap replicates
  compute.rho.adj = T)       

## large sample setting
setting2 <- list(
  n = 80,         # sample size
  tau = 20,       # number of time points
  EY = 0.5,       # marginal infection probability under H0Y
  B = 200,        # inner bootstrap replicates
  R = NA,         # outer bootstrap replicates
  compute.rho.adj = F)

## varied parameters
H0Y <- c(T, F)
H0A <- c(T, F)
#H0Y <- T; H0A <- T
#H0Y <- F; H0A <- F

## number of simulation replicates
n.rep <- 1

## test run one simulation
if (FALSE) {
  setting <- setting2
  sim1.function(
    n = setting$n,
    tau = setting$tau,
    EY = setting$EY,
    B = setting$B,
    R = setting$R,
    compute.rho.adj = F,
    print.progress = T,
    H0A = F,
    H0Y = T,
    seed = sample(1E6, 1))
}

# run simulations ---------------------------------------------------------

## simulation 
sim.in <- expand.grid(
  setting = c(1, 2),
  H0Y = c(T, F),
  H0A = c(T, F),
  sim.id = 1:n.rep + base.seed) %>% 
  filter(!(H0Y == F & H0A == T))

## run simulations
sim.out <- pblapply(
  X = seq_len(nrow(sim.in)),
  FUN = function(ii) {
    
    setting <- get(paste0("setting", sim.in$setting[ii]))
    
    sim1.function(
      n = setting$n,
      tau = setting$tau,
      EY = setting$EY,
      B = setting$B,
      R = setting$R,
      compute.rho.adj = setting$compute.rho.adj,
      H0A   = sim.in$H0A[ii],
      H0Y   = sim.in$H0Y[ii],
      seed  = sim.in$sim.id[ii],
      print.progress = T)
    
    }) %>% 
  bind_rows()

## save results
write.csv(sim.out, row.names = F,
          paste0("simulation/sim_data/sim1/sd",
                 as.integer(cluster.id), ".csv"))

