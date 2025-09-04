###############################################################################
###############################################################################

# RINI Investigatory Simulation

# Brian Richardson

# 2025-09-02

# Purpose: analyze results of simulation 2

###############################################################################
###############################################################################

# setup -------------------------------------------------------------------

rm(list = ls())
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(devtools)
library(kableExtra)
library(scales)
#setwd(dirname(getwd()))

setwd("C:/Users/brich/OneDrive - University of North Carolina at Chapel Hill/Desktop/CIRL/RINI/rILI")

# load results ------------------------------------------------------------

## load simulation results from each of 10 clusters
sim.out.list <- lapply(
  X = 0:9,
  FUN = function(clust) {
    cbind(clust,
          read.csv(paste0("simulation/sim_data/sim2/sd",
                          clust, ".csv")))
  })

## combine simulation results into 1 data frame
sim.res <- bind_rows(sim.out.list) %>% 
  mutate(pval = value,
         Stat1 = factor(substr(tstat, 2, 2),
                        levels = c(3, 2, 1),
                        labels = c("Either", "Current", "Previous")),
         Stat2 = factor(substr(tstat, 3, 3),
                        levels = 1:3,
                        labels = paste0(c("Between", "From", "To"),
                                        "\nInfected")),
         Test = gsub("rho.", "", ptype),
         Test = factor(Test,
                       levels = c("sharp", "oracle", "plug", "adj")),
         H0 = factor(paste0(as.numeric(H0A), as.numeric(H0Y)),
                     levels = c("11", "01", "00")),
         null = case_when(
           Test == "sharp" & H0 %in% c("01", "00") ~ 1,
           Test %in% c("oracle", "plug", "adj") & H0 == "00" ~ 1,
           .default = 0),
         n.tau = paste0(n, ".", tau))

# format data -------------------------------------------------------------

H0_labels <- c(
  "11" = "H[0]^{`#`}",
  "01" = "bar(H)[0]^A * \"∩\" * H[0]^Y",
  "00" = "bar(H)[0]^A * \"∩\" * bar(H)[0]^Y")

test_labels <- c(
  "sharp" = "rho[B]^{`#`}",
  "oracle" = "rho[B]^Y",
  "plug" = "hat(rho)[B]^Y",
  "adj" = "tilde(rho)[M*','*B]^Y")

n.tau_labels <- c(
  "24.5" = "n == 24 ~ ',' ~ tau == 5",
  "48.10" = "n == 48 ~ ',' ~ tau == 10")

n.reps <- n_distinct(sim.res$seed)

# check for errors --------------------------------------------------------

na.res <- sim.res %>% 
  group_by(n.tau, Test, H0, null) %>% 
  summarise(prop.na = mean(is.na(pval)))

na.res %>% 
  arrange(-prop.na)

# plot results ------------------------------------------------------------

## plots by test statistic
for (stat in levels(factor(sim.res$tstat))) {
  
  ## summarize power
  power_df <- sim.res %>%
    filter(tstat == stat) %>%
    group_by(n.tau, Test, H0, null) %>%
    summarise(power = mean(pval < 0.05, na.rm = T),
              .groups = "drop") %>% 
    mutate(label = paste0("power = ",
                          sprintf("%.1f", 100 * power), "%"),
           x = 0.5, 
           y = Inf)
  
  ## shaded rectangles for true nulls
  bg_rects <- sim.res %>%
    distinct(n.tau, Test, H0, null) %>%
    filter(null == 0) %>%
    mutate(xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf)
  
  ## histogram of p-values
  ggplot(filter(sim.res, tstat == stat),
         aes(x = pval)) +
    geom_rect(data = bg_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              fill = "gray") +  # Light gray background
    geom_histogram(bins = 20,
                   color = "black",
                   fill = "black") +
    facet_nested(Test ~ n.tau + H0,
                 labeller = labeller(
                   H0 = as_labeller(H0_labels, label_parsed),
                   Test = as_labeller(test_labels, label_parsed),
                   n.tau = as_labeller(n.tau_labels, label_parsed)),
                 scales = "free") + 
    geom_vline(xintercept = 0.05,
               color = "blue") +
    geom_text(data = power_df,
              size = 4,
              aes(x = x,
                  y = y,
                  label = label),
              vjust = 1.5,
              inherit.aes = FALSE) +  
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    labs(x = "p-value") +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_cartesian(xlim = c(0,1)) +
    ggtitle(paste0("Distribution of P-Values using Test Statistic ", stat),
            subtitle = paste0("Based on ", n.reps, " Simulations"))
  
  ## save image
  ggsave(paste0("simulation/sim_figures/sim2/sim2_hist_", stat, ".png"),
         dpi = 300, width = 10, height = 6)
  
  ## cdf of p-values
  ggplot(filter(sim.res, tstat == stat),
         aes(x = pval)) +
    geom_rect(data = bg_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              fill = "gray") +  
    stat_ecdf(geom = "step") +
    facet_nested(Test ~ n.tau + H0,
                 labeller = labeller(
                   H0 = as_labeller(H0_labels, label_parsed),
                   Test = as_labeller(test_labels, label_parsed),
                   n.tau = as_labeller(n.tau_labels, label_parsed)),
                 scales = "free") + 
    geom_abline(slope = 1,
                intercept = 0.05,
                color = "blue") +
    geom_text(data = power_df,
              size = 4,
              aes(x = x,
                  y = y,
                  label = label),
              vjust = 1.5,
              inherit.aes = FALSE) +  
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    labs(x = "p-value") +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_cartesian(xlim = c(0, 1),
                    ylim = c(0, 1)) +
    ggtitle(paste0("Empirical CDF of P-Values using Test Statistic ", stat),
            subtitle = paste0("Based on ", n.reps, " Simulations"))
  
  ## save image
  ggsave(paste0("simulation/sim_figures/sim2/sim2_ecdf_", stat, ".png"),
         dpi = 300, width = 10, height = 6)
  
  
}
