###############################################################################
###############################################################################

# RINI Investigatory Simulation

# Brian Richardson

# 2025-09-03

# Purpose: analyze results of simulation 3

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
rho.list <- lapply(
  X = 0:9,
  FUN = function(clust) {
    cbind(clust,
          read.csv(paste0("simulation/sim_data/sim3/rho_data",
                          clust, ".csv")))
  })

theta.list <- lapply(
  X = 0:9,
  FUN = function(clust) {
    cbind(clust,
          read.csv(paste0("simulation/sim_data/sim3/theta_data",
                          clust, ".csv")))
  })

## combine simulation results into data frames
rho.res <- bind_rows(rho.list) %>% 
  mutate(pval = value,
         Stat1 = factor(substr(tstat, 2, 2),
                        levels = c(4, 3, 2, 1),
                        labels = c("Both", "Either", "Current", "Previous")),
         Stat2 = factor(substr(tstat, 3, 3),
                        levels = 1:2,
                        labels = paste0(c("Between", "From"),
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

## combine simulation results into data frames
theta.res <- bind_rows(theta.list) %>% 
  mutate(H0 = factor(paste0(as.numeric(H0A), as.numeric(H0Y)),
                     levels = c("11", "01", "00")),
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

n.reps <- n_distinct(rho.res$seed)

# check for errors --------------------------------------------------------

na.res <- rho.res %>% 
  group_by(n.tau, Test, H0, null) %>% 
  summarise(prop.na = mean(is.na(pval)))

na.res %>% 
  arrange(-prop.na)

# summarize power ---------------------------------------------------------

## power data frame
power_df <- rho.res %>%
  group_by(tstat, n.tau, Test, H0, null, Stat1, Stat2) %>%
  summarise(power = mean(pval < 0.05),
            .groups = "drop") %>% 
  mutate(label = paste0("power = ",
                        sprintf("%.1f", 100 * power), "%"),
         x.hist = 0.5, 
         y.hist = Inf,
         x.ecdf = 0.6,
         y.ecdf = 0.2)

# p-value plots -----------------------------------------------------------

## shaded rectangles for true nulls
bg_rects <- rho.res %>%
  distinct(n.tau, Test, H0, null) %>%
  filter(null == 0) %>%
  mutate(xmin = -Inf, xmax = Inf,
         ymin = -Inf, ymax = Inf)

## plots by test statistic
for (stat in levels(factor(rho.res$tstat))) {

  ## histogram of p-values
  ggplot(filter(rho.res, tstat == stat),
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
    geom_text(data = filter(power_df, tstat == stat),
              size = 3,
              aes(x = x.hist,
                  y = y.hist,
                  label = label),
              vjust = 1.5,
              inherit.aes = FALSE) +  
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank()) +
    labs(x = "p-value") +
    scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
    coord_cartesian(xlim = c(0,1)) #+
    #ggtitle(paste0("Distribution of P-Values using Test Statistic ", stat),
    #        subtitle = paste0("Based on ", n.reps, " Simulations"))
  
  ## save image
  ggsave(paste0("simulation/sim_figures/sim3/sim3_hist_", stat, ".png"),
         dpi = 300, width = 8, height = 5)
  
  ## cdf of p-values
  ggplot(filter(rho.res, tstat == stat),
         aes(x = pval)) +
    geom_rect(data = bg_rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              fill = "gray") +  
    stat_ecdf(geom = "step",
              linewidth = 0.8) +
    facet_nested(Test ~ n.tau + H0,
                 labeller = labeller(
                   H0 = as_labeller(H0_labels, label_parsed),
                   Test = as_labeller(test_labels, label_parsed),
                   n.tau = as_labeller(n.tau_labels, label_parsed)),
                 scales = "free") + 
    geom_abline(slope = 1,
                intercept = 0,
                color = "magenta",
                linetype = "dashed",
                linewidth = 0.6,
                alpha = 1) +
    geom_text(data = filter(power_df, tstat == stat),
              size = 3,
              aes(x = x.ecdf,
                  y = y.ecdf,
                  label = label),
              vjust = 1.5,
              inherit.aes = FALSE) +  
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank()) +
    labs(x = "p-value") +
    scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
    coord_cartesian(xlim = c(0, 1),
                    ylim = c(0, 1)) #+
    #ggtitle(paste0("Empirical CDF of P-Values using Test Statistic ", stat),
    #        subtitle = paste0("Based on ", n.reps, " Simulations"))
  
  ## save image
  ggsave(paste0("simulation/sim_figures/sim3/sim3_ecdf_", stat, ".png"),
         dpi = 300, width = 8, height = 5)

}


# power heatmap -----------------------------------------------------------

power_df %>% 
  mutate(
    Power = 100 * power,
    tile.label = paste0(tstat, " (",
                        sprintf("%.1f", Power), "%)")) %>% 
  ggplot(aes(y = Stat1,
             x = Stat2,
             fill = Power)) +
  geom_tile(color = "white") +
  geom_text(aes(label = tile.label,
                color = Power > 50),
            size = 3,
            show.legend = F) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  scale_fill_viridis_c(option = "C",
                       limits = c(0, 100)) +
  facet_nested(Test ~ n.tau + H0,
               labeller = labeller(
                 H0 = as_labeller(H0_labels, label_parsed),
                 Test = as_labeller(test_labels, label_parsed),
                 n.tau = as_labeller(n.tau_labels, label_parsed))) + 
  labs(
    title = "Power of Test by Statistic, Test Procedure, and Null Hypothesis",
    subtitle = paste0("Based on ", n.reps, " Simulations"),
    x = "Definition of Neighbor",
    y = "Definition of Contact",
    fill = "Power") +
  theme_bw() +
  theme(legend.position = "bottom")
  

## save image
ggsave("simulation/sim_figures/sim3/sim3_heatmap.png",
       dpi = 300, width = 12, height = 6)


# STERGM parameter plots --------------------------------------------------

## parameter labels
param.labs <- c(
  "p1" = expression(theta["edges"]^"+"),
  "p2" = expression(theta[Z]^"+"),
  "p3" = expression(theta[Y]^"+"),
  "p4" = expression(theta[ZY]^"+"),
  "m1" = expression(theta["edges"]^"-"),
  "m2" = expression(theta[Z]^"-"),
  "m3" = expression(theta[Y]^"-"),
  "m4" = expression(theta[ZY]^"-"))

## summarize sterm param results
theta.summary <- theta.res %>%
  mutate(param = factor(
    param,
    levels = c(paste0("p", 1:4), paste0("m", 1:4)))) %>% 
  group_by(H0, n.tau, param) %>%
  summarise(bias = mean(est - true),
            ESE = sd(est),
            ASE = mean(sd),
            .groups = "drop")

## plot STERGM parameter estimates
ggplot(theta.res,
       aes(x = param,
           y = est - true)) +
  geom_boxplot(outlier.size = 0.5,
               alpha = 0.6) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2,
    fill = "white",
    color = "black") +
  facet_grid(H0 ~ n.tau,
             labeller = labeller(
               H0 = as_labeller(H0_labels, label_parsed),
               n.tau = as_labeller(n.tau_labels, label_parsed))) +
  labs(
    #title = "Empirical Distribution of STERGM Parameter Estimators"
    x = "Parameter",
    y = "Empirical Bias") +
  geom_hline(yintercept = 0,
             color = "blue",
             linetype = "dashed") +
  scale_x_discrete(
    labels = param.labs) +
  theme_bw()

## save image
ggsave("simulation/sim_figures/sim3/sim3_thetadistr.png",
       dpi = 300, width = 6, height = 6)

## plot ASE vs ESE
ggplot(theta.summary,
       aes(x = ESE,
           y = ASE,
           color = param)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  geom_point(size = 3) +
  facet_grid(H0 ~ n.tau,
             labeller = labeller(
               H0 = as_labeller(H0_labels, label_parsed),
               n.tau = as_labeller(n.tau_labels, label_parsed)),
             scales = "free_x") +
  labs(
    #title = "Average Estimated Standard Errors vs Empirical standard Errors",
    x = "Empirical SE (ESE)",
    y = "Average Estimated SE (ASE)") +
  scale_color_manual(
    values = RColorBrewer::brewer.pal(8, "Set1"),  # or choose your own colors
    labels = param.labs,
    name = "Parameter") +
  theme_bw() +
  theme(legend.position = "bottom")

## save image
ggsave("simulation/sim_figures/sim3/sim3_thetaase.png",
       dpi = 300, width = 6, height = 6)
  
