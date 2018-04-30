#Author: Phuc Nguyen
#Date: 03/18/2018
#Title: Stopping Rules Analysis (Part 3)

#import modules and packages
#========================================================================
setwd("~/Documents/Macalester/Senior/mcmc_se")
# source("utils.R")
source("metropolis.R")
source("pdfs.R")
library(ggplot2)
library(mcmcse)
library(tidyverse)
library(Matching)
library(coda)
library(cramer)
library(gridExtra)
set.seed(2018)


#declaration of helper functions
#========================================================================
#' Does pairwise comparison of chains 
#' using the KS-test to evaluate if
#' the chains have all converged to 
#' the same distribution
#' @param chains: list of at least two MCMC samples
compare_chains <- function (chains) {
  n <- length(chains)
  matrix_stats <- matrix(nrow = n, ncol = n)
  matrix_pval <- matrix(nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      result <- ks.test(chains[[i]], chains[[j]])
      matrix_stats[i, j] <- result$statistic
      matrix_pval[i, j] <- result$p.value
    }
  }
  return(list(matrix_stats, matrix_pval))
}


#' Get the maximum KS-stats and its p-value
#' @param stats: output of compare_chains function
get_max_stats <- function (stats) {
  current_ks <- stats[[1]]
  current_pval <- stats[[2]]
  ks <- max(current_ks, na.rm = TRUE)
  ks_pval <- current_pval[which(current_ks == ks)[1]]
  return(list(ks, ks_pval))
}


#' Plots and saves the max KS-stats and its p-value evalutated
#' at given intervals as MCMC samples increases
#' Native standard error of expected value of param
#' also plotted as comparision
#' @param ks: vector of max KS-stats at each interval
#' @param pval: p-value of max KS-stats
#' @param naive: native s.e. as implemented in package `mcmcse`
#' @param ...: length of each MCMC sample, inteval, number of 
#' parallel chains for naming plot
plot_stats <- function (values, values_name, stops, chain_length, interval, n_chains, target="") {
  for (i in 1:length(values)) {
    ggplot(NULL, aes(x=seq(1, chain_length, interval))) +
      geom_line(aes(y=values[[i]])) +
      geom_vline(xintercept = stops[i], linetype=3, colour='red') +
      ylab(values_name[i]) +
      xlab('chain length')
    ggsave(paste("plots/multichain",
                 Sys.Date(),
                 chain_length + interval + n_chains,
                 target,
                 i,
                 ".png", sep = "_"))
  }
}

#' Plots and saves density of MCMC chains
#' just a visual to see if they seem to 
#' have converged to the same distribution
#' @param all_chains: list of arrays of MCMC samples
#' @param inits: vector of initial values for the chains
plot_densities <- function (all_chains, inits, target="") {
  n_chains <- length(all_chains)
  chain_length <- length(all_chains[[1]])
  data <- data.frame(id = 1:chain_length)
  for (i in 1: n_chains) {
    data[as.character(i)] = all_chains[[i]]
  }
  data <- data %>% gather(key=chain, value=sample, -id)
  
  ggplot(data, aes(x=sample, colour=chain)) + geom_density()
  ggsave(paste("plots/densityOfChains",
            Sys.Date(),
            chain_length + n_chains,
            target,
            ".png", sep = "_"))
}


#' Get name of the target distribution
get_target_name <- function (target) {
  if (target(1) == df(1)) {
    target_name = "normal"
  } else if (target(1) == de(1)) {
    target_name = "exponential"
  } else if (target(1) == dbimod(1)) {
    target_name = "bimodal"
  }
  return (target_name)
}


# Analysis (main simulation)
#========================================================================
#
# Using the KS tests on multiple chains
# Target distribution is a standard normal
# ***Proposal distribution is a normal, mean 0, sd 2***

diagnose <- function (target, n_chains, chain_length, interval, init_vals, epsilon=0.05) {
  n_stops <- chain_length/interval
  target_name <- get_target_name(target)
  
  all_chains <- list()
  for (i in 1:n_chains) {
    all_chains[[i]] <- Metropolis(target, rq, 
                                  chain_length = chain_length, 
                                  y_0 = init_vals[i])
  }
  
  ks <- rep(0, n_stops)
  ks_pval <- rep(0, n_stops)
  naive_se <- rep(0, n_stops)
  gelman <- rep(0, n_stops)
  stops <- c(0,0,0)
  
  for (i in 1:n_stops) {
    current_chain_length <- interval*i
    current_chains <- list()
    current_mcmc_list <- list()
    for (j in 1:n_chains) {
      current_chains[[j]] <- all_chains[[j]][1:current_chain_length]
    }
    for (j in 1:n_chains) {
      current_mcmc_list[[j]] <- mcmc(all_chains[[j]], end = current_chain_length)
    }
    mcmc_list_obj <- mcmc.list(current_mcmc_list)
    
    stats <- compare_chains(current_chains)
    max_stats <- get_max_stats(stats)
    ks[i] <- max_stats[[1]]
    
    # Check if should stop the simulation
    if (ks[i] <= epsilon & stops[1] == 0) {
      stops[1] <- current_chain_length
    } else if (ks[i] > epsilon & stops[1] != 0) {
      stops[1] <- 0
    }
    
    ks_pval[i] <- max_stats[[2]]
    naive_se[i] <- mcse(current_chains[[1]])[[2]]
    # Check if should stop the simulation
    if (naive_se[i] <= epsilon & stops[2] == 0) {
      stops[2] <- current_chain_length
    } else if (naive_se[i] > epsilon & stops[2] != 0) {
      stops[2] <- 0
    }
    
    gelman[i] <- gelman.diag(mcmc_list_obj)$psrf[1]
    # Check if should stop the simulation
    if ((gelman[i] <= 1.005 | gelman[i] >= 0.995) & stops[3] == 0) {
      stops[3] <- current_chain_length
    } else if ((gelman[i] > 1.005 | gelman[i] < 0.995) & stops[3] != 0) {
      stops[3] <- 0
    }
  }
  
  plot_stats(list(ks, naive_se, gelman),
             c('Kolmogorov-Smirnov', 'Naive s.e.', 'Gelman-Rubin'),
             stops,
             chain_length, interval, n_chains,
             target_name)
  
  plot_densities(all_chains, init_vals, target_name)
}


main <- function () {
  # Since proposal is a normal, mean 0, sd 2,
  # the initial conditions for 4 chains are
  # the 2.5, 50, 97.5 quantile and random number
  
  target <- df
  diagnose(target = target, n_chains = 4,
         chain_length = 20000, interval = 100,
         init_vals = c(-3.919928, 0, 3.919928, 1.652525))
  
  target <- de
  diagnose(target = target, n_chains = 4,
          chain_length = 20000, interval = 100,
          init_vals = qexp(c(0.025, 0.5, 0.75, 0.975), 1))

  target <- dbimod
  diagnose(target = target, n_chains = 4,
           chain_length = 20000, interval = 100,
           init_vals = c(-0.5, 4.5, 5.5, 7))
}


# Run main()
#-------------------------------------
main()



