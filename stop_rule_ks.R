#Author: Phuc Nguyen
#Date: 02/26/2018
#Title: Stopping Rules Analysis (Part 1)

#import modules and packages
#========================================================================
source("utils.R")
library(gsubfn)
library(ggplot2)
library(mcmcse)
library(tidyverse)


#declaration of helper functions
#========================================================================

#f is the target distribution: N(0, 1^2)
#returns P(X=x)
df <- function(x) {
  dnorm(x, mean = 0, sd = 1)
}

#q is the proposal distribution: N(0, 2^2)
#returns a sample from q(x|y) 
rq <- function(x) {
  rnorm(1, mean = x, sd = 2)
}

#stopping rules analysis
#========================================================================
#
#
#Using K-S Statistics between bootstrap resamples
#------------------------------------------------------------------------
#calculate K-S stat (both methods), K-S p-value (both methods), and naive s.e. 
#for mean as chain length increases (every INTERVAL steps)

set.seed(2018)

#declare variables
TOTAL_CHAIN_LENGTH <- 2000
INTERVAL <- 100
N <- 1000
num_stops <- TOTAL_CHAIN_LENGTH/INTERVAL
ks_org_stats <- rep(0, num_stops)
ks_org_stats_mean <- rep(0, num_stops)
ks_org_pvals <- rep(0, num_stops)
ks_resamp_stats <- rep(0, num_stops)
ks_resamp_stats_mean <- rep(0, num_stops)
ks_resamp_pvals <- rep(0, num_stops)
naive_se <- rep(0, num_stops)
original_chain <- c()
ks_results_list <- list()

for (i in 1:num_stops) {
  if (i < 2) {
    y0 <- 1
  } else {
    y0 <- original_chain[INTERVAL*(i-1)]
  }
  
  #run MCMC for INTERVAL number of steps
  chain <- Metropolis(df, rq, chain_length = INTERVAL, y_0=y0)
  original_chain <- c(original_chain, chain)
  
  #bootstrap 1000 resamples
  boots <- Bootstrap(original_chain, N = N)
  
  #K-S statistics
  ks_results <- my.ks.test(boots=boots, chain=original_chain)
  ks_results_list[[i]] <- ks_results
  
  ks_org_stats[i] <- max(unlist(ks_results[['ks.original.stats']]))
  ks_org_pvals[i] <- ks_results[['ks.original.pval']][which(ks_results[['ks.original.stats']] == ks_org_stats[i])]
  ks_resamp_stats[i] <- max(unlist(ks_results[['ks.resample.stats']]))
  ks_resamp_pvals[i] <- ks_results[['ks.resample.pval']][which(ks_results[['ks.resample.stats']] == ks_resamp_stats[i])]
  
  #mean K-S statistics
  ks_org_stats_mean[i] <- mean(unlist(ks_results[['ks.original.stats']]))
  ks_resamp_stats_mean[i] <- mean(unlist(ks_results[['ks.resample.stats']]))
  
  #naive standard error
  naive_se[i] <- mcse(original_chain)[[2]]
}

#plot the relationships
chain_stops <- seq(INTERVAL, TOTAL_CHAIN_LENGTH, INTERVAL)
dat <- data.frame(chain_stops, 
                  ks_org_stats, 
                  ks_resamp_stats,
                  ks_org_stats_mean, 
                  ks_resamp_stats_mean,
                  naive_se)

dat <- dat %>%
  gather("key","value", names(dat)[names(dat) != "chain_stops"])
dat1 <- dat %>%
  filter(!(key %in% c('ks_org_stats', 'ks_resamp_stats')))
dat2 <- dat %>%
  filter(!(key %in% c('ks_org_stats_mean', 'ks_resamp_stats_mean')))

datas <- c(dat, dat1, dat2)
file_names <- c(paste("stopping_rule_analysis", TOTAL_CHAIN_LENGTH, INTERVAL, Sys.Date(), ".png", sep="_"),
                paste("stopping_rule_analysis_mean", TOTAL_CHAIN_LENGTH, INTERVAL,  Sys.Date(), ".png", sep="_"),
                paste("stopping_rule_analysis_max", TOTAL_CHAIN_LENGTH, INTERVAL,  Sys.Date(), ".png", sep="_"))

for (i in 1:length(datas)) {
  df <- datas[i]
  ggplot(df, aes(x=chain_stops, y=value, colour=key)) +
    geom_point() +
    geom_line()
  ggsave(file_names[i])
}

#save results
saveRDS(ks_results_list, file=paste0("K-S_results_", Sys.Date(), ".rds"))
save(dat, file=paste0("K-S_stats_df_", Sys.Date(), ".Rda"))
