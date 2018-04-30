#Author: Phuc Nguyen
#Date: 02/09/2018
#Title: Preliminary Simulation Study for Bootstraps of MCMC samples


#import modules and packages
#========================================================================
source("utils.R")
library(gsubfn)
library(ggplot2)


#helper functions declaration
#NOTE: get.stats and ks.resamp have been added to utils.R
#      under my.ks.test function as of 2/26/2018
#================================================================

get.stats <- function(ks_list) {
  #' @param ks_list list of result from ks.test
  #' returns list of statistics and p-values
  
  n <- length(ks_list)
  ks_dist <- c()
  ks_pval <- c()
  for (i in 1:n) {
    ks_dist[i] <- ks_list[[i]]$statistic
    ks_pval[i] <- ks_list[[i]]$p.value
  }
  return(list(dist=ks_dist, pval=ks_pval))
}


ks.resamp <- function(i, boots) {
  #' @param i index of reference resample
  #' @param boots data frame of bootstrap resamples
  #' return a list of ks.test results
  
  n <- dim(boots)[2]
  ref <- boots[,i]
  if (i == (n-1)) {
    result <- list(ks.test(boots[,i], boots[,n]))
  } else {
    new_boots <- boots[,(i+1):n]
    result <- apply(new_boots, 2, ks.test, ref)
  }
  return(result)
}


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

#measure of distance
#================================================================
#Using Kolmogorov-Smirnov test statistic to measure distance
#between the boostrap posteriors
#       a. comparing the distance between each boostrap resamples
#          and the original chain.
#       b. comparing the distance between each pair bootstrap
#          resamples.
#find the biggest difference and standard deviation in each case

set.seed(2018)

max_KS <- list()
sd_KS <- list()
pvals <- list()
lens <- c(100, 1000)
N <- 1000

for (l in lens) {
  #Metropolis MCMC
  chain <- Metropolis(df, rq, chain_length = l)
  
  #Take 1000 bootstrap resamples
  boots <- Bootstrap(chain, N = N)
  
  #Calculate KS test statistics between resamples and original
  KS_org <- apply(boots, 2, ks.test, chain)
  list[KS_org_dist, KS_org_pval] <- get.stats(KS_org)
  
  #Calculate KS test statistics between resamples
  i <- 1:(N-1)
  KS_resample <- unlist(lapply(i, ks.resamp, boots), recursive = FALSE)
  list[KS_resample_dist, KS_resample_pval] <- get.stats(KS_resample)
  
  #Get maximum and standard deviation of KS statistics
  name <- paste0('original',l) 
  max_KS[[name]] <- max(KS_org_dist)
  sd_KS[[name]] <- sd(KS_org_dist)
  pvals[[name]] <- KS_org_pval[which(KS_org_pval == max_KS[[name]])]
  
  #Get maximum and standard deviation of KS statistics
  name <- paste0('resample',l)
  max_KS[[name]] <- max(KS_resample_dist)
  sd_KS[[name]] <- sd(KS_resample_dist)
  pvals[[name]] <- KS_org_pval[which(KS_org_pval == max_KS[[name]])]
}
