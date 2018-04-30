#Author: Phuc Nguyen
#Date: 03/04/2018
#Title: Stopping Rules Analysis (Part 2)

#import modules and packages
#========================================================================
source("utils.R")
source("pdfs.R")
library(gsubfn)
library(ggplot2)
library(mcmcse)
library(tidyverse)

#declaration of helper functions
#========================================================================

#wrapper function for mcmcse mcse.q to calculate s.e. at quantile q,
#returning only the s.e. estimated
mcse.q.wrapper <- function(q, x) {
  return(mcmcse::mcse.q(x=x, q=q)$se)
}

#analysis
#========================================================================
#
#
#Using quantile comparisons between different chain length
#------------------------------------------------------------------------
#For a chain of length 2000, after every x steps,
#calculate the standard error for m quantiles (using mcse.q in mcmcse)
#   to approximate how close the posterior estimates at that point is to the truth
#take the maximum of these s.e.

set.seed(2018)

#declare variables that are set by user:
TOTAL_CHAIN_LENGTH <- 4000
INTERVAL <- 50
M_quantile <- 0.001 #step size for quantiles

num_stops <- TOTAL_CHAIN_LENGTH/INTERVAL
stops <- seq(INTERVAL, TOTAL_CHAIN_LENGTH, INTERVAL)
quantiles <- seq(M_quantile, 1.0-M_quantile, M_quantile)

#target functions:
pdfs <- list(df, de, dbimod)
pdfs_names <- c('normal', 'exponential', 'bimodal') #for pdf parameters, see pdfs.R file
targets <- c("pf", "pe", "pbimod")

for (i in 1:length(pdfs)) {
  pdf <- pdfs[[i]]
  pdf_name <- pdfs_names[i]
  target <- targets[i]
  
  #run a MCMC:
  print(pdf_name)
  chain <- Metropolis(pdf, rq, chain_length = TOTAL_CHAIN_LENGTH)
  
  #calculate max s.e.
  max_se <- c()
  max_q <- c()
  naive_se <- c()
  
  #DONT KNOW IF SHOULD INCLUDE THIS OR NO
  ks <- c()
  ks_pval <- c()
  
  for (i in 1:num_stops) {
    current_chain_length <- INTERVAL*i
    current_chain <- chain[0:current_chain_length]
    se <- unlist(lapply(quantiles, mcse.q.wrapper, x=current_chain))
    max_se[i] <- max(se)
    max_q[i] <- quantiles[which(se == max_se[i])]
    naive_se[i] <- mcse(current_chain)[[2]]
    #ks test
    ks_result <- ks.test(current_chain, target)
    ks[i] <- ks_result$statistic
    ks_pval[i] <- ks_result$p.value
  }
  
  data <- data.frame(stops = stops, 
                     max_se = max_se, 
                     naive_se = naive_se)
  
  data <- data %>%
    gather("key","value", names(data)[names(data) != "stops"])
  
  #plot the results
  #plot of max se compared to naive error
  ggplot(data, aes(x=stops, y=value, colour=key)) +
    geom_line() +
    labs(x="chain length", y="s.e.", title=paste("S.e. as chain length increases,", pdf_name, "target distribution"))
  ggsave(paste0("sequantile_vs_semean_", pdf_name, Sys.Date(), ".png"))
  
  #plot of max se, naive se, compared to KS stats 
  #between posterior estimate and true target
  ggplot(NULL, aes(x=stops)) +
    geom_line(aes(y = max_se, colour = "max_se")) +
    geom_line(aes(y = naive_se, colour = "naive_se")) +
    geom_line(aes(y = ks, colour="KS stats")) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "KS stats")) +
    labs(x="chain length", y="s.e.", title=paste("S.e. as chain length increases,", pdf_name, "target distribution"))
  ggsave(paste0("sequantile_vs_semean_w_KS_", pdf_name, Sys.Date(), ".png"))
  
  #plot of max se, naive se, compared to KS stats pvalues
  #between posterior estimate and true target
  ggplot(NULL, aes(x=stops)) +
    geom_line(aes(y = max_se, colour = "max_se")) +
    geom_line(aes(y = naive_se, colour = "naive_se")) +
    geom_line(aes(y = ks_pval, colour = "KS pvals")) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "KS pvals")) +
    labs(x="chain length", y="s.e.", title=paste("S.e. as chain length increases,", pdf_name, "target distribution"))
  ggsave(paste0("sequantile_vs_semean_w_KSpvals_", pdf_name, Sys.Date(), ".png"))
  
  #plot of max se and the quantile that max comes from
  ggplot(NULL, aes(x = stops)) +
    geom_line(aes(y = max_se, colour="max se")) +
    geom_line(aes(y = max_q, colour="quantile of max se")) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "quantile of max se")) +
    labs(x="chain length", y="s.e.", title=paste("S.e. as chain length increases,", pdf_name, "target distribution"))
  ggsave(paste0("quantile_of_maxse_", pdf_name, Sys.Date(), ".png"))
  
  #plot of se for quantiles considered for chain length of 4000 steps
  ggplot(NULL, aes(x=quantiles, y=se)) +
    geom_path() +
    labs(x="quantile", y="s.e.", title=paste("S.e. at different quantiles of chain", TOTAL_CHAIN_LENGTH, "steps long,", pdf_name,"target distribution"))
  ggsave(paste0("se_quantiles", pdf_name, Sys.Date(), ".png"))
}





