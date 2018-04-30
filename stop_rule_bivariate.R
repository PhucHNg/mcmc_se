#' Author: Phuc Nguyen
#' Date: 04/30/2018
#' Description: Test convergence of a bivariate normal with no covariance between
#'              the dimensions using the following statistics:
#'              Multivariate Smirnov: TODO: cite paper
#'              Multivariate Wald-Wolfewitz: TODO: cite paper as in package `GSAR`
#'              Multivariate Cramer: as in package `cramer`
#'              Multivariate Gelman-Rubin: as implemented in package `coda`
#'              Naive error for each dimension: as implemented in package `mcmcse`
#'              
#' Note: the tests from `GSAR` don't seem to work right, maybe find another implementation
#'      of the algorithms


# Load packages
#====================================================================================
setwd("~/Documents/Macalester/Senior/mcmc_se")
library(cramer)
library(GSAR)
library(coda)
library(mcmcse)
library(ggplot2)
library(MASS)
source("bivariate_normal.R")
source("multivariate_metropolis.R")

 
# Helper functions
#==============================================================================================
# Batch mean of chain to speed up Cramer, Smirnov, and Wald statistics
# return a 
bm <- function(x, batch_size=50, min_b=200, max_b=500) {
  if (dim(x)[1] < max_b) {
    return(t(x))
  }
  if (dim(x)[1] >= max_b) {
    batch_size <- floor(dim(x)[1]/min_b)
    splited_x <- matsplitter(x, batch_size, 2)
  }
  if ((dim(x)[1]/batch_size) > max_b) {
    batch_size <- floor(dim(x)[1]/max_b)
    splited_x <- matsplitter(x, batch_size, 2)
  }
  new_x <- sapply(splited_x, colMeans)
  return(new_x)
}

#https://stackoverflow.com/questions/24299171/function-to-split-a-matrix-into-sub-matrices-in-r
# ISSUE: cannot have a batch size of 1 ...
matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  result <- list()
  for (i in 1:N){
    result[[i]] <- cv[1:r,1:c,i]
  }
  return(result)
} 


# Define Target and Proposal
#====================================================================================
# Target 1: chain with two dimensions totally independent
#     mean (0, 0)
#     covariance matrix(c(1, 0, 0, 2), 2, 2)
mean_1 <- c(0, 0)
covariance_1 <- matrix(c(1, 0, 0, 2), 2, 2)
target_1 <- target_wrapper(mean = mean_1, covariance = covariance_1)


# Target 2: chain with correlated dimensions
#     mean (0, 0)
#     covariance matrix(c(1, 0.6, 0.6, 2))
# Note: not used in this simulation
mean_2 <- c(0, 0)
covariance_2 <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
target_2 <- target_wrapper(mean = mean_2, covariance = covariance_2)


# Proposal is a bivariate with covariance diag(4, 4) (each has sd = 2)
proposal <- proposal_wrapper(covariance = diag(c(4, 4)))


# Run MCMC chains
#================================================================================
# using only the nice bivariate normal
# returns a n_max by 2 matrix
n_max <- 15000
n_chains <- 2
chains <- list()
inits <- list(c(1,1), c(-1,2))
for (i in 1:n_chains) {
  chains[[i]] <- Metropolis_m(target_1, proposal, n_max, inits[[i]])
}


# Calculate diagnostics (main simulation)
#==============================================================================
#   Multivariate Cramer Test
#   Multivariate Smirnov using Minimum Spanning Tree
#   Multivariate Wald-Wolfewitz using MST
#   Multivariate Gelman-Rubin diagnostic from `coda`
#   Univariate naive Batch Mean Standard Error
interval <- 300
burn_in <- interval*1
n_checks <- (n_max-burn_in)/interval
cramer <- list()
cramer_stats <- rep(0, n_checks)
cramer_pval <- rep(0, n_checks)
cramer_result <- rep(0, n_checks)
smirnov <- rep(0, n_checks)
smirnov2 <- rep(0, n_checks)
wwolf <- rep(0, n_checks)
gelman <- rep(0, n_checks)
naive_dim1 <- rep(0, n_checks)
naive_dim2 <- rep(0, n_checks)

for (i in 1: n_checks) {
  # get current chain
  current_chain_length <- interval*i + burn_in
  print(paste("current length:", current_chain_length))
  current_chains <- list()
  for (j in 1:n_chains) {
    current_chains[[j]] <- chains[[j]][burn_in: current_chain_length, ]
  }
  x <- current_chains[[1]]
  y <- current_chains[[2]] 
  small_x <- bm(x) # use batch means to reduce current MCMC sample to 600 samples maximum
                  # can help with i.i.d assumption of Cramer and computational time (?)
  small_y <- bm(y)
  
  # calculate pairwise statistics
  # Cramer Test
  cur_cramer <- cramer::cramer.test(small_x, small_y, sim = "eigenvalue")
  cramer_stats[i] <- cur_cramer$statistic
  cramer_pval[i] <- cur_cramer$p.value
  cramer_result[i] <- cur_cramer$result
  print("calculated Cramer statistics")
  
  # Smirnov and Wild-Wolfewitz test
  # xyT <- cbind(t(x), t(y)) # takes way too long to calculate statistics with raw samples
  # groups <- c(rep(1, dim(x)[1]), rep(2, dim(y)[1]))

  xyT <-  cbind(small_x, small_y)
  groups <- c(rep(1, dim(small_x)[2]), rep(2, dim(small_y)[2]))
  smirnov[i] <- GSAR::RKStest(object = xyT, group = groups)
  print("calcualted Smirnov statistics")
  
  wwolf[i] <-GSAR::WWtest(object = xyT, group = groups)
  print("calculated WWolfewitz statistics")
  
  # Gelman-Rubin diagnostic
  current_mcmc_list <- list()
  for (j in 1:n_chains) {
    current_mcmc_list[[j]] <- mcmc(chains[[j]], end = current_chain_length)
  }
  mcmc_obj <- mcmc.list(current_mcmc_list)
  gelman[i] <- gelman.diag(mcmc_obj)$mpsrf
  print("calculated Gelman statistics")

  # Naive standard error
  naive_dim1[i] <- mcmcse::mcse(current_chains[[1]][,1])[[2]]
  naive_dim2[i] <- mcmcse::mcse(current_chains[[1]][,2])[[2]]
  print("calculated Naive statiistics")
  
  # save output incrementally to a file
  if ((i %% 5) == 0) {
    data <- (cbind(1:n_checks, cramer_stats, cramer_pval, cramer_result, smirnov, gelman, naive_dim1, naive_dim2))
    write.matrix(data, file="bivariate_output.txt")
  }
}


# Plot statistics
#==============================================================================
x_axis <- seq(burn_in, n_max-1, interval)

# Cramer test
ggplot(NULL, aes(x = x_axis)) +
  geom_line(aes(y=cramer_stats), linetype=1, colour='red') +
  geom_line(aes(y=cramer_pval), linetype=2) +
  geom_line(aes(y=cramer_result), linetype=3) +
  ylab('Cramer statistics') +
  xlab('chain length')
ggsave(paste("plots/bivariate/Cramer",
             Sys.Date(),
             n_max + interval + n_chains,
             ".png", sep = "_"))


# Smirnov test
ggplot(NULL, aes(x=x_axis, y=smirnov)) +
  geom_line() +
  ylab('Smirnov statistics') +
  xlab('chain length')
ggsave(paste("plots/bivariate/Smirnov",
             Sys.Date(),
             n_max + interval + n_chains,
             ".png", sep = "_"))


# WWolfewitz test
ggplot(NULL, aes(x=x_axis, y=wwolf)) +
  geom_line() +
  ylab('Wild-Wolfewitz statistics') +
  xlab('chain length')
ggsave(paste("plots/bivariate/WWolf",
             Sys.Date(),
             n_max + interval + n_chains,
             ".png", sep = "_"))


# GelmanRubin diagnostic
ggplot(NULL, aes(x=x_axis, y=gelman)) +
  geom_line() +
  ylab('Gelman-Rubin statistics') +
  xlab('chain length')
ggsave(paste("plots/bivariate/Gelman",
             Sys.Date(),
             n_max + interval + n_chains,
             ".png", sep = "_"))


# Naive error
ggplot(NULL, aes(x=x_axis, y=naive_dim1), linetype=1) +
  geom_line() +
  geom_line(aes(y=naive_dim2), linetype=2) +
  ylab('Naive error') +
  xlab('chain length')
ggsave(paste("plots/bivariate/Naive",
             Sys.Date(),
             n_max + interval + n_chains,
             ".png", sep = "_"))
