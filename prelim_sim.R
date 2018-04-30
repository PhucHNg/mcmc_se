#Author: Phuc Nguyen
#Date: 02/09/2018
#Title: Preliminary Simulation Study

#import modules and packages
#========================================================================
source("utils.R")
library(gsubfn)
library(ggplot2)

#simulations
#========================================================================
#Let N(0, 1^2) be the target "posterior" distribution.
#    N(0, 2^2) be the first proposal distribution.
#    N(0, 0.1^2) be the second proposal distribution.
#Bootstrap MCMC chains (length of 100 and 1000) for each proposal.

set.seed(2018)

#f is the target distribution: N(0, 1^2)
#returns P(X=x)
df <- function(x) {
  dnorm(x, mean = 0, sd = 1)
}

#q1 is the proposal distribution: N(0, 2^2)
#returns a sample from q(x|y) 
rq1 <- function(x) {
  rnorm(1, mean = x, sd = 2)
}

#q2 is the proposal distribution: N(0, 0.1^2)
#returns a sample from q2(x|y) 
rq2 <- function(x) {
  rnorm(1, mean = x, sd = 0.1)
}

#MCMC chains of length 100 and 1000 for q1 and q2 proposal distributions
qs <- list(rq1, rq2)
qs_description <- c('norm0-2', 'norm0-0.1')
lens <- c(100, 1000)
N <- 1000

for (i in 1:length(qs)) {
  q <- qs[[i]]
  des <- qs_description[i]
  for (l in lens) {
    #Metropolis MCMC
    chain <- Metropolis(df, q, chain_length = l)
    
    #Take 1000 bootstrap resamples
    boots <- Bootstrap(chain, N = N)
    
    #Visualize the posterior approximation on top of the true distribution
    filename <- paste0('plots/prelim_', des, '_', l,'.pdf')
    g <- plot_resamples(chain, boots, fun=dnorm, args=list(mean=0, sd=1)) #N(0,1^2) is the target distribution
    ggsave(filename, plot=g)
  }
}

