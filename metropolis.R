#Author: Phuc Nguyen
#Date: 02/09/2018
#Description: Univariate Metropolis algorithm

#'Function simulates a Metropolis MCMC
#'
#' @param f The pdf of target "posterior" distribution  
#' @param q The random sampler of the proposal distribution
#' @param chain_length The length of MCMC chain
#' @param y_0 The initial value
#' @return a 1-d array of MCMC samples
Metropolis <- function(f , q, chain_length=100, y_0=1) {
  
  mcmc <- rep(0, chain_length)
  current <- y_0
  for (i in 1:chain_length) {
    #draw a proposal from the q
    proposal <- q(current)
    
    #calculate the proposal acceptance probability -- alpha 
    alpha <- min(1, f(proposal)/f(current))
    
    #with probability alpha, move to the proposal
    next_val <- sample(x=c(current, proposal), size=1, prob = c(1-alpha, alpha))
    
    #store proposals and stops
    current <- next_val
    mcmc[i] <- current
  }
  
  return(mcmc)
}