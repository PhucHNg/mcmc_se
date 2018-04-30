#' Generate an MCMC sample of a bivariate normal with given mean and covariance.
#' Implementation taken from http://www.math.unm.edu/~ghuerta/sta590/handout9.pdf
#' 

#' Metropolis Hasting for bivariate normal
#' @param chain_length: number of iterations of MCMC
#' @param target: probability of x in target distribution
#' @param proposal: draws a random sample from proposal distribution, takes mean
#' @param y_0: initial condition
#' @return matrix of samples, chain_length x 2 in dimension
Metropolis_m <- function(target, proposal, chain_length, y_0) {
  y <- y_0
  chains <- y
  for (i in 1:chain_length) {
    # draw a proposal from proposed distribution centered at current y
    y_new <- proposal(y)
    acceptance <- min(c(1, target(y_new)/target(y)))
    u <- runif(1)
    if (u < acceptance) {
      y <- y_new
    }
    chains <- c(chains, y)
  }
  return(matrix(chains, (chain_length+1), 2, byrow=TRUE))
}