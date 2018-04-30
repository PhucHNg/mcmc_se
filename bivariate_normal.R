#' #Author: Phuc Nguyen
#' Date: 04/20/2018
#' Description: Multivariate normal density functions


#' Returns a target density function
#' @param x: equivalent to f(x) where f is pdf of a bivariate normal
#' @param mean: 2-d vector, mean of dist
#' @param covariance: 2x2 matrix
target_wrapper <- function(mean, covariance) {
  inverse_cov <- solve(covariance)
  target_dnormm <- function(x) {
    exp(-0.5*(t(x-mean)%*%inverse_cov%*%(x-mean)))
  } 
  return(target_dnormm)
}


#' Returns a random sampler from a bivariate normal
#' with given covariance and mean
proposal_wrapper <- function(covariance) {
  chol_covariance <- t(chol(covariance))
  proposal_rnormm <- function(mean) {
    sample <- mean + chol_covariance%*%rnorm(2) # 2 because this is bivariate
    return(sample)
  }
  return(proposal_rnormm)
}



