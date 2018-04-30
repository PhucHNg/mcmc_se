#Author: Phuc Nguyen
#Date: 03/04/2018
#Description: Probability Density Functions to be used in simulation analysis


#function declaration
#=============================================================================

#' f is the target distribution: N(0, 1^2)
#' @param x any real number
#' returns P(X=x)
df <- function(x) {
  dnorm(x, mean = 0, sd = 1)
}

pf <- function(x, lower.tail=TRUE, log.p=FALSE) {
  pnorm(x, mean = 0, sd = 1, lower.tail = lower.tail, log.p = log.p)
}

#' q is the proposal distribution: N(x, 2^2)
#' @param x any real number used as the mean of q
#' returns a sample from q(y|x) 
rq <- function(x) {
  rnorm(1, mean = x, sd = 2)
}

#' e is an exponential target distribution: Exp(1)
#' @param x any real number
#' returns P(X=x)
de <- function(x) {
  dexp(x, rate = 1)
}

pe <- function(x, lower.tail=TRUE, log.p=FALSE) {
  pexp(x, rate = 1, lower.tail = lower.tail, log.p = log.p)
}

#' bimod is an bimodal target distribution: w(N(m1, v1)) + (1-w)(N(m2, v2))
#' where w = 0.3
#'       m1 = 1
#'       m2 = 5
#'       v1 = v2 = 1^2
#' @param x any real number
#' returns P(X=x)
dbimod <- function(x) {
  w = 0.3
  dn1 = dnorm(x, mean = 1, sd = 1)
  dn2 = dnorm(x, mean = 5, sd = 1)
  return( w*dn1 + (1-w)*dn2)
}

pbimod <- function(x, lower.tail=TRUE, log.p=FALSE) {
  w = 0.3
  pn1 = pnorm(x, mean = 1, sd = 1, lower.tail = lower.tail, log.p = log.p)
  pn2 = pnorm(x, mean = 5, sd = 1, lower.tail = lower.tail, log.p = log.p)
  return( w*pn1 + (1-w)*pn2)
}




