#Author: Phuc Nguyen
#Date: 02/09/2018
#Description: Bootstrap methods

#'Function takes N bootstrap samples of x
#'
#' @param x A vector of values to be resampled
#' @param N number of bootstrap samples to return
Bootstrap <- function(x, N=1000) {
  
  size <- length(x)
  boots <- matrix(0, ncol = N, nrow = size)
  for (i in 1:N) {
    boots[,i] <- sample(x, size = size, replace = TRUE)
  }
  
  return(data.frame(boots))
}


#'Function plots density of MCMC estimated posterior, true target distribution, and bootstrap density
#' @param org Vector of the original MCMC chain 
#' @param boots Data frame containing bootstrap resamples
#' @param fun Function of the true target distribution
#' @param args List of additional arguements for the true target distribution
plot_resamples <- function(org, boots, fun, args=list(), xlim=NULL) {
  
  N <- dim(boots)[2]
  #density of original MCMC estimate
  g <- ggplot(data.frame(org=org), aes(x=org))
  #density of resamples
  for (i in seq(1,N,5)) {
    g <- g + geom_density(data=data.frame(X=boots[,i]), aes(x=X), colour='lightskyblue', linetype=2, size=0.6)
  }
  
  g <- g + geom_density(fill='dodgerblue', colour=NA, alpha=0.4)
  #density of true target distribution
  g <- g +
    stat_function(fun = fun, 
                  args = args, 
                  colour = 'red', size = 1)
  
  if (!is.null(xlim)) {
    g <- g + xlim(xlim)
  } else {
    g <- g + xlim(c(-3,3))
  } 
  
  return(g)
}


#' Helper function to calculate ks statistics between resamples
#' @param i index of reference resample
#' @param boots data frame of bootstrap resamples
#' return a list of ks.test results
ks.resamp <- function(i, boots) {
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


#' Helper function to extract ks-statistics and p-values
#' @param ks_list list of result from ks.test
#' returns list of statistics and p-values
get.ks.stats <- function(ks_list) {
  n <- length(ks_list)
  ks_dist <- c()
  ks_pval <- c()
  for (i in 1:n) {
    ks_dist[i] <- ks_list[[i]]$statistic
    ks_pval[i] <- ks_list[[i]]$p.value
  }
  return(list(dist=ks_dist, pval=ks_pval))
}


#'Function calculates Kolmogorov-Smirnov test statistics between resamples
#' @param boots Dataframe containing bootstrap resamples
#' @param chain Original MCMC chain, default is null. If not null, also 
#' calculate KS statistics between resamples and the original chain
my.ks.test <- function(boots, chain=NULL) {
  result <- list()
  
  if (! is.null(chain)) {
    KS_org <- apply(boots, 2, ks.test, chain)
    list[KS_org_dist, KS_org_pval] <- get.ks.stats(KS_org)
    result[['ks.original.stats']] <- KS_org_dist
    result[['ks.original.pval']] <- KS_org_pval
  }
  
  N <- dim(boots)[2]
  i <- 1:(N-1)
  KS_resample <- unlist(lapply(i, ks.resamp, boots), recursive = FALSE)
  list[KS_resample_dist, KS_resample_pval] <- get.ks.stats(KS_resample)
  result[['ks.resample.stats']]<- KS_resample_dist
  result[['ks.resample.pval']] <- KS_resample_pval
  
  return(result)
}
