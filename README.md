# In Search of a Bayesian Stopping Rule for Markov Chain Monte Carlo Samplers

## Introduction
 
Markov Chain Monte Carlo (MCMC) is an essential technique in estimating posteriors of Bayesian models. 
Given a large enough number of sample size, MCMC sampler has been showed to converge to the 
true posterior distribution. However, in practice, it's still difficult to determine
how large is a large sample size.
Existing diagnostics and stopping rules for posterior estimates using MCMC
tend to only focus on the stablization of the posterior mean estimates. In this project, we test 
different approaches to estimate the stablization of the estimate of the whole posterior distribution with toy univariate and multivariate examples.
Specifically, we attempted to measure this "stablization" by:

+ Comparing the distributions of multiple MCMC samples using Kolmogorov-Smirnov test statistics in univariate case, and its variations in bivariate case.

+ Calculating batch mean standard error at multiple quantiles of a single MCMC sample.

+ Comparing the distributions of bootstrap resamples of a single MCMC sample using Kolmogorov-Smirnov test statistics.

## Simulation results

### Convergence of multiple samples

The idea behind this test is that, given two MCMC samples starting at different intital conditions, they will appear to be two samples drawn from the same distribution when they have both converged to the true posterior density. The Kolmogorov-Smirnov (KS) test statistics is an efficient measure to test the null hypothesis that the two given samples are from the same distribution. This test is particularly useful in this context because it does not make any assumption about the underlying distribution.

#### Univariate case

We ran simulations on the following target distributions:

+ Normal: mean 0, standard deviation 1
+ Exponential: rate 1
+ Bimodal: of two overlapping normal distributions

and compared the stopping points given by the approach mentioned above, the Gelman-Rubin strinking factor diagnostic, and naive posterior mean standard error.

In general, KS statistics decreases as the sample size increases. Depending on the target distributions, sometimes our approach terminate the chain before or after the other standard approaches. 

#### Bivariate case

We ran another simulation using a bivariate normal whose dimensions are uncorrelated. We used three test statistics that are analogous to the KS test for multivariate samples:

+ Cramer test
+ Multivariate Smirnov test using Minimum Spanning Tree
+ Wald-Wolfewitz test using Minimum Spanning Tree

These test statistics proved to be computationally expensive to calculate.

### Reduction of standard error across quantiles

### Variation among bootstrap resamples

