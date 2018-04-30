# In Search of a Bayesian Stopping Rule for Markov Chain Monte Carlo Samplers

## Introduction
 
Markov Chain Monte Carlo (MCMC) is an essential technique in estimating posteriors of Bayesian models. 
Given a large enough number of sample size, MCMC sampler has been showed to converge to the 
true posterior distribution. However, in practice, it's still difficult to determine
how large is a large sample size.
Existing diagnostics and stopping rules for posterior estimates using MCMC
tend to only focus on the stablization of the posterior mean estimates. In this project, we test 
different approaches to 
