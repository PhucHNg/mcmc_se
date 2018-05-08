# In Search of a Bayesian Stopping Rule for Markov Chain Monte Carlo Samplers

## Introduction
 
Markov Chain Monte Carlo (MCMC) is an essential technique in estimating parameters of Bayesian models. 
Given a large enough number of sample size, MCMC sampler has been showed to converge to the 
true posterior distribution (CITE STH). However, in practice, it's still difficult to determine
how large is a large sample size.
Many commonly used diagnostics and stopping rules for MCMC chains
tend to only focus on the stablization of the posterior mean estimates (CITE REVIEW PAPER BY COWELS). 
In this project, we test different approaches to estimate the stablization of the estimate of the whole posterior distribution with toy univariate and multivariate examples.
Specifically, we attempted to measure this "stablization" by:

+ Comparing the distributions of multiple MCMC samples using Kolmogorov-Smirnov test statistics in univariate case, and more general multivariate two sample tests in bivariate case.

+ Calculating batch mean standard error at multiple quantiles of a single MCMC sample, for both univariate and bivariate case.

+ Comparing the distributions of bootstrap resamples of a single MCMC sample using Kolmogorov-Smirnov test statistics for different univariate distributions.

## Simulation results

### Convergence of multiple samples

The idea behind this experiment is that two MCMC samples, starting at different (ideally overdispersed) intital conditions, will appear to be two samples drawn from the same distribution when they have both converged to the true posterior density. The Kolmogorov-Smirnov (KS) (LINK TO KS WIKI) test statistics is an efficient measure to test the null hypothesis that two given samples are from the same density function. This test is particularly useful in this context because it does not make any assumption about the underlying distribution. The p-value associated with this test statistics, however, is not useful because it's always significant when the sample size is large.

#### Univariate case

We run simulations on the following target distributions:

+ Normal: $N(0,1)$
+ Exponential: $Exp(1)$
+ Bimodal: of two overlapping normal distributions: $X = 0.3N(1,1) + 0.7N(5,1)$

and compare the stopping points given by the approach mentioned above, the Gelman-Rubin diagnostic (CITE CODA), and naive posterior mean standard error (CITE MCMCSE).

The approach mentioned above is implemented as follows: Determine a desirable value for the KS statistics (in this experiment, we chose 0.05) which equates to how similar MCMC samples' CDF's are. Since KS statistics is between 0 and 1, one might establish a standardized cut off value across parameters of different scales. Four MCMC chains are ran, starting at points spread out across the proposal distributions. A KS statistic is calculated every 100 time steps between every two chains. The maximum KS statistics is then recorded.

Target | $N(0,1)$                  | $Exp(1)$                  | $0.3N(1,1) + 0.7N(5,1)$
------:|:-------------------------:|:-------------------------:|:-----------------------------------:
       |![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_normal_1_.png)  |  ![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_exponential_1_.png) | ![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_bimodal_1_.png)
       |![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_normal_2_.png)  |  ![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_exponential_2_.png)  |  ![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_bimodal_2_.png)
       |  ![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_normal_3_.png)  |  ![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_exponential_3_.png)  | ![](https://github.com/PhucHNg/mcmc_se/blob/master/plots/univariate_ks/multichain_2018-04-16_20104_bimodal_3_.png)

[SHOWS THREE PLOTS SUPERIMPOSED STATISTICS FOR THREE TARGET DISTRIBUTION]

In general, KS statistics decreases as the sample size increases. Depending on the target distributions, sometimes the KS approach terminates the chains before or after the other standard approaches. [SOME MORE OBSERVATIONS HERE]

#### Bivariate case

We extend the above experiment to higher dimension using a bivariate normal whose dimensions are uncorrelated and one whose dimensions are correlated. We use three test statistics that are analogous to the KS test for multivariate samples:

+ Cramer test (CITE PACKAGE OR PAPER)
+ Multivariate Smirnov test using Minimum Spanning Tree (CITE PAPER AND PACKAGE)
+ Wald-Wolfewitz test using Minimum Spanning Tree (CITE PAPER AND PACKAGE)

These test statistics prove to be computationally expensive to calculate. We use only two chains to calculate the statistics. The Cramer test starts taking a noticable amount of time to calculate statistics for sample size of more than 600, and RStudio crashes when the sample size is above 2000. The other two tests based on MST take even a longer to compute than the Cramer test so are eliminated early on from further experimentation. To at least have a proof of concept for our KS approach in multidimension, thinning is used to reduce the MCMC sample size. Only one sample is kept for every 20 samples. Only Cramer test is performed on the two chains. (Thinning makes the smaller samples approximately independent, which fits the i.i.d. assumption of Cramer test a little better.)

:-------------------------:|:-------------------------:


[SHOWS THREE PLOTS SUPERIMPOSED STATISTICS FOR THREE TARGET DISTRIBUTION]

#### Discussion

Though the idea of using two sample test statistics to compare chains' distribution seem reasonable, it is currently not scalable to higher dimension examples. KS-like tests for multivariate examples. This direction can be further explored if there were other generalizations of KS test in higher dimension with efficient implementations. Another drawback of this approach is that running multiple chains is expensive if run sequentially since one can run one chain for twice as long and get a better estimate of the posterior density. However, (I think) parallel computing can be used to sample independent chains simulatneously at not much more computational cost [R-blogger](https://www.r-bloggers.com/post-10-multicore-parallelism-in-mcmc/).

### Reduction of standard error across quantiles

A common stopping rule is determining a level of accuracy and stop the chains whenever the batch mean standard error estimate (CITE PAPER AND PACKAGE) is below this level. This standard error is usually calculated for the posterior mean, though there also exists methods to estimate standard error of user-specified quantile. Since we are interested in the convergence of the whole posterior density, we propose extending the above stopping rule to include standard error estimates of all quantiles. We implemented our idea as follows: Determine a desired level of accuracy specific to the application domain. Determine the quantiles at which s.e. of posterior estimate is calculated. For every 100 time steps, we calculate the s.e. for all quantiles determined in previous step and record the maximum s.e. MCMC chain can in theory be stopped after the maximum s.e. across the whole posterior estimate is below the accuracy level.

[DISCUSSION OF THE RESULT HERE]

#### Discussion:

Calculating the standard error for every quantiles is computationally intensive. As seen above, for a simple standard normal toy example, we need to check ____ quantiles to get a good--stable across different runs of the experiments--estimate of the maximum standard error across quantiles. Depending on the problem, we might need to check an even bigger number of quantiles to ensure the resulting maximum s.e. curve is somewhat "smooth". We plan to experiment with an adaptive scheme--checking more quantiles where the s.e. estimates are large--to reduce the amount of computation.

In addition, one might standardize the level of accuracy across parameters of different scaltes by using a standardized standard error as proportion of the true standard deviation (CITE PAPER) instead of
a consistent batch mean s.e. Implementation of this mcmc s.e. estimate for the mean and quantiles doesn't exist yet to my knowledge. Similarly, an s.e. estimate that accounts for auto-correlation between dimensions can further improve quality of our stopping rule.

### Variation among bootstrap resamples

Instead of running multiple chains and comparing them with KS test, we experiment with taking bootstrap resamples of a single chain and calculating a KS statistics between every pair for every 100 time steps. The maximum KS statisitics is recorded. However, since all the resamples are drawn from the same distribution, the KS test ends up saying that all our bootstrap resamples have the same underlying distribution. 

[PICTURE OF FAILED BOOTSTRAP HERE]

## Reference


