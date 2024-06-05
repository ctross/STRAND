STRAND
========
 Social network analysis and simulation in R using Stan 
 ------

<img align="right" src="https://github.com/ctross/STRAND/blob/main/logo.png" alt="logo" width="200"> 


**STRAND** is an R package designed for simulation and analysis of network data.  The package can be used to simulate network data under stochastic block models, social relations models, and latent network models. These tools allow for simulation not only of realistic human and animal social networks, but also allow researchers to simulate the effects of potential biases—like respondents falsely reporting ties or failing to recall real ties—on network level properties. 

**STRAND** also allows users to conduct data analysis. Users can specify complex Bayesian social network models using a simple, lm() style, syntax. Single-sampled self report data can be modeled using stochastic block models or the social relations model, with or with-out covariates. Double-sampled network data can be modeled using a latent network approach that accounts for inter-respondent disagreement. Here we provide a brief overview of a workflow. For further details on any step, see our full publication in X at <>.
  
[**STRAND**](https://github.com/ctross/STRAND) is part of an ecosystem of tools for modern social network analysis. [**DieTryin**](https://github.com/ctross/DieTryin) is a companion package designed to facilitate the collection of roster-based network data, and to run network-structured economic games. **ResolveR** (currently under development) is a package for semi-supervised data cleaning, de-duplication, and record linkage.

Install:
--------------
Install by running on R:
```{r}
################################### Install and/or load
 library(devtools)
 install_github('ctross/STRAND@measurement_error')
 library(STRAND)
```

You will need to have [**cmdstanr**](https://mc-stan.org/cmdstanr/) and [**rstan**](https://mc-stan.org/users/interfaces/rstan) installed.

Quickstart guides for these packages can be found [**here, for cmdstanr**](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) and [**here, for rstan**](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). 

**STRAND** calls **Stan** models in the background, so you will need a C++ compiler in your toolchain. Users frequently rely on [**RTools**](https://cran.r-project.org/bin/windows/Rtools/).

Finally, we note that some users will see a warning message: "Declaration of arrays by placing brackets after a variable name is deprecated and will be removed in Stan 2.32.0." This message is ignorable, we will update the Stan code inside STRAND once deprication occurs at Stan 2.32.0.

Use:
--------------

Example models, with test data, can be found here. Note that analysis of large network data sets can be quite slow using MCMC. 

Binary outcomes (runs in about 2 hours): [**Bernoulli models**](https://github.com/ctross/STRAND/blob/main/tutorials/Bernoulli_Example.R)

Binomial outcomes (runs in about 2 minutes): [**Binomial models**](https://github.com/ctross/STRAND/blob/main/tutorials/Binomial_Example.R)

Binomial outcomes with measurement error (runs in about 20 minutes): [**Binomial + ME models**](https://github.com/ctross/STRAND/blob/main/tutorials/Binomial_Measurement_Error_Example.R)

Count outcomes (runs in about 20 seconds): [**Poisson models**](https://github.com/ctross/STRAND/blob/main/tutorials/Poisson_Example.R)

Double-sampled binary outcomes (find out for yourself): [**Latent network models**](https://github.com/ctross/STRAND/blob/main/tutorials/LatentNetwork_Example.R)

An example on both simulating and fitting networks (includes interactions): [**Interaction models**](https://github.com/ctross/STRAND/blob/main/tutorials/Interaction_Example.R)

An example on changing default priors: [**Changing priors**](https://github.com/ctross/STRAND/blob/main/tutorials/ChangingPriors_Example.R)

Note: each of the models included in this package have been fit to real emprical datsets, and tested across a wide-range of simulated data to ensure their quality. However, this package is still rather new. If you come across any weird behavior, or notice any bugs, please open an issue, and we will work to address it!
