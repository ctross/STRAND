STRAND
========
 Social network analysis and simulation in R using Stan 
 ------

<img align="right" src="https://github.com/ctross/STRAND/blob/main/logo3.png" alt="logo" width="140"> 
<img align="right" src="https://github.com/ctross/STRAND/blob/main/logo2.png" alt="logo" width="140">
<img align="right" src="https://github.com/ctross/STRAND/blob/main/logo.png" alt="logo" width="140">

**STRAND** is an R package designed for simulation and analysis of network data.  The package can be used to simulate network data under stochastic block models, social relations models, and latent network models. These tools allow for simulation not only of realistic human and animal social networks, but also allow researchers to simulate the effects of potential biases—like respondents falsely reporting ties or failing to recall real ties—on network level properties. 

**STRAND** also allows users to conduct data analysis. Users can specify complex Bayesian social network models using a simple, lm() style, syntax. Single-sampled self report data can be modeled using stochastic block models or the social relations model, with or with-out covariates. Double-sampled network data can be modeled using a latent network approach that accounts for inter-respondent disagreement. Here we provide a brief overview of a workflow. For further details on any step, see our full publication in X at <>.
  
[**STRAND**](https://github.com/ctross/STRAND) is part of an ecosystem of tools for modern social network analysis. [**DieTryin**](https://github.com/ctross/DieTryin) is a companion package designed to facilitate the collection of roster-based network data, and to run network-structured economic games. **ResolveR** (currently under development) is a package for semi-supervised data cleaning, de-duplication, and record linkage.

Use:
--------------
Install by running on R:
```{r}
################################### Install and/or load
 library(devtools)
 install_github('ctross/STRAND', ref = "major-update-ab1")
 library(STRAND)
```

Example models, with test data, can be found here. 

Binary outcomes: [**Bernoulli models**](https://github.com/ctross/STRAND/blob/ctross-patch-1/Bernoulli_Example.R)

Binomial outcomes: [**Binomial models**](https://github.com/ctross/STRAND/blob/ctross-patch-1/Binomial_Example.R)

Count outcomes: [**Poisson models**](https://github.com/ctross/STRAND/blob/ctross-patch-1/Poisson_Example.R)

Double-sampled binary outcomes: [**Latent network models**](https://github.com/ctross/STRAND/blob/ctross-patch-1/LatentNetwork_Example.R)

Note: each of the models included in this package have been fit to real emprical datsets, and tested across a wide-range of simulated data to ensure their quality. However, this package is still rather new. If you come across any weird behavior, or notice any bugs, please open an issue, and we will work to address it!
