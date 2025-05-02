STRAND
========
 Social network analysis and simulation in R using Stan 
 ------

 <img align="left" src="logo.gif" alt="logo" width="200"> 
 <img align="right" src="strand.gif" alt="logo" width="250">

**STRAND** is an R package designed for simulation and analysis of network data.  The package can be used to simulate network data under stochastic block models, social relations models, and latent network models. These tools allow for simulation not only of realistic human and animal social networks, but also allow researchers to simulate the effects of potential biases—like respondents falsely reporting ties or failing to recall real ties—on network level properties. 

**STRAND** also allows users to conduct data analysis. Users can specify complex Bayesian social network models using a simple, lm-style, syntax. Single-sampled self report data can be modeled using stochastic block models or the Social Relations Model, with or with-out covariates. Double-sampled network data can be modeled using a latent network approach that accounts for inter-respondent disagreement. STRAND also provides methods for more rigorously dealing with missing data and measurement error. We have recently added support for multiplex network models and longitudinal network models. Here we provide a brief overview of various STRAND workflows. For further details, see our full publications at [**Psychological Methods**](https://psycnet.apa.org/record/2023-51200-001) (for latent network models), at [**Journal of Animal Ecology**](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.14021) (for simple network models), and at [**Methods in Ecology and Evolution**](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.70017) (for censoring-robust methods in ecology). 


  
[**STRAND**](https://github.com/ctross/STRAND) is part of an ecosystem of tools for modern social network analysis. [**DieTryin**](https://github.com/ctross/DieTryin) is a companion package designed to facilitate the collection of [**RICH economic games**](https://journals.sagepub.com/doi/abs/10.1177/1525822X16643709), dyadic peer ratings, and roster-based network data in human communities. [**XLSFormulatoR**](https://github.com/ADR1993/XLSFormulatoR) is a package for automatically building name-generator network surveys for KoboToolbox.



Install:
--------------
Install by running on R:
```{r}
################################### Install the latest release
 library(devtools)
 install_github('ctross/STRAND@phosphorescent_desert_buttons')
 library(STRAND)
```

You will need to have [**cmdstanr**](https://mc-stan.org/cmdstanr/) installed.

Quickstart guides can be found [**here.**](https://mc-stan.org/cmdstanr/articles/cmdstanr.html). 

**STRAND** calls **Stan** models in the background, so you will need a C++ compiler in your toolchain. Users frequently rely on [**RTools**](https://cran.r-project.org/bin/windows/Rtools/).

Uses:
========
Example models, with test data, can be found here. Note that analysis of large network data sets can be quite slow using MCMC.  An old adage, however, is that we often spend months or years collecting our data, and so we should be happy waiting a few hours for our models to fit. Our tutorials below are designed to accomplish a few things: (1) we teach users how to fit different kinds of network models using example data sets, (2) we compare STRAND to other software packages for network analysis, and comment on similarities and differences, and (3) we try to clarify any misconceptions potential users might have about Bayesian network analysis models---e.g, we show that Bayesian models with proper priors are always identifiable, that probit versus logit models are essentially equivalent, that STRAND's implementation of the SRM (Social Relations Model) is a measurement-error robust generalization of the original SRM model.


Basic single-layer network models
--------------
The most basic network analysis models are generalizations of the SRM that account for dyadic reciprocity (i.e., correlation in dyadic connections) and generalized reciprocity (i.e., correlations in nodal degree) through correlated dyadic and node-level random effects. Users can further account for focal, target, dyadic, and block-level covariates using simple lm-style syntax. The tutorial code below can be adapted as needed for users' own projects:

Binary/Bernoulli outcomes (using human friendship data from rural Colombia): [**Bernoulli models**](https://github.com/ctross/STRAND/blob/main/tutorials/Bernoulli_Example.R)

Binomial outcomes (using data from baboons): [**Binomial models**](https://github.com/ctross/STRAND/blob/main/tutorials/Binomial_Example.R)

Binomial outcomes with censoring bias (tested using simulated data): [**Binomial + ME models**](https://github.com/ctross/STRAND/blob/main/tutorials/Binomial_Measurement_Error_Example.R)

Poisson outcomes (using data from vampire bats): [**Poisson models**](https://github.com/ctross/STRAND/blob/main/tutorials/Poisson_Example.R)

Gaussian outcomes (using simulated data): [**Gaussian models**](https://github.com/ctross/STRAND/blob/main/tutorials/Gaussian_Example.R)

Undirected networks (tested using simulated data): [**Undirected networks**](https://github.com/ctross/STRAND/blob/main/tutorials/Undirected_Example.R)



Double-sampled single-layer network models
--------------
Respondents often disagree on ties. Alice might report giving to Bob, but Bob might not report receiving from Alice. To make sense of network data with such potential discordance, we integrate a reporting model into the SRM. Users can explicitly model the predictors of false positive reports, the recall rate of true ties, and inter-question duplication rates,  as well as all of the standard predictors using lm-style syntax.

Double-sampled binary outcomes (self-reported food/money sharing from rural Colombia): [**Latent network models**](https://github.com/ctross/STRAND/blob/main/tutorials/LatentNetwork_Example.R)



Multiplex network models
--------------
Networks typically influence each other. Users may wish to study if outgoing ties in one network layer are predictive of incoming ties in a different layer. To model such structure, we provide multiplex generalizations of the SRM for various types of outcomes. Users can still estimate the effects of focal, target, dyadic, and block predictors within each layer, while also estimating residual correlations in random effects within and across network layers at both a generalized and dyadic level (see our application paper at [**Communications Psychology**](https://www.nature.com/articles/s44271-024-00098-1) for a primer on these models). 

Multiplex Binary/Bernoulli outcomes (using RICH economic games): [**Multiplex Bernoulli models**](https://github.com/ctross/STRAND/blob/main/tutorials/Multiplex_Bernoulli_Example.R)

Multiplex Binomial outcomes (using baboon data): [**Multiplex Binomial models**](https://github.com/ctross/STRAND/blob/main/tutorials/Multiplex_Binomial_Example.R)

Multiplex Poisson outcomes (using simulated data): [**Multiplex Poisson models**](https://github.com/ctross/STRAND/blob/main/tutorials/Multiplex_Poisson_Simulation_Example.R)

Multiplex Gaussian outcomes (using simulated data): [**Multiplex Gaussian models**](https://github.com/ctross/STRAND/blob/main/tutorials/Multiplex_Gaussian_Simulation_Example.R)

Multiplex network modeling requires estimation of a highly-structured dyadic correlation matrix with special symmetries. We have two approaches for constructing such block-structured matrices, one based on an $\ell^{2}$
 norm penalty, and one based on constructing a Cholesky factor with a special set of constraints. Here we show that both methods are equally effective: [**Different methods for estimating dyadic reciprocity**](https://github.com/ctross/STRAND/blob/main/tutorials/Multiplex_Bernoulli_Example_Dyadic_Modes.R)


 Longitudinal network dynamics (experimental)
--------------
STRAND now supports longitudinal network analysis. These models draw on a multiplex network analysis framework to study network evolution over time. Longitudinal network models can be thought of as multiplex models with additional symmetries that arise from assuming transportability of effects across time-steps. See details in our preprint. These methods are currently passing the basic test fits, but are considered experimental. We appreciate any feedback or test reports. 

Longitudinal Binary/Bernoulli outcomes (using Colombian friendship data): [**Longitudinal Bernoulli models**](https://github.com/ctross/STRAND/blob/main/tutorials/Longitudinal_Bernoulli_Example.R)

Longitudinal Binomial outcomes (using baboon data): [**Longitudinal Binomial models**](https://github.com/ctross/STRAND/blob/main/tutorials/Longitudinal_Binomial_Example.R)

Longitudinal simulation analysis (using simulated data): [**Longitudinal Generative Simulations**](https://github.com/ctross/STRAND/blob/main/tutorials/Longitudinal_Simulation_Example.R)


Automatic Bayesian imputation of missing data (experimental)
--------------
STRAND now support automatic Bayesian imputation for continuous predictor variables, and automatically slices missing outcomes out of the likelihood. Block predictors are currently imputed deterministically prior to model fitting via columnwise resampling (a more rigorous method may be integrated in the future). These methods are currently passing the basic unit tests, but are considered experimental. We appreciate any feedback or test reports. Once these models go through a decent burn-in period, we will push this functionality to the base functions.

Binary/Bernoulli outcomes with missings: [**Single-layer models with imputation**](https://github.com/ctross/STRAND/blob/main/tutorials/Bernoulli_Example_With_Imputation.R)

Binomial outcomes with missings: [**Binomial + censoring models with imputation**](https://github.com/ctross/STRAND/blob/main/tutorials/Binomial_Measurement_Error_Example_With_Imputation.R)

Double-sampled latent network models with missings: [**Latent network models with imputation**](https://github.com/ctross/STRAND/blob/main/tutorials/LatentNetwork_Example_With_Imputation.R)

Multiplex Poisson outcomes with missings: [**Multiplex models with imputation**](https://github.com/ctross/STRAND/blob/main/tutorials/Multiplex_Poisson_Example_With_Imputation.R)

Longitudinal Bernoulli outcomes with missings: [**Longitudinal models with imputation**](https://github.com/ctross/STRAND/blob/main/tutorials/Longitudinal_Bernoulli_Example_With_Imputation.R)







Miscellaneous features and toy examples
--------------
Below we provide a few more pointed tutorials showing how to deal with specific issues, like structural zeros, prior specification, data simulation, and comparison of STRAND to other tools. We also show a few other things about STRAND models, e.g., probit and logit links yield equivalent inference, binary SRM models are well-specified in a Bayesian framework, etc. We will also include some minimum working examples to address many common questions we get via email.

An example on both simulating and fitting networks (includes interactions): [**Simulating data and fitting interaction models**](https://github.com/ctross/STRAND/blob/main/tutorials/Interaction_Example.R)

An example on both simulating and fitting multiplex networks: [**Simulating data and fitting multiplex models**](https://github.com/ctross/STRAND/blob/main/tutorials/Multiplex_Binomial_Simulation_Example.R)

An example on accounting for structural-zeros/censored-data. For example, there may be N groups of individuals in a dataset, where each group is in a separate enclosure, and thus only with-group ties can be modeled: [**Structural Zeros**](https://github.com/ctross/STRAND/blob/main/tutorials/Bernoulli_Callithrix_Example.R).

An example on interactions between focal (sender), target (receiver), and dyadic effect: [**Between-level interaction**](https://github.com/ctross/STRAND/blob/main/tutorials/Bernoulli_Between_Level_Interaction_Example.R)

An example on changing default priors: [**Changing priors**](https://github.com/ctross/STRAND/blob/main/tutorials/ChangingPriors_Example.R)

Probit versus logit links for binary outcomes: [**Does anything depend on the choice of link function?**](https://github.com/ctross/STRAND/blob/main/tutorials/Bernoulli_Probit_v_Logit_Example.R)

Probit versus logit links for multiplex binary outcomes: [**Does anything depend on the choice of link function?**](https://github.com/ctross/STRAND/blob/main/tutorials/Multiplex_Bernoulli_Probit_v_Logit_Example.R)

STRAND's Gaussian SRM has both dyad-level random effects and dyad-level error. Can we still recover parameters?  [**Yes.**](https://github.com/ctross/STRAND/blob/main/tutorials/Gaussian_Measurement_Example.R)

STRAND's probit/logit SRM has both dyad-level random effects and dyad-level error. Can we still recover parameters?  [**Again, yes.**](https://github.com/ctross/STRAND/blob/main/tutorials/Probit_Measurement_Example.R)



Citations:
--------------
If you use STRAND, please cite us using the most relevant paper:
```{bibtex}
@article{ross2024modelling,
  title={Modelling animal network data in R using STRAND},
  author={Ross, Cody T and McElreath, Richard and Redhead, Daniel},
  journal={Journal of Animal Ecology},
  volume={93},
  number={3},
  pages={254--266},
  year={2024},
  publisher={Wiley Online Library}
}
```

```{bibtex}
@article{redhead2023reliable,
  title={Reliable network inference from unreliable data: A tutorial on latent network modeling using STRAND},
  author={Redhead, Daniel and McElreath, Richard and Ross, Cody T},
  journal={Psychological methods},
  year={2023},
  volume={29},
  number={6},
  pages={1100--1122},
  publisher={American Psychological Association}
}
```

```{bibtex}
@article{sosa2024robust,
 author = {Sosa, Sebastian and McElreath, Mary B. and Redhead, Daniel and Ross, Cody T.},
 title = {Robust Bayesian analysis of animal networks subject to biases in sampling intensity and censoring},
 journal = {Methods in Ecology and Evolution},
 volume = {},
 number = {},
 pages = {1--22},
 doi = {https://doi.org/10.1111/2041-210X.70017}
}
```

Note:
========
Each of the models included in this package have been fit to real empirical datasets, and tested across a wide-range of simulated data to ensure their quality. However, this package is still rather new. If you come across any weird behavior, or notice any bugs, please open an issue, and we will work to address it! 

Additionally, anyone can issue pull requests. If you notice any typos in the documentation, or feel like you can add something useful, please submit a pull request with your proposed changes. We will inspect the changes closely and integrate them if they are helpful.




