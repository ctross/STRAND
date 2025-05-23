% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_NBDA_model.R
\name{fit_NBDA_model}
\alias{fit_NBDA_model}
\title{A function to run a network-based diffusion analysis using the STRAND framework}
\usage{
fit_NBDA_model(
  long_data,
  individual_focal_regression,
  social_block_regression,
  social_focal_regression,
  social_target_regression,
  social_dyad_regression,
  mode = "mcmc",
  stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1,
    iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL,
    init = NULL),
  priors = NULL
)
}
\arguments{
\item{long_data}{A list of data objects of class STRAND prepared using the make_strand_data() function. The data objects must include all covariates and trait diffusion data used in the formulas listed below.}

\item{individual_focal_regression}{A formula for the effects of focal predictors on individual learning rate. This should be specified as in lm(), e.g.: ~ Age * Education.}

\item{social_block_regression}{A formula for the block-level predictors of social attention weights. This should be specified as in lm(), e.g.: ~ Group + Sex. Dont use interactions, however.}

\item{social_focal_regression}{A formula for the effects of focal predictors on social attention weights. This should be specified as in lm(), e.g.: ~ Age * Education.}

\item{social_target_regression}{A formula for the effects of target predictors on social attention weights. This should be specified as in lm(), e.g.: ~ Age * Education.}

\item{social_dyad_regression}{A formula for the predictors of dyadic relationships on social attention weights. This should be specified as in lm(), e.g.: ~ Kinship + Friendship.}

\item{mode}{A string giving the mode that stan should use to fit the model. "mcmc" is default and recommended, and STRAND has functions to make processing the mcmc samples easier. Other options are "optim", to
use the optimizer provided by Stan, and "vb" to run the variational inference routine provided by Stan. "optim" and "vb" are fast and can be used for test runs. To process their output, however,
users must be familar with [cmdstanr](https://mc-stan.org/users/interfaces/cmdstan). We recommmend that users refer to the [Stan user manual](https://mc-stan.org/users/documentation/) for more information about the different modes that Stan can use.}

\item{stan_mcmc_parameters}{A list of Stan parameters that often need to be tuned. Defaults set to: list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
We recommend 1000 sampling and warmup iterations on a single chain for exploratory model fitting. For final runs, we recommend running 2 to 4 chains for twice as long. Be sure to check r_hat, effective sample size, and traceplots.}

\item{priors}{A labeled list of priors for the model. User are only permitted to edit the values. Distributions are fixed.}
}
\value{
A STRAND model object containing the data used, and the Stan results.
}
\description{
This function allows users to analyse empirical or simulated data using a NBDA model in Stan. The user must supply a list of STRAND data objects,
and a series of formulas following standard lm() style syntax.
}
\examples{
\dontrun{
fit = fit_NBDA_model(long_data=model_dat,
                     individual_focal_regression = ~ Age * NoFood,
                     social_block_regression = ~ Ethnicity,
                     social_focal_regression = ~ Age * NoFood,
                     social_target_regression = ~ Age * NoFood,
                     social_dyad_regression = ~ Relatedness + Friends * SameSex,
                     mode="mcmc",
                     stan_mcmc_parameters = list(seed = 1, chains = 1, 
                       parallel_chains = 1, refresh = 1, 
                       iter_warmup = 100, iter_sampling = 100,
                       max_treedepth = NULL, adapt_delta = NULL)
                      )
}

}
