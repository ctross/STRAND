#' A function to run a network-based diffusion analysis using the STRAND framework
#' 
#' This function allows users to analyse empirical or simulated data using a NBDA model in Stan. The user must supply a list of STRAND data objects,
#' and a series of formulas following standard lm() style syntax. 
#'
#' @param long_data A list of data objects of class STRAND prepared using the make_strand_data() function. The data objects must include all covariates and trait diffusion data used in the formulas listed below.
#' @param individual_focal_regression A formula for the effects of focal predictors on individual learning rate. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_block_regression A formula for the block-level predictors of social attention weights. This should be specified as in lm(), e.g.: ~ Group + Sex. Dont use interactions, however.
#' @param social_focal_regression A formula for the effects of focal predictors on social attention weights. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_target_regression A formula for the effects of target predictors on social attention weights. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_dyad_regression A formula for the predictors of dyadic relationships on social attention weights. This should be specified as in lm(), e.g.: ~ Kinship + Friendship.
#' @param mode A string giving the mode that stan should use to fit the model. "mcmc" is default and recommended, and STRAND has functions to make processing the mcmc samples easier. Other options are "optim", to
#' use the optimizer provided by Stan, and "vb" to run the variational inference routine provided by Stan. "optim" and "vb" are fast and can be used for test runs. To process their output, however,
#' users must be familar with [cmdstanr](https://mc-stan.org/users/interfaces/cmdstan). We recommmend that users refer to the [Stan user manual](https://mc-stan.org/users/documentation/) for more information about the different modes that Stan can use. 
#' @param stan_mcmc_parameters A list of Stan parameters that often need to be tuned. Defaults set to: list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
#' We recommend 1000 sampling and warmup iterations on a single chain for exploratory model fitting. For final runs, we recommend running 2 to 4 chains for twice as long. Be sure to check r_hat, effective sample size, and traceplots.
#' @param priors A labeled list of priors for the model. User are only permitted to edit the values. Distributions are fixed. 
#' @return A STRAND model object containing the data used, and the Stan results.
#' @export
#' @examples
#' \dontrun{
#' fit = fit_NBDA_model(long_data=model_dat,
#'                      individual_focal_regression = ~ Age * NoFood,
#'                      social_block_regression = ~ Ethnicity,
#'                      social_focal_regression = ~ Age * NoFood,
#'                      social_target_regression = ~ Age * NoFood,
#'                      social_dyad_regression = ~ Relatedness + Friends * SameSex,
#'                      mode="mcmc",
#'                      stan_mcmc_parameters = list(seed = 1, chains = 1, 
#'                        parallel_chains = 1, refresh = 1, 
#'                        iter_warmup = 100, iter_sampling = 100,
#'                        max_treedepth = NULL, adapt_delta = NULL)
#'                       )
#' }
#' 

fit_NBDA_model = function(long_data,
                          individual_focal_regression,
                          social_block_regression,
                          social_focal_regression,
                          social_target_regression,
                          social_dyad_regression,
                          mode="mcmc",
                          stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL,
                                                      iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL, init = NULL),
                          priors=NULL
                          ){
    ############################################################################ Build data
    if(is.null(names(long_data))) stop("long_data must be a named list. Please add names for each time-step. e.g., names(long_data)=paste('Time', 1:T)")
    data = make_longitudinal_data(long_data = long_data,
                                  block_regression = social_block_regression,
                                  focal_regression = social_focal_regression,
                                  target_regression = social_target_regression,
                                  dyad_regression = social_dyad_regression
                                  )

    data2 = make_longitudinal_data(long_data = long_data,
                                   block_regression = ~ 1,
                                   focal_regression = individual_focal_regression,
                                   target_regression = ~ 1,
                                   dyad_regression = ~ 1
                                   )

    data$ind_focal_set = data2$focal_set
    data$N_params =  c(data$N_params, data2$N_params[1])

    ############################################################################# Check inputs
    if(attributes(data)$class != "STRAND Data Object"){
        stop("fit_NBDA_model() requires a named list of data objects of class: STRAND Data Object. Please use make_strand_data() to build the elements of your data list.")
    }

    if(!("NBDA" %in% attributes(data)$supported_models)){
        stop("The supplied data are not appropriate for an NBDA. Please ensure that an appropriate data list is provided.")
    }
    
    ############### Priors
    data$export_network = 0

    if(is.null(priors)){
      data$priors =  make_priors()
      } else{
    data$priors = priors
      }


    ############################################################################# Fit model
    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","network_based_diffusion_analysis.stan"))

     data$individual_predictors = NULL
     data$dyadic_predictors = NULL
     data$block_predictors = NULL

    if(mode=="mcmc"){
      fit = model$sample(
        data = unclass(data),
        seed = stan_mcmc_parameters$seed,
        chains = stan_mcmc_parameters$chain,
        parallel_chains = stan_mcmc_parameters$parallel_chains,
        refresh = stan_mcmc_parameters$refresh,
        iter_warmup = stan_mcmc_parameters$iter_warmup,
        iter_sampling = stan_mcmc_parameters$iter_sampling,
        max_treedepth = stan_mcmc_parameters$max_treedepth,
        adapt_delta = stan_mcmc_parameters$adapt_delta,
        init = stan_mcmc_parameters$init
        )
       }

    if(mode=="vb"){
     print("Variational inference is fast, but not always dependable. We recommend using vb only for test runs.")   
     fit = model$pathfinder(data = unclass(data))
     }

    if(mode=="optim"){
     print("Optimazation is fast, but not always dependable. We recommend using optim only for test runs.") 
     fit = model$optimize(data = unclass(data))
     }

    if(! mode %in% c("mcmc", "vb", "optim") ){
     stop("Must supply a legal mode value: mcmc, vb, or optim.")
    }

    bob = list(data=data, fit=fit, return_predicted_network = NA )
    attr(bob, "class") = "STRAND Model Object"
    attr(bob, "fit_type") = mode
    attr(bob, "model_type") = "NBDA"
    
    return(bob)
}
