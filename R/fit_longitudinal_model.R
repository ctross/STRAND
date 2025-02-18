#' A function to run a longitudinal combined stochastic block and social relations models using the STRAND framework
#' 
#' This function allows users to analyse empirical or simulated data using a Bayesian stochastic block and social relations model in Stan. The user must supply a list of STRAND data objects,
#' and a series of formulas following standard lm() style syntax. 
#'
#' It is important to note that all individual block (or group) assignment must be supplied as data.  Latent groups will be supported in future releases of STRAND.
#'
#' @param long_data A list of data objects of class STRAND prepared using the make_strand_data() function. The data objects must include all covariates used in the formulas listed below.
#' @param block_regression A formula for the block-level predictors. This should be specified as in lm(), e.g.: ~ Ethnicity + Sex. Dont use interactions, however.
#' @param focal_regression A formula for the predictors of out-degree (i.e., focal effects, or the effects of individual covariates on outgoing ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param target_regression A formula for the predictors of in-degree (i.e., target effects, or the effects of individual covariates on incoming ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param dyad_regression A formula for the predictors of dyadic relationships. This should be specified as in lm(), e.g.: ~ Kinship + Friendship
#' @param gaussian_error_priors Prior estimate for the measurement error in a Gaussian model. Base settings assume low very error. Should be a 2-vector for the mean and standard deviation of a zero-truncateed normal. The 2-vector will be cloned to each layer. To supply layer-specific values, supply a list of 2-vectors. One vector per layer.
#' @param coefficient_mode Set to "fixed" if all time-points should have the same coeffcients. Set to "varying" if each time-point should have its own coefficient.
#' @param random_effects_mode Set to "fixed" if all time-points should have the same coeffcients. Set to "varying" if each time-point should have its own coefficient. Setting to "varying essentially" reduces the longitudinal model to a multiplex model.
#' @param mode A string giving the mode that stan should use to fit the model. "mcmc" is default and recommended, and STRAND has functions to make processing the mcmc samples easier. Other options are "optim", to
#' use the optimizer provided by Stan, and "vb" to run the variational inference routine provided by Stan. "optim" and "vb" are fast and can be used for test runs. To process their output, however,
#' users must be familar with [cmdstanr](https://mc-stan.org/users/interfaces/cmdstan). We recommmend that users refer to the [Stan user manual](https://mc-stan.org/users/documentation/) for more information about the different modes that Stan can use. 
#' @param bandage_penalty A parameter that controls how tightly stiched together correlation structure parameters are. Default is 0.01 for tight stiching of parameters that should be equal. Relaxing to 0.05, or 0.1 can sometimes aid model performance. Setting "bandage_penalty=-1" deploys a different Stan model, which fixes the dyadic matrix perfectly. You must set init=0 below, for this model to initialize. 
#' @param eta Prior on LKJ Cholesky factor, if bandage_penalty = -1.
#' @param stan_mcmc_parameters A list of Stan parameters that often need to be tuned. Defaults set to: list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
#' We recommend 1000 sampling and warmup iterations on a single chain for exploratory model fitting. For final runs, we recommend running 2 to 4 chains for twice as long. Be sure to check r_hat, effective sample size, and traceplots.
#' @param priors A labeled list of priors for the model. User are only permitted to edit the values. Distributions are fixed. 
#' @return A STRAND model object containing the data used, and the Stan results.
#' @export
#' @examples
#' \dontrun{
#' fit = fit_longitudinal_model(data=model_dat,
#'                              block_regression = ~ Ethnicity,
#'                              focal_regression = ~ Age * NoFood,
#'                              target_regression = ~ Age * NoFood,
#'                              dyad_regression = ~ Relatedness + Friends * SameSex,
#'                              coefficient_mode="varying",
#'                              random_effects_mode="fixed",
#'                              mode="mcmc",
#'                              stan_mcmc_parameters = list(seed = 1, chains = 1, 
#'                                parallel_chains = 1, refresh = 1, 
#'                                iter_warmup = 100, iter_sampling = 100,
#'                                max_treedepth = NULL, adapt_delta = NULL)
#'                               )
#' }
#' 

fit_longitudinal_model = function(long_data,
                                  block_regression,
                                  focal_regression,
                                  target_regression,
                                  dyad_regression,
                                  gaussian_error_priors = c(0, 0.5),
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  bandage_penalty = 0.01,
                                  eta = 4,
                                  stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL,
                                                              iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL, init = NULL),
                                  priors=NULL
                                  ){
    ############################################################################ Build data
    if(is.null(names(long_data))) stop("long_data must be a named list. Please add names for each time-step. e.g., names(long_data)=paste('Time', 1:T)")
    data = make_longitudinal_data(long_data = long_data,
                                  block_regression = block_regression,
                                  focal_regression = focal_regression,
                                  target_regression = target_regression,
                                  dyad_regression = dyad_regression
                                  )

    ############################################################################# Check inputs
    if(attributes(data)$class != "STRAND Data Object"){
        stop("fit_longitudinal_model() requires a named list of data objects of class: STRAND Data Object. Please use make_strand_data() to build the elements of your data list.")
    }

    if(!("Longitudinal" %in% attributes(data)$supported_models)){
        stop("The supplied data are not appropriate for a longitudinal block plus social relations model. Please ensure that an appropriate data list is provided.")
    }

    if(!coefficient_mode %in% c("fixed", "varying")){
        stop("coefficient_mode must be set to 'fixed' or 'varying'.")
    }

    if(!random_effects_mode %in% c("fixed", "varying")){
        stop("random_effects_mode must be set to 'fixed' or 'varying'.")
    }
    
    ############### Priors
    data$export_network = 0

    if(is.null(priors)){
      data$priors =  make_priors()
      } else{
    data$priors = priors
      }

    data$bandage_penalty = bandage_penalty

    if(random_effects_mode == "fixed"){
        data$random_effects_mode = 2
    }

    if(random_effects_mode == "varying"){
        data$random_effects_mode = 1
    }

    if(coefficient_mode == "fixed"){
        data$coefficient_mode = 2
    }

    if(coefficient_mode == "varying"){
        data$coefficient_mode = 1
    }

    prior_error_mu = rep(NA, data$N_responses)
    prior_error_sigma = rep(NA, data$N_responses)

    if(is.vector(gaussian_error_priors)){
        for(l in 1:data$N_responses){
           prior_error_mu[l] = gaussian_error_priors[1]             
           prior_error_sigma[l] = gaussian_error_priors[2]   
           }          
    }

    if(is.list(gaussian_error_priors)){
        for(l in 1:data$N_responses){
           prior_error_mu[l] = gaussian_error_priors[[l]][1]             
           prior_error_sigma[l] = gaussian_error_priors[[l]][2]   
           }          
    }

    ############################################################################# Fit model
    if(bandage_penalty == -1){
      data = build_multiplex_bindings_dr_multiplex(data)
      data = build_multiplex_bindings_dr_longitudinal(data)
      data = build_multiplex_bindings_sr_longitudinal(data)
      data$eta = eta
      model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","block_plus_social_relations_model_longitudinal_pinkney.stan"))
        } else{
      model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","block_plus_social_relations_model_longitudinal.stan"))      
        }

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
    attr(bob, "model_type") = "Longitudinal"
    attr(bob, "random_effects_mode") = random_effects_mode
    attr(bob, "coefficient_mode") = coefficient_mode
    
    return(bob)
}
