#' A function to run a downstream analysis on a network fit with combined stochastic block and social relations models using the STRAND framework
#' 
#' @param fit A model fit with STRAND using the fit_block_plus_social_relations_model() function. This generates a latent network, from which node-level random effects are extracted. The effects estimated here (and used to predict downstream variables) are the nodal random effects remaining *after* accounting for all covariates in the model.
#' @param data The data object of class STRAND, prepared using the make_strand_data() function, and used in the fit object above. The data object must include all covariates (and outcomes) used in the formulas listed below.
#' @param outcome A string giving the name of the individual-level variable to be predicted. This must be valid variable from the individual-level slot of the data object.
#' @param exposure If the outcome is binomial, this is the exposure or sample size for each observation.
#' @param mask An indicator if the outcome was masked/censored. If mask[i]==1, then the Stan model skips over modeling the outcome for element i.
#' @param downstream_regression A formula for the predictors of the downstream outcomes (i.e., effects of individual covariates). This should be specified as in lm(), e.g.: ~ Age * Education. This should include all desired predictors except the nodal random effects from the fit object above.
#' @param nodal_effects A string: 'out', 'in', 'both', or 'none'. This determines which nodal predictors get used in the regression predicting the downstream outcome.
#' @param standardized A Boolean, should the random effects from the fit object be standardized unit normals? If TRUE, then the raw random effects are not scaled using the SD parameter.
#' @param outcome_mode A string:  'bernoulli', 'binomial', 'poisson', 'gaussian', 'beta', 'gamma', 'negative_binomial'. The outcome type for the model.
#' @param link_mode A string: 'logit', 'probit', 'log', or 'identity'. Must correspond to the outcome mode.
#' @param mode A string giving the mode that stan should use to fit the model. "mcmc" is default and recommended, and STRAND has functions to make processing the mcmc samples easier. Other options are "optim", to
#' use the optimizer provided by Stan, and "vb" to run the variational inference routine provided by Stan. "optim" and "vb" are fast and can be used for test runs. To process their output, however,
#' users must be familar with [cmdstanr](https://mc-stan.org/users/interfaces/cmdstan). We recommmend that users refer to the [Stan user manual](https://mc-stan.org/users/documentation/) for more information about the different modes that Stan can use. 
#' @param mcmc_parameters A list of Stan parameters that often need to be tuned. Defaults set to: list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
#' We recommend 1000 sampling and warmup iterations on a single chain for exploratory model fitting. For final runs, we recommend running 2 to 4 chains for twice as long. Be sure to check r_hat, effective sample size, and traceplots.
#' @param priors A labeled list of priors for the model. User are only permitted to edit the values. Distributions are fixed. 
#' @return A STRAND model object containing the data used, and the Stan results.
#' @export
#' @examples
#' \dontrun{
#' fit = fit_downstream_nodal_model(
#'    fit = fit1, 
#'    data = model_dat,
#'    outcome = "RS",
#'    downstream_regression = ~ Age + Mass,
#'    nodal_effects = "both",
#'    standardized = TRUE,
#'    outcome_mode = "gamma",
#'    link_mode = "log",
#'    mode = "mcmc",
#'    mcmc_parameters = list(
#'      chains = 1,
#'      iter_warmup = 1500 ,
#'      iter_sampling = 1500 ,
#'      max_treedepth = 13,
#'      refresh = 1,
#'      adapt_delta = 0.98)
#'  )
#' }
#'

fit_downstream_nodal_model = function(
                                fit, 
                                data,
                                outcome,
                                exposure=NULL,
                                mask=NULL,
                                downstream_regression,
                                nodal_effects = "out",
                                standardized = TRUE,
                                outcome_mode="gaussian",
                                link_mode="identity",
                                mode="mcmc",
                                mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, 
                                                      iter_warmup = 500, iter_sampling = 500, 
                                                      max_treedepth = 12, adapt_delta = 0.95, 
                                                      chain_method = "vectorized", cores=1, init = 2),
                                priors=NULL
                                    ){

    ############################################################################# Check inputs
    if(attr(fit,"model_type") != "SRM+SBM"){
        stop("The fit object is not a valid type. Please only apply fit_downstream_nodal_model to models fit with fit_block_plus_social_relations_model().")
    }

    if(data$N_individual_predictors==0 & downstream_regression != ~ 1){
        stop("No individual covariate data has been provided, downstream_regression must equal ~ 1.")
    }

    if(!is.null(outcome_mode) & (!outcome_mode %in% c("bernoulli", "binomial", "poisson", "gaussian", "beta", "gamma", "negative_binomial"))){
          stop("outcome_mode must be either bernoulli, binomial, poisson, beta, gamma, negative_binomial, or gaussian.")
         }

    if(!is.null(link_mode) & (!link_mode %in% c("logit", "probit", "log", "identity"))){
          stop("link_mode must be either logit, probit, log, or identity.")
         }

    ###################### Outcome mode details
         if(is.null(outcome_mode)) stop(" 'outcome_mode' must be declared.")
         if(is.null(link_mode)) stop(" 'link_mode' must be declared.")

         outcome_mode_numeric = NULL
         link_mode_numeric = NULL

         if(outcome_mode=="bernoulli"){
          outcome_mode_numeric = 1

          if(!link_mode %in% c("logit", "probit")){stop("If outcome_mode is 'bernoulli', you must set link_mode to 'logit' or 'probit'.")}

          if(link_mode == "logit"){
            link_mode_numeric = 1
          } 

          if(link_mode == "probit"){
            link_mode_numeric = 2
          } 
         }

         if(outcome_mode=="binomial"){
          if(is.null(exposure)){stop("If outcome is binomial, an exposure variable must be provided.")}
          outcome_mode_numeric = 2

          if(!link_mode %in% c("logit", "probit")){stop("If outcome_mode is 'binomial', you must set link_mode to 'logit' or 'probit'.")}

          if(link_mode == "logit"){
            link_mode_numeric = 1
          } 

          if(link_mode == "probit"){
            link_mode_numeric = 2
          } 
         }

         if(outcome_mode=="poisson"){
          if(link_mode != "log"){stop("If outcome_mode is 'poisson', you must set link_mode to 'log'.")}
          outcome_mode_numeric = 3
          link_mode_numeric = 3
         }

         if(outcome_mode=="negative_binomial"){
          if(link_mode != "log"){stop("If outcome_mode is 'negative_binomial', you must set link_mode to 'log'.")}
          outcome_mode_numeric = 4
          link_mode_numeric = 3
         }

         if(outcome_mode=="gaussian"){
          if(link_mode != "identity"){stop("If outcome_mode is 'gaussian', you must set link_mode to 'identity'.")}
          outcome_mode_numeric = 5
          link_mode_numeric = 4
         }


         if(outcome_mode=="beta"){
          outcome_mode_numeric = 6

          if(!link_mode %in% c("logit", "probit")){stop("If outcome_mode is 'beta', you must set link_mode to 'logit' or 'probit'.")}

          if(link_mode == "logit"){
            link_mode_numeric = 1
          } 

          if(link_mode == "probit"){
            link_mode_numeric = 2
          } 
         }

         if(outcome_mode=="gamma"){
          if(link_mode != "log"){stop("If outcome_mode is 'gamma', you must set link_mode to 'log'.")}
          outcome_mode_numeric = 7
          link_mode_numeric = 3
         }


         if(is.null(outcome_mode_numeric)) stop("outcome_mode not supported")

    ###################################################### Outcome data
     temp_outcome = data$individual_predictors[,which(colnames(data$individual_predictors)==outcome)]
      if(length(temp_outcome)==0){
        stop("Outcome variable not found in individual_predictor set of data object.")
      }
      data$outcomes = temp_outcome
      data$outcomes_real = temp_outcome

    ###################################################### Exposure data
     if(is.null(exposure)){
      data$exposure = rep(1, data$N_id)  
      } else{
       temp_exposure = data$individual_predictors[,which(colnames(data$individual_predictors)==exposure)]
       if(length(temp_exposure)==0){
        stop("Exposure variable not found in individual_predictor set of data object.")
         }
        data$exposure = temp_exposure
      } 

    ###################################################### Mask data
     if(is.null(mask)){
      data$mask = rep(0, data$N_id)  
      } else{
       temp_mask = data$individual_predictors[,which(colnames(data$individual_predictors)==mask)]
       if(length(temp_mask)==0){
        stop("Mask variable not found in individual_predictor set of data object.")
         }
        data$mask = temp_mask
      } 


    ###################################################### Create nodal predictor vectors 
    if(attr(fit, "fit_type")=="numpyro"){
        numpyro = TRUE
    }else{
        numpyro = FALSE
    }

    ################### Network model parameters
    if(numpyro==FALSE){
     stanfit = posterior::as_draws_rvars(fit$fit$draws())
     sr_sigma = posterior::draws_of(stanfit$"sr_sigma")
     sr_L = posterior::draws_of(stanfit$"sr_L") 
     sr_raw = posterior::draws_of(stanfit$"sr_raw")
    }

    if(numpyro==TRUE){
     samps = convert_posterior(fit$fit$get_samples())
     sr_sigma = samps$sr_sigma
     sr_L = samps$sr_L 
     sr_raw = aperm(samps$sr_raw, perm = c(1, 3, 2))
    }

    sr = sr_raw

    for(q in 1:dim(sr_raw)[1]){
      for(i in 1:dim(sr_raw)[2]){
        if(standardized == TRUE){
        sr[q,i,] = (sr_L[q,,] %*% sr_raw[q,i,])
         } else{
        sr[q,i,] = sr_sigma[q,] * (sr_L[q,,] %*% sr_raw[q,i,])
         }
      }    
    }

    data$sr_mus = apply(sr, 2:3, mean)
    data$sr_sds = apply(sr, 2:3, sd)

    #################################### Determine which nodal predictors get used in regression
     if(nodal_effects %in% c("out", "in", "both", "none")){
        
       if(nodal_effects=="out"){
        data$Z = c(1,0)
       }

       if(nodal_effects=="in"){
        data$Z = c(0,1)
       }

       if(nodal_effects=="both"){
        data$Z = c(1,1)
       }

       if(nodal_effects=="none"){
        data$Z = c(0,0)
       }
        } else {
            stop("nodal_effects must be: 'out', 'in', 'both', or 'none'.")
        }


    ############################################################################# Prepare data and parse formulas
     ind_names = colnames(data$individual_predictors)

     ################################################################ Individual model matrix
     if(data$N_individual_predictors>0){
      data$focal_set = model_matrix_strand(downstream_regression , data$individual_predictors, "downstream_regression")
     } else{
      data$focal_set = matrix(1, nrow=data$N_id, ncol=1)
     }
    
    data$N_params = ncol(data$focal_set)

    ############### Priors
    data$export_network = 0

    if(is.null(priors)){
      data$priors =  make_priors()
      } else{
      data$priors = priors
      }

     data$outcome_mode = outcome_mode_numeric
     data$link_mode = link_mode_numeric


    ############################################################################# Fit model
    mcmc_parameters = merge_mcmc_parameters(mcmc_parameters)
    
    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","downstream_block_plus_social_relations_model.stan"))

    data$individual_predictors = NULL
    data$dyadic_predictors = NULL
    data$block_predictors = NULL

    if(mode=="mcmc"){
      fit = model$sample(
        data = unclass(data),
        seed = mcmc_parameters$seed,
        chains = mcmc_parameters$chains,
        parallel_chains = mcmc_parameters$parallel_chains,
        refresh = mcmc_parameters$refresh,
        iter_warmup = mcmc_parameters$iter_warmup,
        iter_sampling = mcmc_parameters$iter_sampling,
        max_treedepth = mcmc_parameters$max_treedepth,
        adapt_delta = mcmc_parameters$adapt_delta,
        init = mcmc_parameters$init
        )
       }

    if(mode=="vb"){
     print("Variational inference is fast, but not always dependable. We recommend using vb only for test runs.")   
     fit = model$pathfinder(data = unclass(data))
     }

    if(mode=="optim"){
     print("Optimazation is fast, but not always dependable. We recommend using optim only for test runs.") 
     fit = model$optimize(data = unclass(data), seed = 123)
     }

    if(! mode %in% c("mcmc", "vb", "optim") ){
     stop("Must supply a legal mode value: mcmc, vb, or optim.")
    }

    bob = list(data=data, fit=fit, return_predicted_network=NULL )
    attr(bob, "class") = "STRAND Model Object"
    attr(bob, "fit_type") = mode
    attr(bob, "model_type") = "Downstream"
    
    return(bob)
}
