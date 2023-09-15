#' A function to run combined stochastic block and social relations models using the STRAND framework
#' 
#' This function allows users to analyse empirical or simulated data using a Bayesian stochastic block and social relations model in Stan. The user must supply a STRAND data object,
#' and a series of formulas following standard lm() style syntax. 
#'
#' It is important to note that all individuals' block (or group) assignment must be supplied as data.  Latent groups will be supported in future releases of STRAND.
#'
#' @param 
#' data A data object of class STRAND, prepared using the make_strand_data() function. The data object must include all covariates used in the formulas listed below.
#' @param 
#' block_regression A formula for the block-level predictors. This should be specified as in lm(), e.g.: ~ Ethnicity + Sex. Dont use interactions, however.
#' @param 
#' focal_regression A formula for the predictors of out-degree (i.e., focal effects, or the effects of individual covariates on outgoing ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param 
#' target_regression A formula for the predictors of in-degree (i.e., target effects, or the effects of individual covariates on incoming ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param 
#' dyad_regression A formula for the predictors of dyadic relationships. This should be specified as in lm(), e.g.: ~ Kinship + Friendship
#' @param 
#' mode A string giving the mode that stan should use to fit the model. "mcmc" is default and recommended, and STRAND has functions to make processing the mcmc samples easier. Other options are "optim", to
#' use the optimizer provided by Stan, and "vb" to run the variational inference routine provided by Stan. "optim" and "vb" are fast and can be used for test runs. To process their output, however,
#' users must be familar with [cmdstanr](https://mc-stan.org/users/interfaces/cmdstan). We recommmend that users refer to the [Stan user manual](https://mc-stan.org/users/documentation/) for more information about the different modes that Stan can use. 
#' @param 
#' return_predicted_network Should predicted tie probabilities be returned? Requires large memory overhead, but can be used to check model fit.
#' @param 
#' stan_mcmc_parameters A list of Stan parameters that often need to be tuned. Defaults set to: list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
#' We recommend 1000 sampling and warmup iterations on a single chain for exploratory model fitting. For final runs, we recommend running 2 to 4 chains for twice as long. Be sure to check r_hat, effective sample size, and traceplots.
#' @param 
#' priors A labeled list of priors for the model. User are only permitted to edit the values. Distributions are fixed. 
#' @return A STRAND model object containing the data used, and the Stan results.
#' @export
#' @examples
#' \dontrun{
#' fit = fit_block_plus_social_relations_model(data=model_dat,
#'                                             block_regression = ~ Ethnicity,
#'                                             focal_regression = ~ Age * NoFood,
#'                                             target_regression = ~ Age * NoFood,
#'                                             dyad_regression = ~ Relatedness + Friends * SameSex,
#'                                             mode="mcmc",
#'                                             stan_mcmc_parameters = list(seed = 1, chains = 1, 
#'                                               parallel_chains = 1, refresh = 1, 
#'                                               iter_warmup = 100, iter_sampling = 100,
#'                                               max_treedepth = NULL, adapt_delta = NULL)
#'                                              )
#' }
#' 

fit_block_plus_social_relations_model = function(data,
                                    block_regression,
                                    focal_regression,
                                    target_regression,
                                    dyad_regression,
                                    mode="mcmc",
                                    return_predicted_network=FALSE,
                                    stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL,
                                                                iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL),
                                    priors=NULL
                                    ){

    ############################################################################# Check inputs
    if(attributes(data)$class != "STRAND Data Object"){
        stop("fit_block_plus_social_relations_model() requires a data object of class: STRAND Data Object. Please use make_strand_data() to build your data list.")
    }

    if(!("SRM+SBM" %in% attributes(data)$supported_models)){
        stop("The supplied data are not appropriate for a block plus social relations model. Please ensure that self_report data are single sampled and a group variable is provided.")
    }

    if(data$N_individual_predictors==0 & focal_regression != ~ 1){
        stop("No individual covariate data has been provided. focal_regression must equal ~ 1 ")
    }

    if(data$N_individual_predictors==0 & target_regression != ~ 1){
        stop("No individual covariate data has been provided. target_regression must equal ~ 1 ")
    }

    if(data$N_dyadic_predictors==0 & dyad_regression != ~ 1){
        stop("No individual covariate data has been provided. dyad_regression must equal ~ 1 ")
    }

    if(data$N_block_predictors==0 & block_regression != ~ 1){
        stop("No block covariate data has been provided. block_regression must equal ~ 1 ")
    }
    
    ############################################################################# Prepare data and parse formulas
     ind_names = colnames(data$individual_predictors)
     dyad_names = names(data$dyadic_predictors)

     ################################################################ Dyad model matrix
     if(data$N_dyadic_predictors>0){
     dyad_dims = c(data$N_id, data$N_id, length(dyad_names))

     dyad_dat = list()
     for(i in 1:dyad_dims[3]){
      dyad_dat[[i]] = c(data$dyadic_predictors[[i]])  
     }

     #dyad_dat = do.call(rbind.data.frame, dyad_dat)
     dyad_dat = as.data.frame(do.call(cbind, dyad_dat))
     colnames(dyad_dat) = dyad_names
     dyad_model_matrix = model.matrix( dyad_regression , dyad_dat )

     dyad_dat_out = array(NA, c(dyad_dims[1], dyad_dims[2], ncol(dyad_model_matrix)))
     for(i in 1:ncol(dyad_model_matrix)){
      dyad_dat_out[,,i] = matrix(dyad_model_matrix[,i], nrow=dyad_dims[1], ncol=dyad_dims[2])  
     }

     dimnames(dyad_dat_out)[[3]] = colnames(dyad_model_matrix)
     data$dyad_set = dyad_dat_out
     } else{
      data$dyad_set = array(1, c(data$N_id, data$N_id, 1))
     }

     ################################################################ Individual model matrix
     if(data$N_individual_predictors>0){
      data$focal_set = model.matrix( focal_regression , data$individual_predictors )
      data$target_set = model.matrix( target_regression , data$individual_predictors )
     } else{
      data$focal_set = matrix(1,nrow=data$N_id, ncol=1)
      data$target_set = matrix(1,nrow=data$N_id, ncol=1)
     }
    
    data$N_params = c(ncol(data$focal_set), ncol(data$target_set), dim(data$dyad_set)[3])

    ################################################################ Block model matrix
     if(data$N_block_predictors>0){
      data$block_set = model.matrix( block_regression , data$block_predictors )
     } else{
      data$block_set = as.array(matrix(1, nrow=data$N_id, ncol=1))
     }

     data$N_group_vars = ncol(data$block_set) 
     data$N_groups_per_var = rep(NA, data$N_group_vars)

     for(i in 1:data$N_group_vars){
      data$N_groups_per_var[i] = length(unique(data$block_set[,i]))  
     }

     data$N_groups_per_var = as.array(data$N_groups_per_var)
  
     data$max_N_groups = max(data$N_groups_per_var)

    ############### Priors
    data$export_network = ifelse(return_predicted_network==TRUE, 1, 0)

    if(is.null(priors)){
      data$priors =  make_priors()
      } else{
    data$priors = priors
      }

    ############################################################################# Fit model

    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","block_plus_social_relations_model.stan"))

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
        adapt_delta = stan_mcmc_parameters$adapt_delta
        )
       }

    if(mode=="vb"){
     print("Variational inference is fast, but not always dependable. We recommend using vb only for test runs.")   
     fit = model$variational(data = unclass(data), seed = 123, output_samples = 2000)
     }

    if(mode=="optim"){
     print("Optimazation is fast, but not always dependable. We recommend using optim only for test runs.") 
     fit = model$optimize(data = unclass(data), seed = 123)
     }

    if(! mode %in% c("mcmc", "vb", "optim") ){
     stop("Must supply a legal mode value: mcmc, vb, or optim.")
    }

    bob = list(data=data, fit=fit, return_predicted_network=return_predicted_network )
    attr(bob, "class") = "STRAND Model Object"
    attr(bob, "fit_type") = mode
    attr(bob, "model_type") = "SRM+SBM"
    
    return(bob)
}

