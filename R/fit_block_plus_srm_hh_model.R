#' A function to run combined stochastic block and social relations models using the STRAND framework with household effects
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
#' hh_focal_regression A formula for the household-level predictors of out-degree (i.e., focal effects, or the effects of individual covariates on outgoing ties). This should be specified as in lm(), e.g.: ~ Wealth 
#' @param 
#' hh_target_regression A formula for the household-level predictors of in-degree (i.e., target effects, or the effects of individual covariates on incoming ties). This should be specified as in lm(), e.g.: ~ Wealth * Residents
#' @param 
#' hh_within_regression A formula for the household-level predictors of within HH relationships. This should be specified as in lm(), e.g.: ~ Wealth
#' @param 
#' hh_between_regression A formula for the household-level predictors of between HH relationships. This should be specified as in lm(), e.g.: ~ Kinship + Friendship
#' @param 
#' mode A string giving the mode that stan should use to fit the model. "mcmc" is default and recommended, and STRAND has functions to make processing the mcmc samples easier. Other options are "optim", to
#' use the optimizer provided by Stan, and "vb" to run the variational inference routine provided by Stan. "optim" and "vb" are fast and can be used for test runs. To process their output, however,
#' users must be familar with [cmdstanr](https://mc-stan.org/users/interfaces/cmdstan). We recommmend that users refer to the [Stan user manual](https://mc-stan.org/users/documentation/) for more information about the different modes that Stan can use. 
#' @param 
#' model_version Should the slower unit-level random effect model be run? Or should STRAND use the faster bivariate Bernoulli outcome model? Options are "ulre" or "fast_bb"
#' @param 
#' return_predicted_network Should predicted tie probabilities be returned? Requires large memory overhead, but can be used to check model fit.
#' @param 
#' stan_mcmc_parameters A list of Stan parameters that often need to be tuned. Defaults set to: list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
#' @param 
#' priors A labeled list of priors for the model. User are only permitted to edit the values. Distributions are fixed. 
#' @return A STRAND model object containing the data used, and the Stan results.
#' @export
#' @examples
#' \dontrun{
#' fit = fit_block_plus_social_relations_hh_model(data=model_dat,
#'                                                block_regression = ~ 1,
#'                                                focal_regression = ~ Age,
#'                                                target_regression = ~ Age,
#'                                                dyad_regression = ~ Relatedness + Friends,
#'                                                hh_focal_regression = ~ Wealth,
#'                                                hh_target_regression = ~ Wealth,
#'                                                hh_within_regression = ~ Wealth,
#'                                                hh_between_regression = ~ Distance,
#'                                                mode="mcmc",
#'                                                stan_mcmc_parameters = list(seed = 1, chains = 1, 
#'                                                parallel_chains = 1, refresh = 1, 
#'                                                iter_warmup = 100, iter_sampling = 100,
#'                                                max_treedepth = NULL, adapt_delta = NULL)
#'                                                 )
#' }
#' 

fit_block_plus_social_relations_hh_model = function(data,
                                    block_regression,
                                    focal_regression,
                                    target_regression,
                                    dyad_regression,
                                    hh_focal_regression,
                                    hh_target_regression,
                                    hh_within_regression,
                                    hh_between_regression,
                                    mode="mcmc",
                                    engine="rstan",
                                    model_version="ulre",
                                    return_predicted_network=FALSE,
                                    stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL,
                                                                iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL),
                                    priors=NULL
                                    ){

    ############################################################################# Check inputs
    if(attributes(data)$class != "STRAND Data Object"){
        stop("fit_block_plus_social_relations_hh_model() requires a data object of class: STRAND Data Object. Please use make_strand_data() to build your data list.")
    }

    if(!("HH_SRM+SBM" %in% attributes(data)$supported_models)){
        stop("The supplied data are not appropriate for a household-level block plus social relations model. ")
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
    
     if(data$N_individual_predictors_hh==0 & hh_focal_regression != ~ 1){
        stop("No household covariate data has been provided. hh_focal_regression must equal ~ 1 ")
    }

    if(data$N_individual_predictors_hh==0 & hh_target_regression != ~ 1){
        stop("No household covariate data has been provided. hh_target_regression must equal ~ 1 ")
    }

    if(data$N_individual_predictors_hh==0 & hh_within_regression != ~ 1){
        stop("No household covariate data has been provided. hh_within_regression must equal ~ 1 ")
    }

    if(data$N_dyadic_predictors_hh==0 & hh_between_regression != ~ 1){
        stop("No household covariate data has been provided. hh_between_regression must equal ~ 1 ")
    }

    if(model_version=="fast_bb" & data$outcome_mode !=  1){
        stop("The faster bivariate Bernoulli model is only avaible for Bernoulli outcome models.")
    }

    if(model_version=="no_dr" & data$outcome_mode !=  1){
        stop("The faster no dyadic reciprocity model is only avaible for Bernoulli outcome models.")
    }


    ############################################################################# Prepare data and parse formulas
     ind_names = colnames(data$individual_predictors)
     dyad_names = names(data$dyadic_predictors)

     ind_names_hh = colnames(data$hh_individual_predictors)
     dyad_names_hh = names(data$hh_dyadic_predictors)

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

     ################################################################ HH Dyad model matrix
     if(data$N_dyadic_predictors_hh>0){
     dyad_dims_hh = c(data$N_hh, data$N_hh, length(dyad_names_hh))

     dyad_dat_hh = list()
     for(i in 1:dyad_dims_hh[3]){
      dyad_dat_hh[[i]] = c(data$hh_dyadic_predictors[[i]])  
     }

     dyad_dat_hh = as.data.frame(do.call(cbind, dyad_dat_hh))
     colnames(dyad_dat_hh) = dyad_names_hh
     dyad_model_matrix_hh = model.matrix( hh_between_regression , dyad_dat_hh )

     dyad_dat_out_hh = array(NA, c(dyad_dims_hh[1], dyad_dims_hh[2], ncol(dyad_model_matrix_hh)))
     for(i in 1:ncol(dyad_model_matrix_hh)){
      dyad_dat_out_hh[,,i] = matrix(dyad_model_matrix_hh[,i], nrow=dyad_dims_hh[1], ncol=dyad_dims_hh[2])  
     }

     dimnames(dyad_dat_out_hh)[[3]] = colnames(dyad_model_matrix_hh)
     data$hh_between_set = dyad_dat_out_hh
     } else{
      data$hh_between_set = array(1, c(data$N_hh, data$N_hh, 1))
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

    ################################################################ Individual HH model matrix
     if(data$N_individual_predictors_hh>0){
      data$hh_focal_set = model.matrix( hh_focal_regression , data$hh_individual_predictors )
      data$hh_target_set = model.matrix( hh_target_regression , data$hh_individual_predictors )
      data$hh_within_set = model.matrix( hh_within_regression , data$hh_individual_predictors )
     } else{
      data$hh_focal_set = matrix(1,nrow=data$N_hh, ncol=1)
      data$hh_target_set = matrix(1,nrow=data$N_hh, ncol=1)
      data$hh_within_set = matrix(1,nrow=data$N_hh, ncol=1)
     }
    
    data$N_params_hh = c(ncol(data$hh_focal_set), ncol(data$hh_target_set), ncol(data$hh_within_set), dim(data$hh_between_set)[3])

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
   if(model_version=="ulre"){
     if(engine=="cmdstanr"){
    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","block_plus_social_relations_model_hh.stan"))
      }
     if(engine=="rstan"){
    model = rstan::stan_model(file=paste0(path.package("STRAND"),"/","block_plus_social_relations_model_hh.stan"))
      }
   }

   if(model_version=="fast_bb"){
    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","block_plus_social_relations_model_hh_fast.stan"))
   }

   if(model_version=="no_dr"){
    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","block_plus_social_relations_model_hh_nodr.stan"))
   }


    if(mode=="mcmc"){
        if(engine=="cmdstanr"){
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

    if(engine=="rstan"){
       fit = rstan::sampling(model, data = unclass(model_dat),
          seed = stan_mcmc_parameters$seed,
          chains = stan_mcmc_parameters$chain,
          warmup = stan_mcmc_parameters$iter_warmup,
          iter = (stan_mcmc_parameters$iter_sampling + stan_mcmc_parameters$iter_warmup),
          refresh = stan_mcmc_parameters$refresh,
          control = list(
           max_treedepth = stan_mcmc_parameters$max_treedepth,
           adapt_delta = stan_mcmc_parameters$adapt_delta
           )
         )
       }
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
    attr(bob, "model_type") = "HH_SRM+SBM"
    attr(bob, "model_version") = model_version
    attr(bob, "engine") = engine
    
    return(bob)
}
