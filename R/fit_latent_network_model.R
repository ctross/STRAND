#' A function to run latent network models using the STRAND framework
#' 
#' This function allows users to analyse empirical or simulated data using a Bayesian latent network model in Stan. The user must supply a STRAND data object,
#' and a series of formulas following standard lm() style syntax. 
#'
#' It is important to note that all individuals block (or group) assignment must be supplied as data.  Latent groups will be supported in future releases of STRAND.
#'
#' @param data A data object of class STRAND, prepared using the make_strand_data() function. The data object must include all covariates used in the formulas listed below.
#' @param fpr_regression A formula for the predictors of false positive rate. Specified as in lm(), e.g.: ~ Age + Education
#' @param rtt_regression A formula for the predictors of the recall rate of true ties. Specified as in lm(), e.g.: ~ Age + Education
#' @param theta_regression A formula for the predictors of theta, the probability that a given individual duplicates a response from layer 1 into layer 2. Specified as in lm(), e.g.: ~ 1
#' @param block_regression A formula for the block-level predictors. This should be specified as in lm(), e.g.: ~ Ethnicity + Sex. Dont use interactions, however.
#' @param focal_regression A formula for the predictors of out-degree (i.e., focal effects, or the effects of individual covariates on outgoing ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param target_regression A formula for the predictors of in-degree (i.e., target effects, or the effects of individual covariates on incoming ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param dyad_regression A formula for the predictors of dyadic relationships. This should be specified as in lm(),, e.g.: ~ Kinship + Friendship
#' @param stop_reflection_invariance If TRUE we add a penalty to the target: target += normal_lpdf(sum(p) | S, penalty). This prevents the backwards loading of p as 1-p, by forcing p to be sparse. The cost is a slower runtime.
#' @param mode A string giving the mode stan should use to fit the model. "mcmc" is default and recommended, and STRAND has functions to make processing the mcmc samples easier. Other options are "optim", to
#' use the optimizer provided by Stan, and "vb" to run the variational inference routine provided by Stan. "optim" and "vb" are fast and can be used for test runs. To process their output, however,
#' users must be familar with [cmdstanr](https://mc-stan.org/users/interfaces/cmdstan). We recommmend that users refer to the [Stan user manual](https://mc-stan.org/users/documentation/) for more information about the different modes that Stan can use. 
#' @param return_predicted_network Should predicted tie probabilities be returned? Requires large memory overhead, but can be used to check model fit.
#' @param mcmc_parameters A list of Stan parameters that often need to be tuned. Defaults set to: list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
#' We recommend 1000 sampling and warmup iterations on a single chain for exploratory model fitting. For final runs, we recommend running 2 to 4 chains for twice as long. Be sure to check r_hat, effective sample size, and traceplots.
#' @param priors A labeled list of priors for the model. Only edits of the values are permitted. Distributions are fixed. 
#' @return A STRAND model object containing the data used, and the Stan results.
#' @export
#' @examples
#' \dontrun{
#' fit = fit_latent_network_model(data=model_dat,
#'                                fpr_regression = ~ Age + Education,
#'                                rtt_regression = ~ Age + Education,
#'                                theta_regression = ~ 1,
#'                                block_regression = ~ Ethnicity,
#'                                focal_regression = ~ Age * NoFood,
#'                                target_regression = ~ Age * NoFood,
#'                                dyad_regression = ~ Relatedness + Friends * SameSex,
#'                                mode="mcmc",
#'                                mcmc_parameters = list(seed = 1, chains = 1, 
#'                                  parallel_chains = 1, 
#'                                  refresh = 1, iter_warmup = 100, iter_sampling = 100,
#'                                  max_treedepth = NULL, adapt_delta = NULL)
#'                                )
#' }
#' 

fit_latent_network_model = function(data,
                                    fpr_regression,
                                    rtt_regression,
                                    theta_regression,
                                    block_regression,
                                    focal_regression,
                                    target_regression,
                                    dyad_regression,
                                    stop_reflection_invariance = TRUE,
                                    mode="mcmc",
                                    return_predicted_network=FALSE,
                                    mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, 
                                                      iter_warmup = 500, iter_sampling = 500, 
                                                      max_treedepth = 12, adapt_delta = 0.95, 
                                                      chain_method = "vectorized", cores=1, init = 2),
                                    priors=NULL
                                    ){
    ############################################################################# Check inputs
    if(attributes(data)$class != "STRAND Data Object"){
        stop("fit_latent_network_model() requires a data object of class: STRAND Data Object. Please use make_strand_data() to build your data list.")
    }

    if(!("LNM" %in% attributes(data)$supported_models)){
        stop("The supplied data are not appropriate for a latent network model. Please ensure that self_report data are double sampled.")
    }

    if(data$N_individual_predictors==0 & focal_regression != ~ 1){
        stop("No individual covariate data has been provided. focal_regression must equal ~ 1 ")
    }

    if(data$N_individual_predictors==0 & target_regression != ~ 1){
        stop("No individual covariate data has been provided. target_regression must equal ~ 1 ")
    }

    if(data$N_individual_predictors==0 & fpr_regression != ~ 1){
        stop("No individual covariate data has been provided. fpr_regression must equal ~ 1 ")
    }

    if(data$N_individual_predictors==0 & rtt_regression != ~ 1){
        stop("No individual covariate data has been provided. rtt_regression must equal ~ 1 ")
    }

    if(data$N_individual_predictors==0 & theta_regression != ~ 1){
        stop("No individual covariate data has been provided. theta_regression must equal ~ 1 ")
    }

    if(data$N_dyadic_predictors==0 & dyad_regression != ~ 1){
        stop("No individual covariate data has been provided. dyad_regression must equal ~ 1 ")
    }

    if(data$N_block_predictors==0 & block_regression != ~ 1){
        stop("No block covariate data has been provided. block_regression must equal ~ 1 ")
    }
    
    if(!all(data$mask[,,1] == t(data$mask[,,2]))){
        stop("A censoring mask layer is only supported in latent network models if the same mask is used for both layers.")
    }
    
    if(data$outcome_mode != 1){
        stop("Latent network models must use: outcome_mode='bernoulli'. ")
    }

    if(data$link_mode != 1){
        stop("Latent network models must use: link_mode='logit'.")
    }

    if(data$outcome_mode==4){
        stop("Gaussian outcomes not supported for this model type.")
    }

    if(attributes(data)$directed == "undirected"){
        stop("You have an undirected outcome. This model is not supported.")
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
     #dyad_dat = as.data.frame(do.call(cbind, dyad_dat))
     dyad_dat = do.call(data.frame, dyad_dat)
     colnames(dyad_dat) = dyad_names
     dyad_model_matrix = model_matrix_strand(dyad_regression , dyad_dat, "dyad_regression")

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
      data$fpr_set = model_matrix_strand(fpr_regression , data$individual_predictors, "fpr_regression")
      data$rtt_set = model_matrix_strand(rtt_regression, data$individual_predictors, "rtt_regression")
      data$theta_set = model_matrix_strand(theta_regression , data$individual_predictors, "theta_regression")
      data$focal_set = model_matrix_strand(focal_regression, data$individual_predictors, "focal_regression")
      data$target_set = model_matrix_strand(target_regression , data$individual_predictors, "target_regression")
     } else{
      data$focal_set = matrix(1,nrow=data$N_id, ncol=1)
      data$target_set = matrix(1,nrow=data$N_id, ncol=1)
      data$fpr_set = matrix(1,nrow=data$N_id, ncol=1)
      data$rtt_set = matrix(1,nrow=data$N_id, ncol=1)
      data$theta_set = matrix(1,nrow=data$N_id, ncol=1)
     }

    data$N_params = c(ncol(data$focal_set), ncol(data$target_set), ncol(data$fpr_set), ncol(data$rtt_set), ncol(data$theta_set), dim(data$dyad_set)[3])

    ################################################################ Block model matrix
     if(data$N_block_predictors>0){
      data$block_set = model_matrix_strand(block_regression , data$block_predictors, "block_regression")
     } else{
      data$block_set = matrix(1, nrow=data$N_id, ncol=1)
     }

     data$N_group_vars = ncol(data$block_set) 
     data$N_groups_per_var = rep(NA, data$N_group_vars)

     for(i in 1:data$N_group_vars){
      data$N_groups_per_var[i] = length(unique(data$block_set[,i]))  
     }
  
     data$max_N_groups = max(data$N_groups_per_var)

    ############### Priors
    data$export_network = ifelse(return_predicted_network==TRUE, 1, 0)

    if(is.null(priors)){
      data$priors =  make_priors()
      } else{
    data$priors = priors
      }

    if(stop_reflection_invariance==TRUE){
     data$stop_reflection_invariance = 1
    } else{
     data$stop_reflection_invariance = 0    
    }

    ############################################################################# Fit model
    mcmc_parameters = merge_mcmc_parameters(mcmc_parameters)
    
    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","latent_network_model.stan"))

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

    bob = list(data=data, fit=fit, return_predicted_network=return_predicted_network )
    attr(bob, "class") = "STRAND Model Object"
    attr(bob, "fit_type") = mode
    attr(bob, "model_type") = "LNM"
    
    return(bob)
}
