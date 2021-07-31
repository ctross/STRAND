#' A function to run STRAND latent network models
#' 
#' This function allows a user to supply emprical or simulated data, and analyse it using a Bayesian block plus social relations model in Stan. The user must supply a STRAND data object,
#' and a series of formulas following standard lm() style syntax. The group of each individual must be provided as data. Latent groups will be supported in later STRAND releases.
#'
#' @param 
#' data A data object of class STRAND, prepared using the make_strand_data() function. Must include the covariates used in the formulas listed below.
#' @param 
#' focal_regression A formula for the predictors of out-degree. Aka focal effects, or the effects of individual covariates on outgoing ties. Specified as in lm(), e.g.: ~ Age * Education
#' @param 
#' target_regression A formula for the predictors of in-degree. Aka target effects, or the effects of individual covariates on incoming ties. Specified as in lm(), e.g.: ~ Age * Education
#' @param 
#' dyad_regression A formula for the predictors of dyadic relationships. Specified as in lm(), e.g.: ~ Kinship + Friendship
#' @param 
#' mode A string giving the mode stan should use to fit the model. "mcmc" is recommended, and STRAND has functions to make processing the mcmc samples easier. Other options are "optim", to
#' use the optimizer provided by Stan, and "vb" to run the variational inference routine provided by Stan. "optim" and "vb" are fast and can be used for test runs. To process their output, however,
#' users must be familar with cmdstanr and its methods.
#' @param 
#' stan_mcmc_parameters A list of Stan parameters that often need to be tuned. Defaults to: list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL, iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
#' @return A STRAND model object containing the data used, and the Stan results.
#' @export
#' @examples
#' \dontrun{
#' fit = fit_block_plus_social_relations_model( data=model_dat,
#'                                              focal_regression = ~ Age * NoFood,
#'                                              target_regression = ~ Age * NoFood,
#'                                              dyad_regression = ~ Relatedness + Friends * SameSex,
#'                                              mode="mcmc",
#'                                              stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 100,
#'                                                                          iter_sampling = 100, max_treedepth = NULL, adapt_delta = NULL)
#'                                              )
#' }
#' 

fit_block_plus_social_relations_model = function(data=model_dat,
                                    focal_regression,
                                    target_regression,
                                    dyad_regression,
                                    mode="mcmc",
                                    stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = NULL,
                                                                iter_sampling = NULL, max_treedepth = NULL, adapt_delta = NULL)
                                    ){

    ############################################################################# Check inputs
    if(attributes(data)$class != "STRAND Data Object"){
        stop("fit_block_plus_social_relations_model() requires a data object of class: STRAND Data Object. Please use make_strand_data() to build your data list.")
    }

    if(!("SRM+SBM" %in% attributes(data)$supported_models)){
        stop("The supplied data are not appropriate for a block plus social relations model. Please ensure that self_report data are single sampled and a group variable is provided.")
    }
    
    ############################################################################# Prepare data and parse formulas
    ind_names = colnames(data$individual_predictors)
    dyad_names = names(data$dyadic_predictors)


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

    data$focal_set = model.matrix( focal_regression , data$individual_predictors )
    data$target_set = model.matrix( target_regression , data$individual_predictors )
    data$dyad_set = dyad_dat_out
    
    data$N_params = c(ncol(data$focal_set), ncol(data$target_set), ncol(dyad_model_matrix))

    data$export_network = 0

    ############################################################################# Fit model

    model = cmdstan_model(paste0(path.package("STRAND"),"/","block_plus_social_relations_model.stan"))

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

    bob = list(data=data, fit=fit, return_latent_network=NA )
    attr(bob, "class") = "STRAND Model Object"
    attr(bob, "fit_type") = mode
    attr(bob, "model_type") = "SRM+SBM"
    
    return(bob)
}

