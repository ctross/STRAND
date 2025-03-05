#' A function to simulate a network-based diffusion proccess using the STRAND framework
#' 
#' This function allows users to simulate data using a NBDA model. The user must supply a list of STRAND data objects, a set of parameters,
#' and a series of formulas following standard lm() style syntax. 
#'
#' @param long_data A list of data objects of class STRAND prepared using the make_strand_data() function. The data objects must include all covariates and trait diffusion data used in the formulas listed below.
#' @param individual_focal_regression A formula for the effects of focal predictors on individual learning rate. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_focal_regression A formula for the effects of focal predictors on social attention weights. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_target_regression A formula for the effects of target predictors on social attention weights. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_dyad_regression A formula for the predictors of dyadic relationships on social attention weights. This should be specified as in lm(), e.g.: ~ Kinship + Friendship.
#' @param individual_focal_parameters A formula for the effects of focal predictors on individual learning rate. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_focal_parameters A formula for the effects of focal predictors on social attention weights. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_target_parameters A formula for the effects of target predictors on social attention weights. This should be specified as in lm(), e.g.: ~ Age * Education.
#' @param social_dyad_parameters A formula for the predictors of dyadic relationships on social attention weights. This should be specified as in lm(), e.g.: ~ Kinship + Friendship.
#' @param base_rates The intercept parameters for the individual and then social learing rates on logg-odds scale.
#' @param ces_parameters A named list: alpha is the share of social influence owing to i to j ties (which implies 1-alpha is the share due to j to i ties), sigma is the elastisity of substitution, eta is returns to scale.
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

simulate_diffusion = function(long_data,
                              individual_focal_regression,
                              social_focal_regression,
                              social_target_regression,
                              social_dyad_regression,
                              individual_focal_parameters,
                              social_focal_parameters,
                              social_target_parameters,
                              social_dyad_parameters,
                              base_rates,
                              ces_parameters = list(alpha=0.95, sigma=100, eta=1)
                          ){
    ############################################################################ Build data
    if(is.null(names(long_data))) stop("long_data must be a named list. Please add names for each time-step. e.g., names(long_data)=paste('Time', 1:T)")
    data = make_longitudinal_data(long_data = long_data,
                                  block_regression = ~ 1,
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

    A = data$outcomes
    E = data$exposure
    SRI = A/E

    ############################################################################# Add parameters
    N_params = data$N_params - 1
    focal_set = data$focal_set[,-1,] 
    target_set = data$target_set[,-1,] 
    dyad_set = data$dyad_set[,,-1,]
    ind_focal_set = data$ind_focal_set[,-1,] 

    for(v in 1:N_params[1]){
     focal_set[,v,] = social_focal_parameters[v]*focal_set[,v,]
    }

    for(v in 1:N_params[2]){
     target_set[,v,] = social_target_parameters[v]*target_set[,v,]
    }

    for(v in 1:N_params[3]){
     dyad_set[,,v,] = social_dyad_parameters[v]*dyad_set[,,v,]
    }

    for(v in 1:N_params[4]){
     ind_focal_set[,v,] = individual_focal_parameters[v]*ind_focal_set[,v,]
    }

    ############################################################################# Diffusion model
    psi = rep(NA, data$N_id)
    diffusion_outcomes = data$diffusion_outcomes

    for(t in 1:data$N_responses){
     for(i in 1:data$N_id){
        # Individual learning rate
        theta = inv_logit(base_rates[1] + sum(ind_focal_set[i,,t])) 

        if(t==1){
        diffusion_outcomes[i, t] = rbinom(1, prob=theta, size=1) 
        } else{

        # Check if i has trait
        if(diffusion_outcomes[i, t-1]==1){
         diffusion_outcomes[i, t] = 1
        }  else{

        # Social learning weights for each alter
        for(j in 1:data$N_id){
         psi[j] = inv_logit(base_rates[2] + sum(focal_set[i,,t]) + sum(target_set[j,,t]) + sum(dyad_set[i,j,,t]))
        }

        # Now integrate all information
        attention_weighted_network =  psi * CES(SRI[i,,t]*diffusion_outcomes[,t-1], SRI[,i,t]*diffusion_outcomes[,t-1], ces_parameters$alpha, ces_parameters$sigma, ces_parameters$eta)
        latent_prediction = theta + (1-theta)*inv_logit_shifted( sum(attention_weighted_network) )
        diffusion_outcomes[i, t] = rbinom(1, prob=latent_prediction, size=1)
      }}                   
     }}
    
    ########## Send out
    for(t in 1:data$N_responses){
     long_dat[[t]]$diffusion_outcomes = diffusion_outcomes[,t]
    }

    return(long_dat)
}
