#' A function to simulate single layer directed networks using a combined stochastic block and social relations model
#'
#' This is a function to simulate single layer network data with a stochastic block structure, sender-receiver effects, and dyadic reciprocity. This function
#' is essentially the union of a social relations model and a stochastic block model.
#'
#' @param N_id Number of individuals.
#' @param B List of matrices that hold intercept and offset terms. Log-odds. The first matrix should be  1 x 1 with the value being the intercept term.
#' @param V Number of blocking variables in B.
#' @param groups Dataframe of the block IDs of each individual for each variable in B.
#' @param sr_mu Mean vector for sender and receivier random effects. In most cases, this should be c(0,0).
#' @param dr_mu Mean vector for dyadic random effects. In most cases, this should be c(0,0).
#' @param exposure_mu Intercept term for log-odds of encounter.
#' @param censoring_mu Intercept term for log-odds of censoring.
#' @param sr_sigma A standard deviation vector for sender and receivier random effects. The first element controls node-level variation in out-degree, the second in in-degree.
#' @param dr_sigma Standard deviation for dyadic random effects.
#' @param exposure_sigma Standard deviation for exposure random effects.
#' @param censoring_sigma Standard deviation for censoring random effects.
#' @param sr_rho Correlation of sender-receiver effects (i.e., generalized reciprocity).
#' @param dr_rho Correlation of dyad effects (i.e., dyadic reciprocity).
#' @param exposure_max Max sample size of observations for a given focal.
#' @param N_trials Number of binomial trials in follow-up detectability experiment.
#' @param outcome_mode Outcome mode: only "binomial" is supported.
#' @param link_mode Link mode: can be "logit", "probit".
#' @param individual_predictors An N_id by N_individual_parameters matrix of covariates.
#' @param dyadic_predictors An N_id by N_id by N_dyadic_parameters array of covariates.
#' @param exposure_predictors An N_id by N_individual_parameters matrix of covariates.
#' @param censoring_predictors An N_id by N_individual_parameters matrix of covariates.
#' @param individual_effects A 2 by N_individual_parameters matrix of slopes. The first row gives effects of focal characteristics (on out-degree). 
#' The second row gives effects of target characteristics (on in-degree).
#' @param dyadic_effects An N_dyadic_parameters vector of slopes.
#' @param exposure_effects An N_parameters vector of slopes.
#' @param censoring_effects An N_parameters vector of slopes.
#' @return A list of objects including: network (an adjacency matrix of binary outcomes), tie_strength (an adjacency matrix with probability weights), 
#' group_ids (a vector of length N_id, giving the group of each individual), individual_predictors (the supplied covariate data is saved along with the network data), 
#' and dyadic_predictors (the supplied covariate data is saved along with the network data).
#' @export
#'
                                                                                               
simulate_sbm_plus_srm_network_with_measurement_bias = function(N_id = 30,
                                                                                                                           
                                                               B = NULL,
                                                               V = 3,
                                                               groups=NULL,

                                                               sr_mu = c(0,0),
                                                               sr_sigma = c(0.3, 1.5),
                                                               sr_rho = 0.6,

                                                               dr_mu = 0,
                                                               dr_sigma = 1,
                                                               dr_rho = 0.7,

                                                               exposure_mu = 1.9,
                                                               exposure_sigma = 0.01,
                                                               exposure_max = 50,
                                                               
                                                               censoring_mu = 1.9,
                                                               censoring_sigma = 0.01,
                                                               N_trials = 20,

                                                               outcome_mode = "binomial",
                                                               link_mode = "logit",              

                                                               individual_predictors = NULL,
                                                               dyadic_predictors = NULL,                                                               
                                                               exposure_predictors = NULL,
                                                               censoring_predictors = NULL,

                                                               individual_effects = NULL,
                                                               dyadic_effects = NULL,
                                                               exposure_effects = NULL,
                                                               censoring_effects = NULL){  
###############################
####### Run some checks #######
###############################
   if(!is.null(individual_predictors)){
     if(is.null(individual_effects)){
      stop("If individual_predictors is supplied, a matching matrix of individual_effects must be supplied.")
     }

    if(length(individual_effects[1,]) != length(individual_predictors[1,])){
      stop("The number of columns of individual_effects must match that of individual_predictors.")
    }

    if(nrow(individual_effects) != 2 ){
      stop("The number of rows of individual_effects must be 2: one for sender effects, the other for receiver. Un-needed slopes can be set to 0.")
    }

    if(nrow(individual_predictors) != N_id ){
     stop("The number of rows of individual_predictors must be N_id.")
    }
   }

   if(outcome_mode != "binomial"){
    stop("Only binomial output supported.")
   }

###############################
####### Model true network ####
###############################
  # True networks make reference to the network without bias.
  # Create correlation matrices (aka matrixes). It defines the reciprocity of interactions.
  Rho_sr = Rho_dr = Rho_int = diag(c(1,1))
  Rho_sr[1,2] = Rho_sr[2,1] = sr_rho
  Rho_dr[1,2] = Rho_dr[2,1] = dr_rho

  # Varying effects on individuals
  # ## Determine for each dyad its baseline interaction frequencies (rmvnorm2) +
  # ## individuals' characteristics effect on interactions sum(individual_effects[1,] * individual_predictors[i,])
  sr = matrix(NA, nrow=N_id, ncol=2)
  for(i in 1:N_id){
    sr[i,] = rmvnorm2(1 , Mu=sr_mu, sigma=sr_sigma, Rho= Rho_sr)

    if(!is.null(individual_predictors)){
      sr[i,1] = sr[i,1] + sum(individual_effects[1,]*individual_predictors[i,])
      sr[i,2] = sr[i,2] + sum(individual_effects[2,]*individual_predictors[i,])
    }
  }

  # Build network
  dr = p = y = matrix(0, N_id, N_id)
  # Loop over upper triangle and create ties from i to j, and j to i
  # Determine for each interaction its baseline interaction frequencies (rmvnorm2) +
  # dyads' characteristics effect on interactions sum(dyadic_effects * dyadic_predictors[i, j,])
  for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
      # Dyadic effects
       dr_scrap = rmvnorm2(1, Mu=rep(dr_mu,2), sigma=rep(dr_sigma,2), Rho=Rho_dr)

       if(!is.null(dyadic_predictors)){
         dr_scrap[1] = dr_scrap[1] + sum(dyadic_effects*dyadic_predictors[i,j,])
         dr_scrap[2] = dr_scrap[2] + sum(dyadic_effects*dyadic_predictors[j,i,])
        }

       # ## If subgroups are declared, determine within and between group link frequencies.
       B_i_j = B_j_i = c()
       for(v in 1:V){
          B_i_j[v] =  B[[v]][groups[i,v] , groups[j,v] ]
          B_j_i[v] =  B[[v]][groups[j,v] , groups[i,v] ]
       }
       dr[i,j] = dr_scrap[1] + sum(B_i_j)
       dr[j,i] = dr_scrap[2] + sum(B_j_i)
    }
  }

  # Sum dyad and interaction weights and create tie probability matrix.
  for (i in 1:(N_id-1) ){
    for (j in (i+1):N_id){
      if(link_mode=="logit"){
      p[i,j] = inv_logit(sr[i,1] + sr[j,2] + dr[i,j])
      p[j,i] = inv_logit(sr[j,1] + sr[i,2] + dr[j,i])
      }

      if(link_mode=="probit"){
      p[i,j] = pnorm(sr[i,1] + sr[j,2] + dr[i,j])
      p[j,i] = pnorm(sr[j,1] + sr[i,2] + dr[j,i])
      }
    }
  }

##################################
#######  Model exposure bias   ###
##################################
  ideal_samps = matrix(NA, nrow=N_id, ncol=N_id) # Sample without bias.
  true_samps = matrix(NA, nrow=N_id, ncol=N_id)  # Sample with bias.
  exposure_prob = rep(NA, N_id)                  # The probability of sample based on individual characteristics.
  true_exposure = rep(NA, N_id)                  # Observation time with bias.
  exposure_offset = rep(NA, N_id)                # Observation
  exposure_factors = rep(0, N_id)               
  diag(true_samps) = 0
  diag(ideal_samps) = 0
  
  # For each individual, determine:
  for(i in 1:N_id){
        
    if(!is.null(exposure_predictors)){
    exposure_factors[i] = exposure_factors[i] + sum(exposure_effects*exposure_predictors[i,])
    }
    
    exposure_offset[i] = rnorm(1, 0, exposure_sigma)                                             # Observation bias.
    
    if(link_mode=="logit"){
    exposure_prob[i] = inv_logit(exposure_mu + exposure_factors[i] + exposure_offset[i])         # Its observation probabilities based on individual characteristics.
      }

    if(link_mode=="probit"){
    exposure_prob[i] = pnorm(exposure_mu + exposure_factors[i] + exposure_offset[i])             # Its observation probabilities based on individual characteristics.
      }

    true_exposure[i] = rbinom(1, size=exposure_max, prob=exposure_prob[i])                       # Its observations with bias.
  }

  # For each dyad, determine its sampling with and without bias according to previous steps. 
  for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
      ideal_samps[i,j] = sum(exposure_max + exposure_max)
      ideal_samps[j,i] = ideal_samps[i,j]

      true_samps[i,j] = sum(true_exposure[i] + true_exposure[j])
      true_samps[j,i] = true_samps[i,j]
    }
  }
  
###################################
#######  Model censoring bias   ###
###################################
   censoring_offset = rep(NA, N_id)                
   censoring_prob = rep(0, N_id) 
   censoring_factors = rep(0, N_id)               

   # For each individual, determine:
   for(i in 1:N_id){

    if(!is.null(censoring_predictors)){
    censoring_factors[i] = censoring_factors[i] + sum(censoring_effects*censoring_predictors[i,])
    }

    censoring_offset[i] = rnorm(1, 0, censoring_sigma) 
   if(link_mode=="logit"){
    censoring_prob[i] = inv_logit(censoring_mu + censoring_factors[i] + censoring_offset[i]) 
      }        

   if(link_mode=="probit"){
    censoring_prob[i] = pnorm(censoring_mu + censoring_factors[i] + censoring_offset[i]) 
      }                                        
    }
    
  eta = 1 - censoring_prob             # Flip to represent NOT censoring
  
###############################
#######  Model outcomes #######
###############################
  # Create an interaction matrix according to the ties probability matrix and observation bias.
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
      y[i,j] = rbinom( 1 , size=true_samps[i,j], prob = p[i,j]*eta[i]*eta[j] ) 
      y[j,i] = rbinom( 1 , size=true_samps[j,i], prob = p[j,i]*eta[j]*eta[i] )
    }
  }
  
  trials = rep(N_trials, N_id)
  detected = NULL
  for(i in 1: N_id){
    detected[i]  = rbinom(1, size = trials[i], prob = eta[i])
  }

  diag(y) = 0
  diag(p) = 0
  diag(dr) = 0
  
  colnames(y) = rownames(y) = 1:ncol(y)
  
  return(list(network= y,
              tie_strength=p,
              group_ids=groups,
              individual_predictors=individual_predictors,
              dyadic_predictors=dyadic_predictors,
              exposure_predictors=exposure_predictors,
              sr=sr,
              dr=dr,
              true_samps=true_samps,
              ideal_samps=ideal_samps,
              true_exposure=true_exposure,
              censoring_prob = censoring_prob,
              detected = detected,
              trials = trials))
}
