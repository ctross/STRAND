#' A function to simulate longitudinal directed networks using a combined stochastic block and social relations model
#'
#' This is a function to simulate multilayer longitudinal network data with a stochastic block structure, sender-receiver effects, and dyadic reciprocity. This function
#' is essentially the union of a social relations model and a stochastic block model.
#'
#' @param N_id Number of individuals.
#' @param N_timesteps Number of network layers, one per timestep.
#' @param B List of List of matrices that hold intercept and offset terms. Log-odds. The first matrix should be 1 x 1 with the value being the intercept term. The first list is over layers/time-steps. The second is over block variables within layers.
#' @param V Number of blocking variables in each layer.
#' @param groups A list of dataframes of the block IDs of each individual for each variable in B.
#' @param sr_mu A vector for sender and receivier random effects. In most cases, this should be a vector of 0s 2*N_timesteps long.
#' @param dr_mu A vector for dyadic random effects. In most cases, this should be a vector of 0s N_timesteps long.
#' @param sr_sigma A standard deviation vector for sender and receivier random effects. The first N_timesteps elements control node-level variation in out-degree, the second N_timesteps elements control in-degree.
#' @param dr_sigma Standard deviation for dyadic random effects. This should be a vector N_timesteps long.
#' @param sr_Rho Correlation of sender-receiver effects (i.e., generalized reciprocity). Needs to be a valid correlation matrix, 2*N_timesteps by 2*N_timesteps. 
#' @param dr_Rho Correlation of dyad effects (i.e., dyadic reciprocity). Needs to be a valid correlation matrix, 2*N_timesteps by 2*N_timesteps. 
#' @param outcome_mode Outcome mode: can be "bernoulli", "poisson", or "binomial".
#' @param link_mode Link mode: can be "logit", "probit", or "log". For pois, you must use log.
#' @param individual_predictors An N_timesteps list of N_id x N_individual_parameters array of covariates.
#' @param dyadic_predictors An N_timesteps list of N_id x N_id x N_dyadic_parameters array of covariates.
#' @param individual_effects A list of 2 by N_individual_parameters matrix of slopes. The list runs over layers/timesteps. In each element, the first row gives effects of focal characteristics (on out-degree). 
#' The second row gives effects of target characteristics (on in-degree).
#' @param dyadic_effects A list of N_dyadic_parameters vectors of slopes.
#' @return A list of objects including: network (an adjacency tensor of binary outcomes), tie_strength (an adjacency tensor with probability weights), 
#' group_ids (a vector of length N_id, giving the group of each individual), individual_predictors (the supplied covariate data is saved along with the network data), 
#' and dyadic_predictors (the supplied covariate data is saved along with the network data).
#' @export
                                                                                               
simulate_longitudinal_network = function(N_id = 99,                        # Number of respondents
                                         N_timesteps = 3,                  # Number of layers slash timesteps
                                         B = NULL,                         # Tie probabilities
                                         V = 3,                            # Blocking variables
                                         groups=NULL,                      # Group IDs
                                         sr_mu = c(0,0),                   # Average sender and reciever effect log odds
                                         sr_sigma,                         # Sender and reciever effect variances 
                                         sr_Rho,                           # Correlation of sender and reciever effects
                                         dr_mu = 0,                        # Average i to j dyad effect and j to i dyad effect log odds
                                         dr_sigma,                         # Variance of dyad effects 
                                         dr_Rho,                           # Correlation of i to j dyad effect and j to i dyad effect 
                                         outcome_mode="bernoulli",         # outcome mode
                                         link_mode = "logit",              # link mode
                                         individual_predictors = NULL,     # An N_timesteps list of N x N_individual_parameters matrices of covariates
                                         dyadic_predictors = NULL,         # An N_timesteps list of N x N x N_dyadic_parameters array of covariates
                                         individual_effects = NULL,        # The effects of predictors on sender effects (row 1) and receiver effects (row 2)
                                         dyadic_effects = NULL             # The effects of predictors on dyadic ties
                                       ){
# Varying effects on individuals
 sr = array(NA, c(N_timesteps, N_id, 2))

for(i in 1:N_id){
  sr_scrap_full = rmvnorm2(1 , Mu=sr_mu, sigma=sr_sigma, Rho=sr_Rho ) 

 for(l in 1:N_timesteps){
  sr_scrap = sr_scrap_full[c(l, N_timesteps+l)]

  sr[l,i,1] = sr_scrap[1]
  sr[l,i,2] = sr_scrap[2] 

 if(!is.null(individual_predictors)){
  sr[l,i,1] = sr[l,i,1] + sum(individual_effects[[l]][1,]*individual_predictors[[l]][i,]) 
  sr[l,i,2] = sr[l,i,2] + sum(individual_effects[[l]][2,]*individual_predictors[[l]][i,]) 
   }
  }
}

# Build true network
 samps = dr = p = y_true = array(NA, c(N_timesteps, N_id, N_id))
 for(l in 1:N_timesteps){
 samps[l,,] = matrix(rpois(N_id^2,50),nrow=N_id,ncol=N_id)
 }

# Loop over upper triangle and create ties from i to j, and j to i
for ( i in 1:(N_id-1) ){
 for ( j in (i+1):N_id){
# Dyadic effects
 dr_scrap_full = rmvnorm2(1, Mu=rep(dr_mu,2), sigma=rep(dr_sigma,2), Rho=dr_Rho)

 # For each layer
 for(l in 1:N_timesteps){
    dr_scrap = dr_scrap_full[c(l, N_timesteps+l)]

 if(!is.null(dyadic_predictors[[l]])){
  dr_scrap[1] = dr_scrap[1] + sum(dyadic_effects[[l]]*dyadic_predictors[[l]][i,j,]) 
  dr_scrap[2] = dr_scrap[2] + sum(dyadic_effects[[l]]*dyadic_predictors[[l]][j,i,]) 
  }

 B_i_j = B_j_i = c()
  for(v in 1:V){
    B_i_j[v] =  B[[l]][[v]][groups[[l]][i,v] , groups[[l]][j,v] ]
    B_j_i[v] =  B[[l]][[v]][groups[[l]][j,v] , groups[[l]][i,v] ]
  }

 dr[l,i,j] = dr_scrap[1] + sum(B_i_j)
 dr[l,j,i] = dr_scrap[2] + sum(B_j_i)

# Simulate outcomes
if(outcome_mode=="bernoulli"){
  if(link_mode=="logit"){
 p[l,i,j] = inv_logit( sr[l,i,1] + sr[l,j,2] + dr[l,i,j])
 y_true[l,i,j] = rbern( 1 , p[l,i,j] )

 p[l,j,i] = inv_logit( sr[l,j,1] + sr[l,i,2] + dr[l,j,i])
 y_true[l,j,i] = rbern( 1 , p[l,j,i] )
  }

  if(link_mode=="probit"){
 p[l,i,j] = pnorm( sr[l,i,1] + sr[l,j,2] + dr[l,i,j])
 y_true[l,i,j] = rbern( 1 , p[l,i,j] )

 p[l,j,i] = pnorm( sr[l,j,1] + sr[l,i,2] + dr[l,j,i])
 y_true[l,j,i] = rbern( 1 , p[l,j,i] )
  }
 }

 if(outcome_mode=="binomial"){
  if(link_mode=="logit"){
 p[l,i,j] = inv_logit( sr[l,i,1] + sr[l,j,2] + dr[l,i,j])
 y_true[l,i,j] = rbinom( 1 , size=samps[l,i,j], prob=p[l,i,j] )

 p[l,j,i] = inv_logit( sr[l,j,1] + sr[l,i,2] + dr[l,j,i])
 y_true[l,j,i] = rbinom( 1 , size=samps[l,j,i], prob=p[l,j,i] )
  }

 if(link_mode=="probit"){
 p[l,i,j] = pnorm( sr[l,i,1] + sr[l,j,2] + dr[l,i,j])
 y_true[l,i,j] = rbinom( 1 , size=samps[l,i,j], prob=p[l,i,j] )

 p[l,j,i] = pnorm( sr[l,j,1] + sr[l,i,2] + dr[l,j,i])
 y_true[l,j,i] = rbinom( 1 , size=samps[l,j,i], prob=p[l,j,i] )
  }
 }

 if(outcome_mode=="poisson"){
  if(link_mode=="log"){
 p[l,i,j] = exp(sr[l,i,1] + sr[l,j,2] + dr[l,i,j])
 y_true[l,i,j] = rpois( 1 , p[l,i,j] )

 p[l,j,i] = exp(sr[l,j,1] + sr[l,i,2] + dr[l,j,i])
 y_true[l,j,i] = rpois( 1 , p[l,j,i] )
 }
 }

 }
 }}

 for(l in 1:N_timesteps){
 for(i in 1:N_id){
    y_true[l,i,i] = 0
    p[l,i,i] = 0
    dr[l,i,i] = 0
 }}

return(list(network=y_true, tie_strength=p,  group_ids=groups, individual_predictors=individual_predictors, dyadic_predictors=dyadic_predictors, sr=sr, dr=dr, exposure=samps))
}

