#' A function to simulate single layer directed networks using a social relations model
#'
#' This is a function to simulate single layer network data with sender-receiver effects and dyadic reciprocity. This function
#' is essentially a social relations model.
#'
#' @param 
#' N_id Number of individuals.
#' @param 
#' B Intercept tie log-odds.
#' @param 
#' sr_mu Mean vector for sender and receivier random effects. In most cases, this should be c(0,0).
#' @param 
#' dr_mu Mean vector for dyadic random effects. In most cases, this should be c(0,0).
#' @param 
#' sr_sigma Standard deviation vector for sender and receivier random effects. The first element controls node-level variation in out-degree, the second in in-degree.
#' @param 
#' dr_sigma Standard deviation for dyadic random effects.
#' @param 
#' sr_rho Correlation of sender-receiver effects: aka. generalized reciprocity.
#' @param 
#' dr_rho Correlation of dyad effects: aka. dyadic reciprocity.
#' @param 
#' mode Outcome mode: can be "bernoulli", "poisson", or "binomial".
#' @param 
#' individual_predictors An N_id by N_individual_parameters matrix of covariates.
#' @param 
#' dyadic_predictors An N_id by N_id by N_dyadic_parameters array of covariates.
#' @param 
#' individual_effects A 2 by N_individual_parameters matrix of slopes. The first row gives effects of focal characteristics (on out-degree). 
#' The second row gives effects of target characteristics (on in-degree).
#' @param 
#' dyadic_effects An N_dyadic_parameters vector of slopes.
#' @return A list of objects including: network (an adjacency matrix of binary outcomes), tie_strength (an adjacency matrix with probability weights), 
#' individual_predictors (the supplied covariate data is saved along with the network data), 
#' and dyadic_predictors (the supplied covariate data is saved along with the network data).
#' @export
#' @examples
#' \dontrun{
#' library(igraph)
#' N_id = 100
#' A = simulate_srm_network(N_id = N_id, B=-7, 
#'                          individual_predictor=matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1), 
#'                          individual_effects=matrix(c(1.2, 1.5), ncol=1, nrow=2),
#'                          sr_sigma = c(1.4, 0.8), sr_rho = 0.5,
#'                          dr_sigma = 1.2, dr_rho = 0.8,
#'                          mode="bernoulli"
#'                          )
#' Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
#' 
#' plot(Net, edge.arrow.size =0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
#' }
#'
                                         
simulate_srm_network = function(N_id = 99,                            # Number of respondents
                                    B = -4,                           # Tie log odds
                                    sr_mu = c(0,0),                   # Average sender (cell 1) and reciever (cell 2) effect log odds
                                    sr_sigma = c(0.3, 1.5),           # Sender (cell 1) and reciever (cell 2) effect variances 
                                    sr_rho = 0.6,                     # Correlation of sender and reciever effects
                                    dr_mu = c(0,0),                   # Average i to j dyad effect (cell 1) and j to i dyad effect (cell 2) log odds
                                    dr_sigma = 1,                     # Variance of dyad effects 
                                    dr_rho = 0.7,                     # Correlation of i to j dyad effect and j to i dyad effect 
                                    mode="bernoulli",                 # outcome mode
                                    individual_predictors = NULL,     # A matrix of covariates
                                    dyadic_predictors = NULL,         # An array of covariates
                                    individual_effects = NULL,        # The effects of predictors on sender effects (row 1) and receiver effects (row 2)
                                    dyadic_effects = NULL             # The effects of predictors on dyadic ties
                                ){
   ##################################### Run some checks
   

   ######### Individual parameters
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

# Create correlation matrices (aka matrixes)
Rho_sr = Rho_dr = diag(c(1,1))
Rho_sr[1,2] = Rho_sr[2,1] = sr_rho
Rho_dr[1,2] = Rho_dr[2,1] = dr_rho

# Varying effects on individuals
sr = matrix(NA, nrow=N_id, ncol=2)
for( i in 1:N_id){
 sr[i,] = rmvnorm2(1 , Mu=sr_mu, sigma=sr_sigma, Rho=Rho_sr ) 

 if(!is.null(individual_predictors)){
  sr[i,1] = sr[i,1] + sum(individual_effects[1,]*individual_predictors[i,]) 
  sr[i,2] = sr[i,2] + sum(individual_effects[2,]*individual_predictors[i,]) 
  }
 }

# Build true network
dr = p = y_true = matrix(NA, N_id, N_id)
samps = matrix(rpois(N_id^2,15),nrow=N_id,ncol=N_id)
# Loop over upper triangle and create ties from i to j, and j to i
for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
# Dyadic effects
 dr_scrap = rmvnorm2(1, Mu=c(dr_mu), sigma=rep(dr_sigma,2), Rho=Rho_dr)

 if(!is.null(dyadic_predictors)){
  dr_scrap[1] = dr_scrap[1] + sum(dyadic_effects*dyadic_predictors[i,j,]) 
  dr_scrap[2] = dr_scrap[2] + sum(dyadic_effects*dyadic_predictors[j,i,]) 
  }

 dr[i,j] = dr_scrap[1] + B
 dr[j,i] = dr_scrap[2] + B

# Simulate outcomes
if(mode=="bernoulli"){
 p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j])
 y_true[i,j] = rbern( 1 , p[i,j] )

 p[j,i] = inv_logit( sr[j,1] + sr[i,2] + dr[j,i])
 y_true[j,i] = rbern( 1 , p[j,i] )
 }

 if(mode=="poisson"){
 p[i,j] = exp( sr[i,1] + sr[j,2] + dr[i,j])
 y_true[i,j] = rpois( 1 , p[i,j] )

 p[j,i] = exp( sr[j,1] + sr[i,2] + dr[j,i])
 y_true[j,i] = rpois( 1 , p[j,i] )
 }

 if(mode=="binomial"){
 p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j])
 y_true[i,j] = rbinom( 1 , size=samps[i,j], prob=p[i,j] )

 p[j,i] = inv_logit( sr[j,1] + sr[i,2] + dr[j,i])
 y_true[j,i] = rbinom( 1 , size=samps[j,i], prob=p[j,i] )
 }
        }
    }

for ( i in 1:N_id ){
    y_true[i,i] = 0
    p[i,i] = 0
    dr[i,i] = 0
 }

return(list(network=y_true, tie_strength=p,  individual_predictors=individual_predictors, dyadic_predictors=dyadic_predictors, sr=sr, dr=dr, samps=samps))
}



