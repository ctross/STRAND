#' A function to simulate single layer directed networks with a stochastic block structure
#'
#' This is a function to simulate single layer network data with a stochastic block structure. 
#'
#' @param N_id Number of individuals.
#' @param B List of matrices that hold intercept and offset terms. Log-odds. The first matrix should be  1 x 1 with the value being the intercept term.
#' @param V Number of blocking variables in B.
#' @param groups Dataframe of the block IDs of each individual for each variable in B.
#' @param error_sigma Standard deviation for measurement error in Gaussian models.
#' @param outcome_mode Outcome mode: can be "bernoulli", "poisson", "binomial", or "gaussian".
#' @param link_mode Link mode: can be "logit", "probit", "log", or "identity". For poisson, you must use log. For Gaussian you must use identity.
#' @param individual_predictors An N_id by N_individual_parameters matrix of covariates.
#' @param dyadic_predictors An N_id by N_id by N_dyadic_parameters array of covariates.
#' @param individual_effects A 2 by N_individual_parameters matrix of slopes. The first row gives effects of focal characteristics (on out-degree). 
#' The second row gives effects of target characteristics (on in-degree).
#' @param dyadic_effects An N_dyadic_parameters vector of slopes.
#' @return A list of objects including: network (an adjacency matrix of binary outcomes), tie_strength (an adjacency matrix with probability weights), 
#' group_ids (a vector of length N_id, giving the group of each individual), individual_predictors (the supplied covariate data is saved along with the network data), 
#' and dyadic_predictors (the supplied covariate data is saved along with the network data).
#' @export
#' @examples
#' \dontrun{
#' library(igraph)
#' V = 1            # One blocking variable
#' G = 3            # Three categories in this variable
#' N_id = 100       # Number of people
#'
#' clique = sample(1:3, N_id, replace=TRUE)
#' B = matrix(-6, nrow=G, ncol=G)
#' diag(B) = -1.5 # Block matrix
#'
#' B[1,3] = -5
#' B[3,2] = -4
#' 
#' A = simulate_sbm_network(N_id = N_id, B=list(B=B), V=V, groups = data.frame(clique=factor(clique)),
#'                          individual_predictor=matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1), 
#'                          individual_effects=matrix(c(1.7, 0.3),ncol=1, nrow=2),
#'                          outcome_mode="bernoulli"
#'                                )
#'
#' Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
#' V(Net)$color = c("turquoise4","gray13", "goldenrod3")[A$group_ids$clique]
#' 
#' plot(Net, edge.arrow.size =0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
#'}
#'
                                         
simulate_sbm_network = function(N_id = 100,                           # N people
                                    B = NULL,                         # Block tie probabilities
                                    V = 3,                            # Blocking variables
                                    groups=NULL,                      # Group IDs
                                    error_sigma = 0.01,               # Error variance in Gaussian model
                                    outcome_mode="bernoulli",         # outcome mode
                                    link_mode = "logit",              # link mode
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

# Varying effects on individuals
sr = matrix(NA, nrow=N_id, ncol=2)
for( i in 1:N_id){
 sr[i,] = c(0,0) 

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
 dr_scrap = c(0,0)

 if(!is.null(dyadic_predictors)){
  dr_scrap[1] = dr_scrap[1] + sum(dyadic_effects*dyadic_predictors[i,j,]) 
  dr_scrap[2] = dr_scrap[2] + sum(dyadic_effects*dyadic_predictors[j,i,]) 
  }

   B_i_j = B_j_i = c()
  for(v in 1:V){
    B_i_j[v] =  B[[v]][groups[i,v] , groups[j,v] ]
    B_j_i[v] =  B[[v]][groups[j,v] , groups[i,v] ]
  }

 dr[i,j] = dr_scrap[1] + sum(B_i_j)
 dr[j,i] = dr_scrap[2] + sum(B_j_i)

# Simulate outcomes
if(outcome_mode=="bernoulli"){
 if(link_mode=="logit"){
 p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j])
 y_true[i,j] = rbern( 1 , p[i,j] )

 p[j,i] = inv_logit( sr[j,1] + sr[i,2] + dr[j,i])
 y_true[j,i] = rbern( 1 , p[j,i] )
 }

 if(link_mode=="probit"){
 p[i,j] = pnorm( sr[i,1] + sr[j,2] + dr[i,j])
 y_true[i,j] = rbern( 1 , p[i,j] )

 p[j,i] = pnorm( sr[j,1] + sr[i,2] + dr[j,i])
 y_true[j,i] = rbern( 1 , p[j,i] )
 }
 }

if(outcome_mode=="binomial"){
 if(link_mode=="logit"){
 p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j])
 y_true[i,j] = rbinom( 1 , size=samps[i,j], prob=p[i,j] )

 p[j,i] = inv_logit( sr[j,1] + sr[i,2] + dr[j,i])
 y_true[j,i] = rbinom( 1 , size=samps[j,i], prob=p[j,i] )
 }

 if(link_mode=="probit"){
 p[i,j] = pnorm( sr[i,1] + sr[j,2] + dr[i,j])
 y_true[i,j] = rbinom( 1 , size=samps[i,j], prob=p[i,j] )

 p[j,i] = pnorm( sr[j,1] + sr[i,2] + dr[j,i])
 y_true[j,i] = rbinom( 1 , size=samps[j,i], prob=p[j,i] )
 }
 }

if(outcome_mode=="poisson"){
 if(link_mode=="log"){
 p[i,j] = exp( sr[i,1] + sr[j,2] + dr[i,j])
 y_true[i,j] = rpois( 1 , p[i,j] )

 p[j,i] = exp( sr[j,1] + sr[i,2] + dr[j,i])
 y_true[j,i] = rpois( 1 , p[j,i] )
 }
 }

 if(outcome_mode=="gaussian"){
 if(link_mode=="identity"){
 p[i,j] = sr[i,1] + sr[j,2] + dr[i,j]
 y_true[i,j] = rnorm( 1 , p[i,j], error_sigma )

 p[j,i] = sr[j,1] + sr[i,2] + dr[j,i]
 y_true[j,i] = rnorm( 1 , p[j,i], error_sigma )
 }
 }
        }
    }

for ( i in 1:N_id ){
    y_true[i,i] = 0
    p[i,i] = 0
    dr[i,i] = 0
 }

return(list(network=y_true, tie_strength=p, group_ids=groups, individual_predictors=individual_predictors, dyadic_predictors=dyadic_predictors, sr=sr, dr=dr, samps=samps))
}
