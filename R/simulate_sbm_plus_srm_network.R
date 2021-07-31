#' A function to simulate single layer directed networks 
#'
#' This is a function to simulate single layer network data with block structure, sender-receiver effects, and dyadic reciprocity. This function
#' is essentially the union of a social relations model and a stochastic block model.
#'
#' @param 
#' N_id Number of individuals.
#' @param 
#' N_groups Number of groups.
#' @param 
#' group_probs A vector of the probabilities of being in each group.
#' @param 
#' B Tie probabilities between group blocks.  If B is null, then in_block and out_block can be specified to create B.
#' @param 
#' in_block Tie probabilities between memebers of the same group. Overridden by a non-NULL B parameter.
#' @param 
#' out_block Tie probabilities between memebers of different groups.  Overridden by a non-NULL B parameter. 
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
#' individual_predictors An N_id by N_individual_parameters matrix of covariates.
#' @param 
#' dyadic_predictors An N_id by N_id by N_dyadic_parameters array of covariates.
#' @param 
#' individual_effects A 2 by N_individual_parameters matrix of slopes. The first row gives effects of focal characteristics (on out-degree). 
#' The second row gives effects of target characteristics (on in-degree).
#' @param 
#' dyadic_effects An N_dyadic_parameters vector of slopes.
#' @return A list of 5 objects: network (an adjacency matrix of binary outcomes), tie_strength (an adjacency matrix with probability weights), 
#' group_ids (a vector of length N_id, giving the group of each individual), individual_predictors (the supplied covariate data is saved along with the network data), 
#' and dyadic_predictors (the supplied covariate data is saved along with the network data).
#' @export
#' @examples
#' \dontrun{
#' B = diag(3)*0.001
#' B[1,3] = 0.0001
#' B[3,2] = 0.0007
#' B[1,1] = 0.002
#' 
#' A = simulate_sbm_plus_srm_network(N_id = 100, B=B, individual_predictor=matrix(rnorm(100, 0, 1), nrow=100,ncol=1), individual_effects=matrix(c(1.2, 0.5),ncol=1,nrow=2))
#' Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
#' V(Net)$color = c("turquoise4","gray13", "goldenrod3")[A$group_ids]
#' 
#' plot(Net, edge.arrow.size =0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
#' }
#'
                                         
simulate_sbm_plus_srm_network = function(N_id = 90,                        # Number of respondents
                                         N_groups = 3,                     # Number of block, aka social groups
                                         group_probs = c(0.2, 0.5, 0.3),   # Density of each group in overall network
                                         B = NULL,                         # Block tie probabilities
                                         in_block = 0.02,                  # Tie probability wthin groups  
                                         out_block = 0.01,                 # Tie probability between groups
                                         sr_mu = c(0,0),                   # Average sender (cell 1) and reciever (cell 2) effect log odds
                                         sr_sigma = c(0.3, 2),             # Sender (cell 1) and reciever (cell 2) effect variances 
                                         sr_rho = 0.6,                     # Correlation of sender and reciever effects
                                         dr_mu = c(0,0),                   # Average i to j dyad effect (cell 1) and j to i dyad effect (cell 2) log odds
                                         dr_sigma = 3,                     # Variance of dyad effects 
                                         dr_rho = 0.7,                     # Correlation of i to j dyad effect and j to i dyad effect 
                                         individual_predictors = NULL,     # A matrix of covariates
                                         dyadic_predictors = NULL,         # An array of covariates
                                         individual_effects = NULL,        # The effects of predictors on sender effects (row 1) and receiver effects (row 2)
                                         dyadic_effects = NULL             # The effects of predictors on dyadic ties
                                         )
{
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

  if(!is.null(B)){
   if(dim(B)[1] != N_groups | dim(B)[2] != N_groups | !is.matrix(B)){
    stop("B must be a matrix with dimensions N_groups by N_groups ")
   }}

# Sample respondents into groups
groups = sample( 1:N_groups , size=N_id , replace=TRUE , prob=group_probs )

# Define interaction matrix across groups
if(is.null(B)){
 B = diag(N_groups)
 for ( i in 1:length(B) ) if ( B[i]==0 ) B[i] = out_block
 for ( i in 1:length(B) ) if ( B[i]==1 ) B[i] = in_block
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
p = y_true = matrix(NA, N_id, N_id)
# Loop over upper triangle and create ties from i to j, and j to i
for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
# Dyadic effects
 dr = rmvnorm2(1, Mu=c(dr_mu), sigma=rep(dr_sigma,2), Rho=Rho_dr)

 if(!is.null(dyadic_predictors)){
  dr[1] = dr[1] + sum(dyadic_effects*dyadic_predictors[i,j,]) 
  dr[2] = dr[2] + sum(dyadic_effects*dyadic_predictors[j,i,]) 
  }

# Simulate outcomes
 p[i,j] = inv_logit( logit(B[ groups[i] , groups[j] ]) + sr[i,1] + sr[j,2] + dr[1])
 y_true[i,j] = rbern( 1 , p[i,j] )

 p[j,i] = inv_logit( logit(B[ groups[j] , groups[i] ]) + sr[j,1] + sr[i,2] + dr[2])
 y_true[j,i] = rbern( 1 , p[j,i] )
        }
    }

for ( i in 1:N_id ){
    y_true[i,i] = 0
    p[i,i] = 0
 }

return(list(network=y_true, tie_strength=p, group_ids=groups, individual_predictors=individual_predictors, dyadic_predictors=dyadic_predictors))
}



