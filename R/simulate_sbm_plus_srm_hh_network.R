#' A function to simulate single layer directed networks using a combined stochastic block and social relations model with household effects
#'
#' This is a function to simulate single layer network data with a stochastic block structure, sender-receiver effects, and dyadic reciprocity. This function
#' is essentially the union of a social relations model and a stochastic block model. Includes SRM structure at both individual and household level.
#'
#' @param 
#' N_id Number of individuals.
#' @param 
#' N_hh Number of households.
#' @param 
#' B List of matrices that hold intercept and offset terms. Log-odds. The first matrix should be  1 x 1 with the value being the intercept term.
#' @param 
#' V Number of blocking variables in B.
#' @param 
#' groups Dataframe of the block IDs of each individual for each variable in B.
#' @param 
#' sr_mu Mean A vector for sender and receivier random effects. In most cases, this should be c(0,0).
#' @param 
#' dr_mu Mean A vector for dyadic random effects. In most cases, this should be c(0,0).
#' @param 
#' sr_sigma A standard deviation vector for sender and receivier random effects. The first element controls node-level variation in out-degree, the second in in-degree.
#' @param 
#' dr_sigma Standard deviation for dyadic random effects.
#' @param 
#' sr_rho Correlation of sender-receiver effects (i.e., generalized reciprocity).
#' @param 
#' dr_rho Correlation of dyad effects (i.e., dyadic reciprocity).
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
#' B = matrix(-8, nrow=G, ncol=G)
#' diag(B) = -4.5 # Block matrix
#'
#' B[1,3] = -5.9
#' B[3,2] = -6.9
#' 
#' A = simulate_sbm_plus_srm_hh_network(N_id = N_id, B=list(B=B), V=V, 
#'                          groups=data.frame(clique=factor(clique)),
#'                          individual_predictor=matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1), 
#'                          individual_effects=matrix(c(1.7, 0.3),ncol=1, nrow=2),
#'                          sr_sigma = c(1.4, 0.8), sr_rho = 0.5,
#'                          dr_sigma = 1.2, dr_rho = 0.8,
#'                          mode="bernoulli"
#'                                )
#'
#' Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
#' V(Net)$color = c("turquoise4","gray13", "goldenrod3")[A$group_ids$clique]
#' 
#' plot(Net, edge.arrow.size =0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
#'}
#'
                                                                                               
simulate_sbm_plus_srm_hh_network = function( N_id = 99,                        # Number of respondents
                                             N_hh = 30,                        # Number of households
                                             HH = NULL,                        # Household ID of each respondent
                                             B = NULL,                         # Tie probabilities
                                             V = 3,                            # Blocking variables
                                             groups=NULL,                      # Group IDs
                                             sr_mu = c(0,0),                   # Average sender (cell 1) and reciever (cell 2) effect log odds
                                             sr_sigma = c(0.3, 1.5),           # Sender (cell 1) and reciever (cell 2) effect variances 
                                             sr_rho = 0.6,                     # Correlation of sender and reciever effects
                                             dr_mu = c(0,0),                   # Average i to j dyad effect (cell 1) and j to i dyad effect (cell 2) log odds
                                             dr_sigma = 1,                     # Variance of dyad effects 
                                             dr_rho = 0.7,                     # Correlation of i to j dyad effect and j to i dyad effect 
                                             hh_sr_mu = c(0,0),                # Average sender (cell 1) and reciever (cell 2) effect log odds
                                             hh_sr_sigma = c(0.3, 1.5),        # Sender (cell 1) and reciever (cell 2) effect variances 
                                             hh_sr_rho = 0.6,                  # Correlation of sender and reciever effects
                                             hh_dr_mu = c(0,0),                # Average i to j dyad effect (cell 1) and j to i dyad effect (cell 2) log odds
                                             hh_dr_sigma = 1,                  # Variance of dyad effects 
                                             hh_dr_rho = 0.7,                  # Correlation of i to j dyad effect and j to i dyad effect
                                             hh_diag_offset = 2,               # Log odds offset for within-household ties
                                             mode="bernoulli",                 # outcome mode
                                             individual_predictors = NULL,     # A matrix of covariates
                                             dyadic_predictors = NULL,         # An array of covariates
                                             individual_effects = NULL,        # The effects of predictors on sender effects (row 1) and receiver effects (row 2)
                                             dyadic_effects = NULL,            # The effects of predictors on dyadic ties
                                             hh_individual_predictors = NULL,  # A matrix of covariates
                                             hh_dyadic_predictors = NULL,      # An array of covariates
                                             hh_individual_effects = NULL,     # The effects of predictors on sender effects (row 1) and receiver effects (row 2)
                                             hh_dyadic_effects = NULL          # The effects of predictors on dyadic ties
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

   ######### Household parameters
   if(!is.null(hh_individual_predictors)){
   if(is.null(hh_individual_effects)){
    stop("If hh_individual_predictors is supplied, a matching matrix of hh_individual_effects must be supplied.")
   }

   if(length(hh_individual_effects[1,]) != length(hh_individual_predictors[1,])){
    stop("The number of columns of hh_individual_effects must match that of hh_individual_predictors.")
   }

   if(nrow(hh_individual_effects) != 2 ){
    stop("The number of rows of hh_individual_effects must be 2: one for sender effects, the other for receiver. Un-needed slopes can be set to 0.")
   }

   if(nrow(hh_individual_predictors) != N_hh ){
    stop("The number of rows of hh_individual_predictors must be N_hh.")
   }
   }

# Create ind-level correlation matrices (aka matrixes)
Rho_sr = Rho_dr = diag(c(1,1))
Rho_sr[1,2] = Rho_sr[2,1] = sr_rho
Rho_dr[1,2] = Rho_dr[2,1] = dr_rho

# Create hh-level correlation matrices (aka matrixes)
hh_Rho_sr = hh_Rho_dr = diag(c(1,1))
hh_Rho_sr[1,2] = hh_Rho_sr[2,1] = hh_sr_rho
hh_Rho_dr[1,2] = hh_Rho_dr[2,1] = hh_dr_rho

# Varying effects on individuals
sr = matrix(NA, nrow=N_id, ncol=2)
for( i in 1:N_id){
 sr[i,] = rmvnorm2(1 , Mu=sr_mu, sigma=sr_sigma, Rho=Rho_sr ) 

 if(!is.null(individual_predictors)){
  sr[i,1] = sr[i,1] + sum(individual_effects[1,]*individual_predictors[i,]) 
  sr[i,2] = sr[i,2] + sum(individual_effects[2,]*individual_predictors[i,]) 
  }
 }

# Varying effects on households
hh_sr = matrix(NA, nrow=N_hh, ncol=2)
for( i in 1:N_hh){
 hh_sr[i,] = rmvnorm2(1 , Mu=hh_sr_mu, sigma=hh_sr_sigma, Rho=hh_Rho_sr ) 

 if(!is.null(hh_individual_predictors)){
  hh_sr[i,1] = hh_sr[i,1] + sum(hh_individual_effects[1,]*hh_individual_predictors[i,]) 
  hh_sr[i,2] = hh_sr[i,2] + sum(hh_individual_effects[2,]*hh_individual_predictors[i,]) 
  }
 }

 # HH dyadic effects. Loop over upper triangle and create ties from i to j, and j to i
 hh_dr = matrix(NA, nrow=N_hh, ncol=N_hh)
 for( i in 1:(N_hh-1) ){
  for( j in (i+1):N_hh){
    hh_dr_scrap = rmvnorm2(1, Mu=c(hh_dr_mu), sigma=rep(hh_dr_sigma,2), Rho=hh_Rho_dr)

   if(!is.null(hh_dyadic_predictors)){
    hh_dr_scrap[1] = hh_dr_scrap[1] + sum(hh_dyadic_effects*hh_dyadic_predictors[i,j,]) 
    hh_dr_scrap[2] = hh_dr_scrap[2] + sum(hh_dyadic_effects*hh_dyadic_predictors[j,i,]) 
    }

    hh_dr[i,j] = hh_dr_scrap[1] 
    hh_dr[j,i] = hh_dr_scrap[2] 
 }}

 for( i in 1:N_hh){
  hh_dr[i,i] = hh_diag_offset + rnorm(1, hh_dr_mu[1], hh_dr_sigma)
 }

# Build true network
dr = p = y_true = matrix(NA, N_id, N_id)
samps = matrix(rpois(N_id^2, 15), nrow=N_id, ncol=N_id)

# Loop over upper triangle and create ties from i to j, and j to i
for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
# Dyadic effects
 dr_scrap = rmvnorm2(1, Mu=c(dr_mu), sigma=rep(dr_sigma,2), Rho=Rho_dr)

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
if(mode=="bernoulli"){
 p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j] + hh_sr[HH[i],1] + hh_sr[HH[j],2] + hh_dr[HH[i],HH[j]])
 y_true[i,j] = rbern( 1 , p[i,j] )

 p[j,i] = inv_logit( sr[j,1] + sr[i,2] + dr[j,i] + hh_sr[HH[j],1] + hh_sr[HH[i],2] + hh_dr[HH[j],HH[i]])
 y_true[j,i] = rbern( 1 , p[j,i] )
 }

 if(mode=="poisson"){
 p[i,j] = exp( sr[i,1] + sr[j,2] + dr[i,j] + hh_sr[HH[i],1] + hh_sr[HH[j],2] + hh_dr[HH[i],HH[j]])
 y_true[i,j] = rpois( 1 , p[i,j] )

 p[j,i] = exp( sr[j,1] + sr[i,2] + dr[j,i] + hh_sr[HH[j],1] + hh_sr[HH[i],2] + hh_dr[HH[j],HH[i]])
 y_true[j,i] = rpois( 1 , p[j,i] )
 }

 if(mode=="binomial"){
 p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j] + hh_sr[HH[i],1] + hh_sr[HH[j],2] + hh_dr[HH[i],HH[j]])
 y_true[i,j] = rbinom( 1 , size=samps[i,j], prob=p[i,j] )

 p[j,i] = inv_logit( sr[j,1] + sr[i,2] + dr[j,i] + hh_sr[HH[j],1] + hh_sr[HH[i],2] + hh_dr[HH[j],HH[i]])
 y_true[j,i] = rbinom( 1 , size=samps[j,i], prob=p[j,i] )
 }
  }
    }

for ( i in 1:N_id ){
    y_true[i,i] = 0
    p[i,i] = 0
    dr[i,i] = 0
 }

return(list(network=y_true, tie_strength=p,  group_ids=groups, individual_predictors=individual_predictors, dyadic_predictors=dyadic_predictors, hh_ids=HH,
            hh_individual_predictors=hh_individual_predictors, hh_dyadic_predictors=hh_dyadic_predictors, sr=sr, dr=dr, hh_sr=hh_sr, hh_dr=hh_dr, samps=samps))
}



