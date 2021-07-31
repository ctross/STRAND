#' Simulate a self-reported network 
#'
#' This function allows the user to simulate a self-reported network, along with a 'true' network, and a network of resource flows. 
#' The function first simulates the true network data using the simulate_sbm_plus_srm_network() function, and then simulate self-reports and observable resource flows over the true network.
#' This allows the user to investigate the effects of response biases, such as false positive rate, on network properties.
#'
#' @param 
#' N_id Number of individuals.
#' @param 
#' N_groups Number of groups.
#' @param 
#' group_probs A vector of the probabilities of individuals being in each group.
#' @param 
#' B Tie probabilities between group blocks.  If B is null, then in_block and out_block can be specified to create B.
#' @param 
#' in_block Tie probabilities between members of the same group. This is overridden by a non-NULL B parameter.
#' @param 
#' out_block Tie probabilities between members of different groups. This is overridden by a non-NULL B parameter. 
#' @param 
#' sr_mu Mean vector for sender and receivier random effects. In most cases, this should be c(0,0).
#' @param 
#' dr_mu Mean vector for dyadic random effects. In most cases, this should be c(0,0).
#' @param 
#' sr_sigma Standard deviation vector for sender and receivier random effects. The first element controls node-level variation in out-degree, the second in in-degree.
#' @param 
#' dr_sigma Standard deviation for dyadic random effects.
#' @param 
#' sr_rho Correlation of sender-receiver effects (i.e., generalized reciprocity).
#' @param 
#' dr_rho Correlation of dyad effects: (i.e., dyadic reciprocity).
#' @param 
#' individual_predictors An N_id by N_individual_parameters matrix of covariates.
#' @param 
#' dyadic_predictors An N_id by N_id by N_dyadic_parameters matrix of covariates.
#' @param 
#' individual_effects A 2 by N_individual_parameters matrix of slopes. The first row gives effects of focal characteristics (on out-degree). 
#' The second row gives effects of target characteristics (on in-degree).
#' @param 
#' dyadic_effects An N_dyadic_parameters vector of slopes.
#' @param 
#' fpr_effects A 3 by N_predictors matrix of slopes. The first row controls the effects of covariates on layer 1 responses, and the second row layer 2 responses. 
#' The third row controls the effects of covariates on false positive in the observation network.  
#' @param 
#' rtt_effects A 3 by N_predictors matrix of slopes. The first row controls the effects of covariates on layer 1 responses, and the second row layer 2 responses. 
#' The third row controls the effects of covariates on false positive in the observation network.  
#' @param 
#' theta_effects An N_predictors vector of slopes controling the effects of covariates on name duplication from layer 1 to layer 2.
#' @param 
#' false_positive_rate The baseline false positive rate. This should be supplied as a 3 vector, with each element controlling the false positive rate for a single layer. 
#' Support is on the unit interval. If covariates are centered, this can be loosly thought of as an average false positive rate.
#' @param 
#' recall_of_true_ties The baseline recall rate of true ties. This should be supplied as a 3 vector, with each element controlling the true tie recall rate for a single layer. 
#' Support is on the unit interval. If covariates are centered, this can be loosly thought of as an average true tie recall rate.
#' @param 
#' theta_mean The baseline probability of name duplication from layer 1 to layer 2. Scalar.
#' Support is on the unit interval. If covariates are centered, this can be loosly thought of as an average true tie recall rate.
#' @param 
#' fpr_sigma Standard deviation 3-vector for false_positive_rate random effects. There should be one value for each layer.
#' @param 
#' rtt_sigma Standard deviation 3-vector for recall_of_true_ties random effects. There should be one value for each layer.
#' @param 
#' theta_sigma Standard deviation scalar for theta random effects. 
#' @param 
#' N_responses Number of self-report layers. =1 for single-sampled, =2 for double-sampled.
#' @param 
#' N_periods Number of time-periods in which observed transfers are sampled. 
#' @param 
#' flow_rate A vector of length N_periods, each element controls the probability that a true tie will result in a transfer in a given time period.
#' @param 
#' decay_curve A vector of length N_periods, each element controls the increment log-odds of recalling a true tie at time T, as a function of a transfer occuring in period t.
#' @return A list of data formatted for use in Stan models.
#' @export
#' @examples
#' \dontrun{
#' B = diag(3)*0.01
#' B[1,3] = 0.001
#' B[3,2] = 0.002
#' B[1,1] = 0.02
#' 
#' A = simulate_selfreport_network(N_id = 100, B=B, individual_predictor=matrix(rnorm(100, 0, 1), nrow=100,ncol=1), individual_effects=matrix(c(1.2, 0.5),ncol=1,nrow=2))
#' Net = graph_from_adjacency_matrix(A$reporting_network[,,1], mode = c("directed"))
#' V(Net)$color = c("turquoise4","gray13", "goldenrod3")[A$group_ids]
#' 
#' plot(Net, edge.arrow.size =0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
#' }
#'

###########################################################################################
simulate_selfreport_network = function(  N_id = 99,                        # Number of respondents
                                         N_groups = 3,                     # Number of block, aka social groups
                                         group_probs = c(0.2, 0.5, 0.3),   # Density of each group in overall network
                                         B = NULL,                         # Block tie probabilities
                                         in_block = 0.02,                  # Tie probability wthin groups  
                                         out_block = 0.01,                 # Tie probability between groups
                                         sr_mu = c(0,0),                   # Average sender (cell 1) and reciever (cell 2) effect log odds
                                         sr_sigma = c(1,1),                # Sender (cell 1) and reciever (cell 2) effect variances 
                                         sr_rho = 0.6,                     # Correlation of sender and reciever effects
                                         dr_mu = c(0,0),                   # Average i to j dyad effect (cell 1) and j to i dyad effect (cell 2) log odds
                                         dr_sigma = 1,                     # Variance of dyad effects 
                                         dr_rho = 0.7,                     # Correlation of i to j dyad effect and j to i dyad effect 
                                         individual_predictors = NULL,     # A matrix of covariates
                                         dyadic_predictors = NULL,         # An array of covariates
                                         individual_effects = NULL,        # The effects of predictors on sender effects (row 1) and receiver effects (row 2)
                                         dyadic_effects = NULL,            # The effects of predictors on dyadic ties
                                         fpr_effects = NULL,               # False positive rate of effect 1: ego->alter, alter->ego, goods from ego -> alter
                                         rtt_effects = NULL,               # Recall of true ties: ego->alter, alter->ego, goods from ego -> alter
                                         theta_effects = NULL,            
                                         false_positive_rate = c(0.01, 0.02, 0.00), 
                                         recall_of_true_ties = c(0.8, 0.6, 0.99),
                                         theta_mean = 0.125, 
                                         fpr_sigma = c(0.3, 0.2, 0.0), 
                                         rtt_sigma = c(0.5, 0.2, 0.0),
                                         theta_sigma = 0.2,
                                         N_responses = 2,
                                         N_periods = 12, 
                                         flow_rate = rbeta(12, 3, 30),
                                         decay_curve = rev(1.5*exp(-seq(1,4, length.out=12)))
                                         ){

# Reports of food transferred from i to j in layer 1 & from j to i in layer 2
statement_flows = array( NA , dim=c( N_id , N_id ,  N_responses ) )
# Food transferred from i to j in every period
goods_flows = array( NA , dim=c( N_id , N_id , N_periods ) )

# y_true is a true network of directed giving ties from ego to alter
G_net = simulate_sbm_plus_srm_network(N_id = N_id, N_groups = N_groups, group_probs = group_probs, B=B,
                                      in_block = in_block, out_block = out_block,
                                      sr_mu = sr_mu,  sr_sigma = sr_sigma, sr_rho = sr_rho,
                                      dr_mu = dr_mu,  dr_sigma = dr_sigma, dr_rho = dr_rho,
                                      individual_predictors = individual_predictors,   
                                      dyadic_predictors = dyadic_predictors,        
                                      individual_effects = individual_effects,        
                                      dyadic_effects = dyadic_effects
                                      )

group_id = G_net$group_ids
y_true = G_net$network


# Link functions for individual-level parameters
theta = rtt_goods = fpr_goods = rtt_rev = fpr_rev = rtt = fpr = rep(NA, N_id)
for(i in 1:N_id){

 if(!is.null(individual_predictors)){
  ef_fpr = sum(fpr_effects[1,]*individual_predictors[i,])
  ef_rtt = sum(rtt_effects[1,]*individual_predictors[i,])
  ef_fpr_rev = sum(fpr_effects[2,]*individual_predictors[i,])
  ef_rtt_rev = sum(rtt_effects[2,]*individual_predictors[i,])
  ef_theta = sum(theta_effects*individual_predictors[i,])
  ef_fpr_goods = sum(fpr_effects[3,]*individual_predictors[i,])
  ef_rtt_goods = sum(rtt_effects[3,]*individual_predictors[i,])
  }
  else{
  ef_fpr = 0
  ef_rtt = 0
  ef_fpr_rev = 0
  ef_rtt_rev = 0
  ef_theta = 0
  ef_fpr_goods = 0
  ef_rtt_goods = 0
  }

fpr[i] = rnorm(1, logit(false_positive_rate[1]) + ef_fpr , fpr_sigma[1])
rtt[i] = rnorm(1, logit(recall_of_true_ties[1]) + ef_rtt , rtt_sigma[1])

fpr_rev[i] = rnorm(1, logit(false_positive_rate[2]) + ef_fpr_rev , fpr_sigma[2])
rtt_rev[i] = rnorm(1, logit(recall_of_true_ties[2]) + ef_rtt_rev , rtt_sigma[2])

theta[i] = inv_logit(rnorm(1, logit(theta_mean) + ef_theta , theta_sigma[1]))

fpr_goods[i] = rnorm(1, logit(false_positive_rate[3]) + ef_fpr_goods , fpr_sigma[3])
rtt_goods[i] = rnorm(1, logit(recall_of_true_ties[3]) + ef_rtt_goods , rtt_sigma[3])
}

fpr_out = cbind(fpr, fpr_rev, fpr_goods)
rtt_out = cbind(rtt, rtt_rev, rtt_goods)

# Simulate flow of goods per time period
# independent random sample from each time period. Need to include a covariate about having food to share? Censoring? -1 no food to share, 1 food to share and shared, 0 food to share and not shared
for (g in 1:N_periods){
for ( i in 1:N_id ){
  for ( j in 1:N_id ){
    if(i != j){
      if(y_true[i,j]==0){
        goods_flows[i, j, g] = rbern(1, inv_logit(fpr_goods[i]))
      } else{
       goods_flows[i, j, g] = rbern(1,  inv_logit(rtt_goods[i]))*rbern(1,flow_rate[g])
      }
      
    }
    }}}

# Simulate statements about flows
for ( i in 1:N_id ){
  for ( j in 1:N_id ){
    if(i != j){
      ###### ego reporting on who ego gave food to 
      weighted_offset_ij = sum(decay_curve[which(goods_flows[i, j, ]==1)])
  
      statement_flows[i, j, 1] = rbern(1,  
                                 inv_logit(fpr[i])*(1-y_true[i,j]) + 
                                 inv_logit(rtt[i] + weighted_offset_ij)*y_true[i,j])
      }
      }}

# (below) code these ties in reverse order since question is reversed. if accuracy is perfect, then statement_flows[j, i, 1] = statement_flows[j, i, 2]? make sure to check
 if(N_responses == 2){

  for ( i in 1:N_id ){
   for ( j in 1:N_id ){
    if(i != j){
        ###### ego reporting on who gave food to ego
        weighted_offset_ji = sum(decay_curve[which(goods_flows[j, i, ]==1)])

        statement_flows[i, j, 2] = rbern(1, 
                                               inv_logit(fpr_rev[i])*(1-y_true[j,i]) + 
                                 inv_logit(rtt_rev[i] + weighted_offset_ji)*y_true[j,i])


     # Question contaminiation bias. That is a bias to over-represent reciprical relations, even when not true
      if(statement_flows[i, j, 1]==1 & statement_flows[i, j, 2]==0)
       statement_flows[i, j, 2] = rbern(1, theta[i])
         }
         }}

  
}

for(i in 1:N_id){
 for (g in 1:N_periods)
 goods_flows[i,i,g] <- 0

 statement_flows[i,i,1] <- 0
 if(N_responses == 2)
 statement_flows[i,i,2] <- 0
}



return(list(true_network=y_true, transfer_network=goods_flows, reporting_network=statement_flows, group_ids=group_id, N_id = N_id,
            N_groups = N_groups, N_periods=N_periods, N_responses=N_responses, sr = G_net$sr, dr=G_net$dr, 
            fpr=fpr_out, rtt=rtt_out, theta=theta ))
}


