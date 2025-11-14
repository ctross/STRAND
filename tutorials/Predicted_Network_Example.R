##############################################################################################
#
#   Understanding the predicted network
#
##############################################################################################
 
# When users fit a block_plus_social_relations_model, or similar models, they can export a 
# predicted network of ties weights. These tie weight are returned, even for dyads that are
# masked out of the regression. The predictions use the supplied covariate data. Sometimes,
# users will want the full set of predictions, other times user will want to post-process the
# predicted network before using it. Lets look at two examples.

# Load libraries
library(STRAND)
library(ggplot2)
library(igraph)

#############################################################################################################
####### Case 1: Here we have a case when different groups of animals are kept in different cages
# Between-cage ties are impossible, and can be treated as strucutural zeros and masked out of the regression

# Load package data
data(Callithrix_Data)

# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[1]]$Aggressed)
mask = list(Aggressed = Callithrix_Data[[1]]$NoOpportunity)               # 0s are possible ties, 1s are censored/impossible. The mask will remove such terms from the likelihood of the model.
dyad = list(RankDiff = standardize_strand(Callithrix_Data[[1]]$RankDiff)) # rank distances can be computed, even between individuals in different cages, but arent necessarily meaningful
indiv =  Callithrix_Data[[1]]$Individual

model_dat_cm = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             mask=mask,                    # if mask != NULL, then it will be used by any model that reads this data.
                             outcome_mode = "bernoulli",
                             link_mode = "logit"
                             )
# Model
fit_cm_1 =  fit_block_plus_social_relations_model(data=model_dat_cm,
                              block_regression = ~ 1,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff, 
                              return_predicted_network =TRUE,
                              mode="mcmc",
                              mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = 11, adapt_delta = 0.95)
)

# Get results. Data behind mask (e.g., in RankDiff) *dont* affect *parameter estimates*. 
res_cm_mask = summarize_strand_results(fit_cm_1)

# Note that the data behind the mask *do* influence the *predicted ties*! 
# Effects of covariates behind the mask are extrapolated to predict for the masked set.
pred_network_1a = apply(res_cm_mask$samples$predicted_network_sample, c(2,3), mean)
pred_network_1b = pred_network_1a*(1-mask[[1]])

# Plot
par(mfrow=c(2,2))
image(pred_network_1a)  # The off-diagonal blocks are meaningless in this case
image(pred_network_1b)  # Remove those predictions using the mask

# Now we can plot binarized network graphs if we want
binarized_network_1a = graph_from_adjacency_matrix(ifelse(pred_network_1a>0.5,1,0), mode="directed") 
binarized_network_1b = graph_from_adjacency_matrix(ifelse(pred_network_1b>0.5,1,0), mode="directed")  

V(binarized_network_1a)$color = c(rep("darkred",5), rep("slateblue",6), rep("grey50",5), rep("orange4",4))
V(binarized_network_1b)$color = c(rep("darkred",5), rep("slateblue",6), rep("grey50",5), rep("orange4",4))

plot(binarized_network_1a, edge.arrow.size = 0.5)    # Not realistic
plot(binarized_network_1b, edge.arrow.size = 0.5)    # Realistic


#######################################################################################################################
####### Case 2: Lets keep the same cage setup as before, but lets assume further that the field-notes on individual-13s
# outgoing ties were lost at random, but we know the relevant covariate data (e.g., sex and rank) for individual 13. 

# Load package data
data(Callithrix_Data)

new_mask = base_mask = Callithrix_Data[[1]]$NoOpportunity
base_outcomes = Callithrix_Data[[1]]$Aggressed

base_outcomes[13,] = 0 # These ties might have happened, but we lost our notes on them. 
new_mask[13,] = 1      # We expand our mask to reflect that we dont know the true outcomes.

# Create the STRAND data object
outcome = list(Aggressed = base_outcomes)
mask = list(Aggressed = new_mask)                                         # 0s are possible ties, 1s are censored/impossible. The mask will remove such terms from the likelihood of the model.
dyad = list(RankDiff = standardize_strand(Callithrix_Data[[1]]$RankDiff)) # rank distances can be computed, even between individuals in different cages, but arent necessarily meaningful
indiv =  Callithrix_Data[[1]]$Individual

model_dat_cm_2 = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             mask=mask,                    # if mask != NULL, then it will be used by any model that reads this data.
                             outcome_mode = "bernoulli",
                             link_mode = "logit"
                             )
# Model
fit_cm_2 =  fit_block_plus_social_relations_model(data=model_dat_cm_2,
                              block_regression = ~ 1,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff,
                              return_predicted_network =TRUE,
                              mode="mcmc",
                              mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = 11, adapt_delta = 0.95)
)

# Get results. Data behind mask dont affect parameter estimates.
res_cm_mask_2 = summarize_strand_results(fit_cm_2)

# Again, lets look at the predicted networks
pred_network_1a = apply(res_cm_mask$samples$predicted_network_sample, c(2,3), mean) # Raw
pred_network_1b = pred_network_1a*(1-mask[[1]])                                     # Mask out everything (including individual 13)  
pred_network_1c = pred_network_1a*(1-base_mask)                                     # Mask out only cross-cage dyads 

par(mfrow=c(2,3))
image(pred_network_1a)  # The off-diagonal blocks are meaningless in this case
image(pred_network_1b)  # This also sets ties involving individual-13 to be zero, not what we want
image(pred_network_1c)  # Cross-cage dyads are masked, but within individual 13s cage, we predict for individual 13

# Now we can plot binarized network graphs if we want
binarized_network_1a = graph_from_adjacency_matrix(ifelse(pred_network_1a>0.5,1,0), mode="directed") 
binarized_network_1b = graph_from_adjacency_matrix(ifelse(pred_network_1b>0.5,1,0), mode="directed")  
binarized_network_1c = graph_from_adjacency_matrix(ifelse(pred_network_1c>0.5,1,0), mode="directed") 

V(binarized_network_1a)$color = c(rep("darkred",5), rep("slateblue",6), rep("grey50",5), rep("orange4",4))
V(binarized_network_1b)$color = c(rep("darkred",5), rep("slateblue",6), rep("grey50",5), rep("orange4",4))
V(binarized_network_1c)$color = c(rep("darkred",5), rep("slateblue",6), rep("grey50",5), rep("orange4",4))

plot(binarized_network_1a, edge.arrow.size = 0.5)    # Not realistic
plot(binarized_network_1b, edge.arrow.size = 0.5)    # Not realistic
plot(binarized_network_1c, edge.arrow.size = 0.5)    # More realistic imputation of network on the basis of rank diff and sex
