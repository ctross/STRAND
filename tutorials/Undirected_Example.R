#########################################################
#
#   Bernoulli Analyses - Simulated undirected data 
#
#########################################################

# Clear working space
rm(list = ls())
set.seed(1)
# Load libraries
library(STRAND)
library(rethinking)
library(ggplot2)

# STRAND is designed to model directed networks. "Undirected" networks can be thought of
# as just a special case of directed networks where reciprocity parameters are unity.

# Make data with this sample size
N_id = 95

# Covariates
Kinship = rlkjcorr( 1 , N_id , eta=1.5 )                 # Dyadic covariates should be symmetric for undirected networks
Dominant = ceiling(rlkjcorr( 1 , N_id , eta=1.5 ) - 0.1) #

Mass = rbern(N_id, 0.4)

# Organize data into arrays
dyadic_preds = array(NA,c(N_id,N_id,2))

dyadic_preds[,,1] = Kinship
dyadic_preds[,,2] = Dominant

# Set effect sizes
sr_mu = c(0,0)  
sr_sigma = c(1.8, 1.8)           # Sender and receiver sds should be the same in "Undirected" networks
sr_rho = 0.999                   # Rho should be 1 in limit in "Undirected" networks
dr_mu = 0 
dr_sigma = 3.5
dr_rho = 0.999                   # Rho should be 1 in limit in "Undirected" networks
sr_effects_1 = c(1.9, 1.9)       # Sender and receiver effects should be the same in "Undirected" networks
dr_effects_1 = c(-1.9, 2.1)       

# Block structure
group_probs_block_size = c(0.25, c(0.5, 0.5)*(1-0.25))

B_1 = matrix(-9,nrow=1,ncol=1)
B_2 = matrix(rnorm(9,0,3),nrow=3,ncol=3)

B_2 = 0.5*B_2 + 0.5*t(B_2)       # Block effects should be symmetric in "Undirected" networks

diag(B_2) = diag(B_2) + 3

B = list(B_1, B_2)
 
groups_1 = rep("Any", N_id) 
groups_2 = sample(c("Red","White","Blue"), size=N_id, replace=TRUE, prob=group_probs_block_size)


groups = data.frame(Intercept=as.numeric(factor(groups_1)), Merica=as.numeric(factor(groups_2)))
groups_f = data.frame(Intercept=factor(groups_1), Merica=factor(groups_2))
individual = data.frame(Mass=Mass)

#################################################### Simulate SBM + SRM network
G = simulate_sbm_plus_srm_network(N_id = N_id, 
                         B = B, 
                         V = 2,
                         groups = groups,                  
                         sr_mu = sr_mu,  
                         sr_sigma = sr_sigma, 
                         sr_rho = sr_rho,
                         dr_mu = dr_mu,  
                         dr_sigma = dr_sigma, 
                         dr_rho = dr_rho,
                         outcome_mode = "binomial", 
                         link_mode = "logit",                 
                         individual_predictors = data.frame(Mass = Mass),
                         dyadic_predictors = dyadic_preds,
                         individual_effects = matrix(sr_effects_1, nrow=2, ncol=1),
                         dyadic_effects = dr_effects_1
                         )  

########################################################################## 
image(G$tie_strength)     # Tie strength is now essentially undirected up to a small amount of rng error      
image(G$network)          # And so is the outcome network, up to some randomness caused by rbinom for weaker ties

# As a hack to make fully "Undirected", we can write:
G$network2 = G$network + t(G$network)  # Fully symmeterize
G$samps2 = G$samps + t(G$samps)        # Fully symmeterize
image(G$network2)  
image(G$samps2)  
image(G$network2/G$samps2)  

################################################### Organize for model fitting
# Add rownames and colnames as needed
name_vec = paste("Individual", 1:N_id)
rownames(G$network2) = colnames(G$network2) = name_vec
rownames(G$samps2) = colnames(G$samps2) = name_vec
rownames(Kinship) = colnames(Kinship) = name_vec
rownames(Dominant) = colnames(Dominant) = name_vec
rownames(groups_f) = name_vec
rownames(individual) = name_vec

model_dat = make_strand_data(outcome=list(Outcome = G$network2),  
                             block_covariates=groups_f, 
                             individual_covariates=individual, 
                             dyadic_covariates=list(Kinship=Kinship, Dominant=Dominant),
                             exposure = list(Outcome = G$samps2),  
                             outcome_mode = "binomial", 
                             link_mode="logit")

# Model the data with STRAND
fit =  fit_block_plus_social_relations_model(data=model_dat,
                            block_regression = ~ Merica,
                              focal_regression = ~ Mass,             # Make sure to include the same individual predictors for both focal and target effects
                              target_regression = ~ Mass,            #
                              dyad_regression = ~ Kinship + Dominant,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = NULL, adapt_delta = 0.95)
)

# Check parameter recovery
res = summarize_strand_results(fit)

# Check that the estimated effects are close to the generative ones

############################### Plots
vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, site="HP", only_slopes=TRUE)
vis_1

################################################################# Reciprocity terms
# Simple correlations and basic 4-way variance partition
 strand_VPCs(fit, n_partitions = 4, include_reciprocity = TRUE, mode="cor")

# In data, the VPCs are about 
 c(sr_sigma, dr_sigma, sqrt(3.14^2 / 3))^2 / sum(c(sr_sigma, dr_sigma, sqrt(3.14^2 / 3))^2)

# note: sqrt(3.14^2 / 3)) is the error variance implied by the logit link


