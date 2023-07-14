#####################################
#
#   Binomial Analyses - Simulated data  
#
########################################

# Clear working space
rm(list = ls())
set.seed(1)
# Load libraries
library(STRAND)
library(rethinking)
library(ggplot2)


# Make data
N_id = 100

# Covariates
Kinship = rlkjcorr( 1 , N_id , eta=1.5 )
Dominant = ceiling(rlkjcorr( 1 , N_id , eta=1.5 ) - 0.1)
Mass = rbern(N_id, 0.4)

# Organize into list
dyadic_preds = array(NA,c(N_id,N_id,3))

dyadic_preds[,,1] = Kinship
dyadic_preds[,,2] = Dominant
dyadic_preds[,,3] = Kinship*Dominant

# Set effect sizes
sr_mu = c(0,0)  
sr_sigma = c(2.2, 1.7) 
sr_rho = 0.55
dr_mu = c(0,0) 
dr_sigma = 1.5
dr_rho= 0.6
sr_effects_1 = c(1.9, 1.3)
dr_effects_1 = c(1.2, 1.7, -2.2)

# Block structure
group_probs_block_size = c(0.25, c(0.25, 0.25)*(1-0.25))

B_1 = matrix(-10,nrow=1,ncol=1)
B_2 = matrix(rnorm(9,0,3),nrow=3,ncol=3)
B_3 = matrix(rnorm(4,0,3),nrow=2,ncol=2)

diag(B_2) = diag(B_2) + 2
diag(B_3) = diag(B_3) + 3.5

B=list(B_1, B_2, B_3)
 
groups_1 = rep("Any",N_id) 
groups_2 = sample( c("Red","White","Blue") , size=N_id , replace=TRUE , prob=group_probs_block_size )
groups_3 = sample( c("Strange", "Charm") , size=N_id , replace=TRUE , prob=c(0.5,0.5) )

groups = data.frame(Intercept=as.numeric(factor(groups_1)), Merica=as.numeric(factor(groups_2)), Quantum=as.numeric(factor(groups_3)))
groups_f = data.frame(Intercept=factor(groups_1), Merica=factor(groups_2), Quantum=factor(groups_3))

#################################################### Simulate SBM + SRM network
G = simulate_sbm_plus_srm_network(N_id = N_id, 
                         B = B, 
                         V=3,
                         groups=groups,                  
                         sr_mu = sr_mu,  
                         sr_sigma = sr_sigma, 
                         sr_rho = sr_rho,
                         dr_mu = dr_mu,  
                         dr_sigma = dr_sigma, 
                         dr_rho = dr_rho,
                         mode="binomial",                  
                         individual_predictors = data.frame(Mass=Mass),
                         dyadic_predictors = dyadic_preds,
                         individual_effects = matrix(sr_effects_1,nrow=2,ncol=1),
                         dyadic_effects = dr_effects_1
                         )        

################################################### Organize for model fitting
model_dat = make_strand_data(outcome=list(G$network),  block_covariates=groups_f, individual_covariates=data.frame(Mass=Mass), 
                           dyadic_covariates=list(Kinship=Kinship, Dominant=Dominant),  outcome_mode = "binomial", exposure=list(G$samps))

# Model the data with STRAND
fit =  fit_block_plus_social_relations_model(data=model_dat,
                            block_regression = ~ Merica + Quantum,
                              focal_regression = ~ Mass,
                              target_regression = ~ Mass,
                              dyad_regression = ~ Kinship*Dominant,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = NULL, adapt_delta = .9)
)

# Check parameter recovery
res = summarize_strand_results(fit)

############################### Plots
vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, site="HP", only_slopes=TRUE)
vis_1


vis_2 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, site="HP", only_technicals=TRUE, only_slopes=FALSE)
vis_2

##### Check all of the block parameters
B_2_Pred = matrix(res$summary_list$`Other estimates`[4:12,2], nrow=3, ncol=3, byrow=TRUE) # Blue, Red, White
plot(B_2_Pred~B_2)

B_3_Pred = matrix(res$summary_list$`Other estimates`[13:16,2], nrow=2, ncol=2, byrow=TRUE) # Charm, Strange
plot(B_3_Pred~B_3)

# NOTE: The block offsets are only identified relative to one-another. Calculate contrasts to look for differences between parameters.

 
