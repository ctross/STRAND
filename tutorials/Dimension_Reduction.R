##############################################
#
#   Binomial Analyses - Simulated data with interaction  
#
########################################

# Clear working space

# Load libraries
library(igraph)
library(rethinking)
library(ggplot2)
library(STRAND)

# Make data
set.seed(1)
N_id = 60

# Covariates
Kinship = STRAND::standardize(rlkjcorr( 1 , N_id , eta=1.5 ))
Dominant = ceiling(rlkjcorr( 1 , N_id , eta=1.5 ) - 0.1)
Mass = rbern(N_id, 0.4)

# Organize into list
dyadic_preds = array(NA,c(N_id,N_id,3))

dyadic_preds[,,1] = Kinship
dyadic_preds[,,2] = Dominant
dyadic_preds[,,3] = Kinship*Dominant

# Set effect sizes
sr_mu = c(0,0)  
sr_sigma = c(2.0, 1.7) 
sr_rho = 0.55
dr_mu = 0 
dr_sigma = 1.8
dr_rho= 0.6
sr_effects_1 = c(1.9, 1.3)
dr_effects_1 = c(0.6, 1.7, -1.2)

# Block structure
group_probs_block_size = c(0.25, c(0.25, 0.25)*(1-0.25))

B_1 = matrix(-11,nrow=1,ncol=1)
B_2 = matrix(rnorm(9,1,3),nrow=3,ncol=3)
B_3 = matrix(rnorm(4,0,3),nrow=2,ncol=2)

diag(B_2) = diag(B_2) + 2
diag(B_3) = diag(B_3) + 1.7

B=list(B_1, B_2, B_3)
 
groups_1 = rep("Any",N_id) 
groups_2 = sample( c("Red","White","Blue") , size=N_id , replace=TRUE , prob=group_probs_block_size )
groups_3 = sample( c("Strange", "Charm") , size=N_id , replace=TRUE , prob=c(0.5,0.5) )

groups = data.frame(Intercept=as.numeric(factor(groups_1)), Merica=as.numeric(factor(groups_2)), Quantum=as.numeric(factor(groups_3)))
groups_f = data.frame(Intercept=factor(groups_1), Merica=factor(groups_2), Quantum=factor(groups_3))
individual = data.frame(Mass=Mass)

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
                         outcome_mode="bernoulli", 
                         link_mode="logit",                 
                         individual_predictors = data.frame(Mass=Mass),
                         dyadic_predictors = dyadic_preds,
                         individual_effects = matrix(sr_effects_1,nrow=2,ncol=1),
                         dyadic_effects = dr_effects_1
                         )        


Net = graph_from_adjacency_matrix(G$network, mode = c("directed"))
V(Net)$color = c("turquoise4","gray13", "goldenrod3")[G$group_ids$Merica]

par(mfrow=c(1,2))
plot(Net, edge.arrow.size = 0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
image(G$network)


############################################# Simulate multiplex networks from the latent network
M = 5 # Network layers
EE = 10 # Samples per dyad
alpha = matrix(NA, nrow=M, ncol=2)

alpha[1,] = c(-3, 6)    # Feeding
alpha[2,] = c(0, 6)     # Territory
alpha[3,] = c(-7, 8)    # Courtship
alpha[4,] = c(2, 0.75)  # Nesting
alpha[5,] = c(4, -6)    # Attacking

for(m in 1:M){
print(paste(inv_logit(alpha[m,1]), ";", inv_logit(alpha[m,1] + alpha[m,2])))
}

Outcomes = array(0, c(N_id, N_id, M))

for(m in 1:M){
 for(i in 1:N_id){
  for(j in 1:N_id){
     Outcomes[i,j,m] = rbinom(1, size=EE, prob=inv_logit(alpha[m,1] + alpha[m,2]*G$tie_strength[i,j]))  
  } 
 }   
}

Exposure = matrix(EE, nrow=N_id, ncol=N_id)

################################## Fit model

# Outcomes stored as a labeled list
outcome = list(
 Feeding = Outcomes[,,1], 
 Territory = Outcomes[,,2], 
 Courtship = Outcomes[,,3], 
 Nesting = Outcomes[,,4],  
 Attacking = Outcomes[,,5]
)

# Exposure stored as a labeled list
exposure = list(
 Feeding = Exposure, 
 Territory = Exposure, 
 Courtship = Exposure,
 Nesting = Exposure, 
 Attacking = Exposure
)

# Individual data in data-frame
individual = data.frame(Mass=Mass)
groups_f = data.frame(Intercept=factor(groups_1), Merica=factor(groups_2), Quantum=factor(groups_3))

# Row names
name_vec = paste("Individual", 1:N_id)

rownames(Kinship) = colnames(Kinship) = name_vec
rownames(Dominant) = colnames(Dominant) = name_vec
rownames(groups_f) = name_vec
rownames(individual) = name_vec

for(m in 1:M){
rownames(outcome[[m]]) = colnames(outcome[[m]]) = name_vec
rownames(exposure[[m]]) = colnames(exposure[[m]]) = name_vec
}

# Merge data
dat = make_strand_data(
 outcome = outcome,
 block_covariates=groups_f, 
 individual_covariates=individual,
 dyadic_covariates=list(Kinship=Kinship, Dominant=Dominant),
 exposure = exposure,
 outcome_mode="binomial",
 link_mode="logit",
 multiplex = TRUE
)

# Model 
fit = fit_multiplex_model_dimension_reduction(
 data=dat,
 block_regression = ~ Merica + Quantum,
 focal_regression = ~ Mass,
 target_regression = ~ Mass,
 dyad_regression = ~ Kinship * Dominant,
 mode="mcmc",
 return_predicted_network = TRUE,
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 500, 
   iter_sampling = 500, 
   max_treedepth = NULL, 
   adapt_delta = 0.98)
)

res = summarize_strand_results(fit)

##########################################
# Check latent network recovery
par(mfrow=c(1,3))
image(apply(res$samples$predicted_network_sample, 2:3, mean))
image(G$tie_strength)
plot(c(apply(res$samples$predicted_network_sample, 2:3, mean)), c(G$tie_strength))


