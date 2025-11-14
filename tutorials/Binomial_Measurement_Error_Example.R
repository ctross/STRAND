##############################################################################################
#
#   Binomial Analyses with sampling biases  
#
##############################################################################################

# Load libraries
library(STRAND)
library(ggplot2)
library(igraph)

# Create data
set.seed(67)

V = 1            # One blocking variable
G = 3            # Three categories in this variable
N_id = 85        # Number of bonobos

Group = sample(1:3, N_id, replace=TRUE)
B = matrix(-12, nrow=G, ncol=G)
diag(B) = -8.2  # Block matrix

B[1,3] = -9.1
B[3,2] = -10.9

Coloration = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
SizeDiff = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
                                                               
A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
                                                   B=list(B=B),
                                                   V=V,
                                                   groups=data.frame(Group=factor(Group)),
                                                   individual_predictors=Coloration,
                                                   individual_effects=matrix(c(1.9, 0.5),ncol=1, nrow=2),
                                                   dyadic_predictors = SizeDiff,
                                                   dyadic_effects = c(-1.5),
                                                   sr_mu = c(0,0),
                                                   sr_sigma = c(1.4, 0.8),
                                                   sr_rho = 0.6,
                                                   dr_mu = 0,
                                                   dr_sigma = 1.0,
                                                   dr_rho = 0.75,
                                                   exposure_mu = 4.5,
                                                   exposure_sigma = 1.9,
                                                   exposure_max = 20,
                                                   censoring_mu = -4.5,
                                                   censoring_sigma = 0.5,
                                                   exposure_predictors = Coloration,
                                                   censoring_predictors = Coloration,
                                                   exposure_effects = c(-2.7),
                                                   censoring_effects = c(2.1),
                                                   outcome_mode = "binomial",
                                                   link_mode = "logit"
                                                   )

# Plot data
image(A$network)

# Prep dyadic data
# Add colnames and rownames
animal_names = paste("Bonobo", 1:N_id)
SizeDiff = as.matrix(SizeDiff[,,1])
colnames(SizeDiff) = rownames(SizeDiff) = animal_names
colnames(A$net) = rownames(A$net) = animal_names
colnames(A$true_samps) = rownames(A$true_samps) = animal_names

# Make lists
grooming = list(Grooming = A$net)
exposure = list(Grooming = A$true_samps)
dyad = list(SizeDiff = SizeDiff)

# Prep individual data
# Make data frames
block = data.frame(Group = as.factor(Group))
indiv =  data.frame(Coloration = Coloration)

# Add colnames and rownames
rownames(block) = animal_names
rownames(indiv) = animal_names

model_dat = make_strand_data(outcome = grooming,
                             individual_covariates = indiv, 
                             block_covariates = block,
                             dyadic_covariates = dyad,
                             outcome_mode = "binomial",
                             link_mode = "logit",
                             exposure = exposure
                             )

# Add in data needed for measurement model supplements
model_dat$sampled = A$true_exposure
model_dat$sampled_exposure = rep(20, N_id)

model_dat$detected = A$detected
model_dat$detected_exposure = A$trials 

# Censoring model with correct specification                                                                
fit =  fit_block_plus_social_relations_model_with_measurement_bias(data=model_dat,
                              block_regression = ~ Group,
                              focal_regression = ~ Coloration,
                              target_regression = ~ Coloration,
                              sampling_regression = ~ Coloration,
                              censoring_regression = ~ Coloration,
                              dyad_regression = ~  SizeDiff,
                              mode="mcmc",
                              mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 500, iter_sampling = 500,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

######################################## Estimates
res = summarize_strand_results(fit)

# Parameter estimates look good, in spite of censoring based on coloration

####################################### Variance partition and reciprocity
VPCs_1 = strand_VPCs(fit, n_partitions = 5, include_reciprocity = TRUE, mode="adj")  # STRAND STYLE, separates error and dyadic effects
VPCs_2 = strand_VPCs(fit, n_partitions = 3, include_reciprocity = TRUE, mode="adj")  # AMEN STYLE, merges error and dyadic effects

VPCs_1
VPCs_2

####################################### Block effects
process_block_parameters(input=fit, focal="2 to 1", base="1 to 1", HPDI=0.9)
process_block_parameters(input=fit, focal="3 to 1", base="1 to 1", HPDI=0.9)
process_block_parameters(input=fit, focal="3 to 2", base="1 to 1", HPDI=0.9)


###############################################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
# Check all the relevant parameters
fit$fit$summary()
fit$fit$summary("focal_effects")
fit$fit$summary("target_effects")
fit$fit$summary("block_effects")

################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit$fit$draws(), pars = c("focal_effects[1]","target_effects[1]","sr_L[2,1]","dr_L[2,1]"))

