################################################# Gaussian models 
# Load libraries
 library(rethinking)
 library(igraph)
 library(STRAND)
 library(ggplot2)

set.seed(1)
# Make data
N_id = 50

# Covariates
Kinship = standardize(rlkjcorr( 1 , N_id , eta=1.5 ))
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
dr_mu = 0 
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
                         error_sigma = 1.0,
                         outcome_mode="gaussian", 
                         link_mode="identity",                 
                         individual_predictors = data.frame(Mass=Mass),
                         dyadic_predictors = dyadic_preds,
                         individual_effects = matrix(sr_effects_1,nrow=2,ncol=1),
                         dyadic_effects = dr_effects_1
                         )      

image(G$network)

################################################################################################ Fit STRAND
# Create the STRAND data object
# Add row and colnames
name_vec = paste("Individual", 1:N_id)
rownames(G$network) = colnames(G$network) = name_vec
rownames(G$samps) = colnames(G$samps) = name_vec
rownames(Kinship) = colnames(Kinship) = name_vec
rownames(Dominant) = colnames(Dominant) = name_vec
rownames(groups_f) = name_vec
rownames(individual) = name_vec

model_dat = make_strand_data(outcome=list(Association = G$network),  
                             block_covariates=groups_f, 
                             individual_covariates=individual, 
                             dyadic_covariates=list(Kinship=Kinship, Dominant=Dominant),  
                             outcome_mode = "gaussian", 
                             link_mode="identity"
                             )


fit = fit_block_plus_social_relations_model(
    data=model_dat,
    block_regression = ~ Merica + Quantum,
    focal_regression = ~ Mass,
    target_regression = ~ Mass,
    dyad_regression = ~ Kinship*Dominant,
    mode="mcmc",
    stan_mcmc_parameters = list(
      chains = 1,
      iter_warmup = 500 ,
      iter_sampling = 500 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )

res = summarize_strand_results(fit)



# Plots
strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE)

# Estimates of dyadic reciprocity correcting for measurement error (i.e., reciprocity in true network)
# dr_rho should be 0.6
strand_VPCs(fit, n_partitions = 5, include_reciprocity = TRUE, mode="cor")

# Estimates of dyadic reciprocity biased by measurement error (i.e., reciprocity of reports)
# dyadic rho in reports is only around 0.4, since measurement error leads to biased reporting. From theory, we know the value should be: 0.6*(1.5^2/(1^2 + 1.5^2)) = 0.41... in this example.
strand_VPCs(fit, n_partitions = 5, include_reciprocity = TRUE, mode="adj")
