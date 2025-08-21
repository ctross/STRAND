###################################### Load packages
# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors
 color_set = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

 library(STRAND)
 library(stringr)
 library(ggplot2)
 library(psych)
 library(rethinking)
 library(Matrix)
 library(igraph)

###################################### First simulate a network
# Clear working space
rm(list = ls())
set.seed(46+2)

# Make data with this sample size
N_id = 70
labels = paste("Ind", 1:N_id)

# Covariates
Kinship = STRAND::standardize(rlkjcorr( 1 , N_id , eta=1.5 ) )               
Dominant = ceiling(rlkjcorr( 1 , N_id , eta=1.5 ) - 0.1) 
rownames(Kinship) = colnames(Kinship) = labels
rownames(Dominant) = colnames(Dominant) = labels

Mass = rbern(N_id, 0.4)
Gold = rnorm(N_id, 0, 1)

# Organize data into arrays
dyadic_preds = array(NA,c(N_id,N_id,2))

dyadic_preds[,,1] = Kinship
dyadic_preds[,,2] = Dominant

# Set effect sizes
sr_mu = c(0,0)  
sr_sigma = c(1.8, 1.8)           
sr_rho = -0.7                  
dr_mu = 0 
dr_sigma = 3.5
dr_rho = -0.8                  
sr_effects_1 = c(1.9, 1.9)     
sr_effects_2 = c(-1.1, -1.1)     
dr_effects_1 = c(-1.9, 2.1)       

# Block structure
group_probs_block_size = c(0.25, c(0.5, 0.5)*(1-0.25))

B_1 = matrix(-12,nrow=1,ncol=1)
B_2 = matrix(rnorm(9,0,3),nrow=3,ncol=3)

diag(B_2) = diag(B_2) + 4.5

B = list(B_1, B_2)
 
groups_1 = rep("Any", N_id) 
groups_2 = sample(c("Red","White","Blue"), size=N_id, replace=TRUE, prob=group_probs_block_size)


groups = data.frame(Intercept=as.numeric(factor(groups_1)), Merica=as.numeric(factor(groups_2)))
groups_f = data.frame(Intercept=factor(groups_1), Merica=factor(groups_2))
individual = data.frame(Mass=Mass, Gold=Gold)
rownames(individual) = labels
rownames(groups_f) = labels

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
                         individual_predictors = data.frame(Mass = Mass, Gold=Gold),
                         dyadic_predictors = dyadic_preds,
                         individual_effects = rbind(sr_effects_1, sr_effects_2),
                         dyadic_effects = dr_effects_1
                         )  


image(G$network/G$samps)

Net = graph_from_adjacency_matrix(G$network/G$samps, mode = c("directed"))
V(Net)$color = c("turquoise4","gray13", "goldenrod3")[G$group_ids$Merica]

plot(Net, edge.arrow.size =0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)

#################################################### Now organize data for diffusion analysis
 long_dat = NULL
 T = 100

# In this example, the network and data are constant over timesteps, but this need not be the case always
 for(t in 1:T){
  dyadic_preds = list(Kinship = Kinship, Dominant = Dominant)

  outcome = G$network
  exposure = G$samps
  rownames(outcome) = colnames(outcome) = labels
  rownames(exposure) = colnames(exposure) = labels

  long_dat[[t]] = make_strand_data(
                       outcome = list(Tolerance=outcome),
                       exposure = list(Tolerance=exposure),
                       block_covariates = groups_f, 
                       individual_covariates = individual, 
                       dyadic_covariates = dyadic_preds,
                       diffusion_outcome = rbinom(N_id, size=1, prob=0.1), 
                       diffusion_exposure = NULL, 
                       diffusion_mask = NULL,
                       longitudinal = TRUE,
                       link_mode = "logit", 
                       outcome_mode="binomial")
 }

 names(long_dat) = paste("Time", c(1:T))



diff_dat = simulate_diffusion(long_dat,
                   base_rates = c(-6, -8),
                   individual_focal_regression = ~ Gold + Mass,
                   individual_focal_parameters = c(0.8, -0.97),
                   social_focal_regression = ~ Gold + Mass,
                   social_focal_parameters = c(-1.4, 1.9),
                   social_target_regression = ~ Gold + Mass,
                   social_target_parameters = c(-0.9, 1.7),
                   social_dyad_regression = ~ Kinship + Dominant,
                   social_dyad_parameters = c(1.9, 1.8),
                   ces_parameters = list(alpha = 0.05, sigma = 100, eta = 1)
                          )


plot_diffusion(diff_dat)


################################################################ Fit using a point estimate of network
fit = fit_NBDA_model(long_data = diff_dat,
    individual_focal_regression = ~ Gold + Mass,
    social_focal_regression = ~ Gold + Mass,
    social_target_regression = ~ Gold + Mass,
    social_dyad_regression = ~ Kinship + Dominant,
    social_block_regression = ~ 1,
    network_treatment = "point",
    ces_settings = "es_inf_rts_1",
    mode="mcmc",
    stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500,
                                iter_sampling = 500, max_treedepth = 12, adapt_delta = 0.95)
    )


res = summarize_strand_results(fit)
