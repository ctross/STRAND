#########################################################################################
#
# Multiplex dimension reduction: Gaussian outcomes - Simulated from binomial network data
#
#########################################################################################

# Load libraries
library(reticulate)
library(posterior)
library(igraph)
library(ggplot2)
library(STRAND)

# Make data
set.seed(46+2)
N_id = 50

# Covariates
Kinship = standardize_strand(rlkjcorr( 1 , N_id , eta=1.5 ))
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

#################################################### Simulate a single latent SBM + SRM network
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
                         outcome_mode="binomial", 
                         link_mode="logit",                 
                         individual_predictors = data.frame(Mass=Mass),
                         dyadic_predictors = dyadic_preds,
                         individual_effects = matrix(sr_effects_1,nrow=2,ncol=1),
                         dyadic_effects = dr_effects_1
                         )        

image(G$network)


############################################# Simulate multiplex layers from the latent network
M = 5   # Network layers
EE = 10 # Samples per dyad
alpha = matrix(NA, nrow=M, ncol=2)

# Loadings
alpha[1,] = c(-3, 6)    # Feeding
alpha[2,] = c(0, 6)     # Territory
alpha[3,] = c(-7, 8)    # Courtship
alpha[4,] = c(2, 0.75)  # Nesting
alpha[5,] = c(4, -6)    # Attacking

for(m in 1:M){
print(paste(inv_logit(alpha[m,1]), ";", inv_logit(alpha[m,1] + alpha[m,2])))
}

Outcomes = array(0, c(N_id, N_id, M))

error_sigma = c(1.0, 0.4, 0.1, 0.9, 2.1)

for(m in 1:M){
 for(i in 1:N_id){
  for(j in 1:N_id){
     Outcomes[i,j,m] = rnorm(1, alpha[m,1] + alpha[m,2]*G$tie_strength[i,j], error_sigma[m])   
  } 
 }   
}

Exposure = matrix(EE, nrow=N_id, ncol=N_id)

################################## Fit model
# Outcomes stored as a labeled list
outcome = list(
 Feeding = Outcomes[,,1],    # Note that the network loading on this layer is forced to be positive to help set the direction of the latent axis!
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
 outcome_mode="gaussian", 
 link_mode="identity", 
 multiplex = TRUE
)

######################################################### Model NumPyro Light
fit_numpyro_0 = fit_multiplex_model_dimension_reduction(
 data=dat,
 block_regression = ~ Merica + Quantum,
 focal_regression = ~ Mass,
 target_regression = ~ Mass,
 dyad_regression = ~ Kinship * Dominant,
 expected_ties = 60,                # Tell the model to expect about this many ties to help minimize the chances that the latent factor loads backwards
 mode="numpyro",
 return_predicted_network = FALSE,
 mcmc_parameters = list(
   seed = 1,
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.98,
   cores=4,
   init = 0.01,
   chain_method = "vectorized")
)

# Check fit and loadings
res_numpyro_0 = summarize_strand_results(fit_numpyro_0)
format(object.size(res_numpyro_0), units = "Mb")
fit_numpyro_0[[6]]

######################################################### Model NumPyro, with network reconstruction
fit_numpyro = fit_multiplex_model_dimension_reduction(
 data=dat,
 block_regression = ~ Merica + Quantum,
 focal_regression = ~ Mass,
 target_regression = ~ Mass,
 dyad_regression = ~ Kinship * Dominant,
 expected_ties = 60,         # Tell the model to expect about this many ties to help minimize the chances that the latent factor loads backwards
 mode="numpyro",
 return_predicted_network = TRUE,
 mcmc_parameters = list(
   seed = 1,
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.98,
   cores=4,
   init = 0.01,
   chain_method = "vectorized")
)

# Check fit and loadings
res_numpyro = summarize_strand_results(fit_numpyro)
format(object.size(res_numpyro), units = "Mb")
fit_numpyro[[6]]


######################################################### Model Stan, with network reconstruction
fit_stan = fit_multiplex_model_dimension_reduction(
 data=dat,
 block_regression = ~ Merica + Quantum,
 focal_regression = ~ Mass,
 target_regression = ~ Mass,
 dyad_regression = ~ Kinship * Dominant,
 expected_ties = 60,                  # Tell the model to expect about this many ties to help minimize the chances that the latent factor loads backwards
 mode="mcmc",
 return_predicted_network = TRUE,
 mcmc_parameters = list(
   seed = 1,
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   init = 0.01,
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.98)
)

# Check fit and loadings
res_stan = summarize_strand_results(fit_stan)
format(object.size(res_stan), units = "Mb")
fit_stan[[6]]

#################################################################### Check loadings
# Be careful to check the loadings, since there is a possibility of 
# the model "loading backwards", with all of the predictors having flipped signs.
# Fit with a different seed, or with init=0.
res_numpyro_0$summary_list$'Layer loadings'
res_numpyro$summary_list$'Layer loadings'
res_stan$summary_list$'Layer loadings'

#################################################################### Compare size and runtime
format(object.size(res_numpyro_0), units = "Mb")
format(object.size(res_numpyro), units = "Mb")
format(object.size(res_stan), units = "Mb")

fit_numpyro_0[[6]]
fit_numpyro[[6]]
fit_stan[[6]]

#################################################################### Compare effects of covariates
res_numpyro_tab  = strand_caterpillar_plot(res_numpyro, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=FALSE, site="NumPyro", export_as_table = TRUE)
res_stan_tab = strand_caterpillar_plot(res_stan, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=FALSE, site="Stan", export_as_table = TRUE)

# Identify block effects
# Numpyro
X = mean(res_numpyro_tab$Median[16])
res_numpyro_tab$Median[16] = res_numpyro_tab$Median[16] - X
res_numpyro_tab$LI[16] = res_numpyro_tab$LI[16] - X
res_numpyro_tab$HI[16] = res_numpyro_tab$HI[16] - X

X = mean(res_numpyro_tab$Median[17:25])
res_numpyro_tab$Median[17:25] = res_numpyro_tab$Median[17:25] - X
res_numpyro_tab$LI[17:25] = res_numpyro_tab$LI[17:25] - X
res_numpyro_tab$HI[17:25] = res_numpyro_tab$HI[17:25] - X

X = mean(res_numpyro_tab$Median[26:29])
res_numpyro_tab$Median[26:29] = res_numpyro_tab$Median[26:29] - X
res_numpyro_tab$LI[26:29] = res_numpyro_tab$LI[26:29] - X
res_numpyro_tab$HI[26:29] = res_numpyro_tab$HI[26:29] - X

# Stan
X = mean(res_stan_tab$Median[16])
res_stan_tab$Median[16] = res_stan_tab$Median[16] - X
res_stan_tab$LI[16] = res_stan_tab$LI[16] - X
res_stan_tab$HI[16] = res_stan_tab$HI[16] - X

X = mean(res_stan_tab$Median[17:25])
res_stan_tab$Median[17:25] = res_stan_tab$Median[17:25] - X
res_stan_tab$LI[17:25] = res_stan_tab$LI[17:25] - X
res_stan_tab$HI[17:25] = res_stan_tab$HI[17:25] - X

X = mean(res_stan_tab$Median[26:29])
res_stan_tab$Median[26:29] = res_stan_tab$Median[26:29] - X
res_stan_tab$LI[26:29] = res_stan_tab$LI[26:29] - X
res_stan_tab$HI[26:29] = res_stan_tab$HI[26:29] - X

True = c(sr_sigma[1], sr_effects_1[1],
         sr_sigma[2], sr_effects_1[2],
         dr_sigma, dr_effects_1, sr_rho, dr_rho, error_sigma, c(t(B_1))-mean(B_1), c(t(B_2))-mean(B_2), c(t(B_3))-mean(B_3)
  )

res_true = res_numpyro_tab

res_true$Median = True
res_true$Site = "True"
res_true$LI = NA
res_true$HI = NA

vis_1 = rbind(res_numpyro_tab, res_stan_tab, res_true)
vis_1 = vis_1[which(!vis_1$Variable %in% c("offset, Any to Any")),]

p = ggplot(vis_1, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Site, color=Site)) + geom_linerange(size = 1, position=position_dodge(width=0.32)) + 
        geom_point(size = 2, position=position_dodge(width=0.32)) + facet_grid(Submodel ~ ., scales = "free", space = "free") + geom_hline(aes(yintercept = 0), 
        color = "black", linetype = "dashed") + labs(y = "Regression parameters", x = "") + 
        theme(strip.text.x = element_text(size = 14, face = "bold"), 
          strip.text.y = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14), 
        axis.title.y = element_text(size = 14, face = "bold"), 
                legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.title.x = element_blank()) + theme(strip.text.y = element_text(angle = 360)) + 
        coord_flip() + theme(panel.spacing = unit(1, 
        "lines")) + scale_color_manual(values=c("NumPyro"="chocolate1", "Stan"="deepskyblue2", "True" = "black")) + theme(legend.position="bottom")
p


##########################################
# Check latent network recovery - numpyro
par(mfrow=c(1,3))
image(apply(res_numpyro$samples$predicted_network_sample, 2:3, mean))
image(G$tie_strength)
plot(c(apply(res_numpyro$samples$predicted_network_sample, 2:3, mean)), c(G$tie_strength))

# Check latent network recovery - stan
par(mfrow=c(1,3))
image(apply(res_stan$samples$predicted_network_sample, 2:3, mean))
image(G$tie_strength)
plot(c(apply(res_stan$samples$predicted_network_sample, 2:3, mean)), c(G$tie_strength))

# Compare between NumPyro and Stan
plot(c(apply(res_stan$samples$predicted_network_sample, 2:3, mean)), c(apply(res_numpyro$samples$predicted_network_sample, 2:3, mean)))


###############################################
# Dont forget to check rhat and ess

# Stan
 fit_stan$fit$summary("dyad_effects")
 fit_stan$fit$summary("focal_effects")
 fit_stan$fit$summary("target_effects")
 fit_stan$fit$summary("block_effects")

# JAX
 samples = fit_numpyro$fit$get_samples()
 jax_summary(samples$dyad_effects)
 jax_summary(samples$focal_effects)
 jax_summary(samples$target_effects)
 jax_summary(samples$block_effects)
