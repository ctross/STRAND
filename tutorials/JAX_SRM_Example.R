#########################################################################
#
#   Single-layer network analyses with JAX 
#
#########################################################################

# Load libraries
 library(reticulate)
 library(posterior)
 library(igraph)
 library(STRAND)
 library(ggplot2)

 set.seed(1)

# Make data
 N_id = 110

# Covariates
 Kinship = standardize_strand(rlkjcorr( 1 , N_id , eta=1.5 ))
 Dominant = matrix(rbern(N_id*N_id), 0.2, nrow=N_id, ncol=N_id)
 Mass = rbern(N_id, 0.4)
 Age = rnorm(N_id, 0, 1)
 Love = rbern(N_id, 0.4)
 Fire = rbern(N_id, 0.6)

# Organize into list
 dyadic_preds = array(NA,c(N_id,N_id,3))

 dyadic_preds[,,1] = Kinship
 dyadic_preds[,,2] = Dominant
 dyadic_preds[,,3] = Kinship*Dominant

# Set effect sizes
 sr_mu = c(0,0)  
 sr_sigma = c(2.2, 1.7) 
 sr_rho = -0.7
 dr_mu = 0 
 dr_sigma = 3.5
 dr_rho= 0.8
 sr_effects_1 = c(2, 1.5)
 sr_effects_2 = c(-1.4, 1.0)
 sr_effects_3 = c(1.4, -1.3)
 sr_effects_4 = c(-0.9, 0)
 dr_effects_1 = c(1.5, -1.9, 2.4)

# Block structure
 group_probs_block_size = c(0.25, c(0.25, 0.25)*(1-0.25))

 B_1 = matrix(-8,nrow=1,ncol=1)
 B_2 = matrix(rnorm(9,0,3),nrow=3,ncol=3)
 B_3 = matrix(rnorm(4,0,3),nrow=2,ncol=2)

 diag(B_2) = diag(B_2) + 3.5
 diag(B_3) = diag(B_3) + 3.5

 B = list(B_1, B_2) #, B_3)
 
 groups_1 = rep("Any",N_id) 
 groups_2 = sample( c("Red","White","Blue") , size=N_id , replace=TRUE , prob=group_probs_block_size )
 groups_3 = sample( c("Strange", "Charm") , size=N_id , replace=TRUE , prob=c(0.5,0.5) )

 groups = data.frame(Intercept=as.numeric(factor(groups_1)), Merica=as.numeric(factor(groups_2)))
 groups_f = data.frame(Intercept=factor(groups_1), Merica=factor(groups_2))
 individual = data.frame(Mass = Mass, Age = Age, Love = Love, Fire = Fire)

#################################################### Simulate SBM + SRM network
 G = simulate_sbm_plus_srm_network(N_id = N_id, 
                         B = B, 
                         V = 2,
                         groups=groups,                  
                         sr_mu = sr_mu,  
                         sr_sigma = sr_sigma, 
                         sr_rho = sr_rho,
                         dr_mu = dr_mu,  
                         dr_sigma = dr_sigma, 
                         dr_rho = dr_rho,
                         error_sigma = 1.0,
                         outcome_mode="bernoulli", 
                         link_mode="logit",                 
                         individual_predictors = individual,
                         dyadic_predictors = dyadic_preds,
                         individual_effects = cbind(sr_effects_1, sr_effects_2, sr_effects_3, sr_effects_4),
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
                             outcome_mode="bernoulli", 
                             link_mode="logit",
                             exposure=list(Association = G$samps)
                             )

####### First test the lighter version which doesnt save dyadic samples or network predictions which are big N x N X N_samples objects
 fit_numpyro_0 = fit_block_plus_social_relations_model(
    data=model_dat,
    block_regression = ~ Merica,
    focal_regression = ~ Mass + Age + Love + Fire,
    target_regression = ~ Mass + Age + Love + Fire,
    dyad_regression = ~ Kinship*Dominant,
    return_predicted_network=FALSE,
    mode="numpyro",
    mcmc_parameters = list(
      seed = 1,  
      chains = 1,
      iter_warmup = 1000 ,
      iter_sampling = 1000 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.95,
      cores=2,
      chain_method = "vectorized")
  )

 res_numpyro_0 = summarize_strand_results(fit_numpyro_0)
 format(object.size(res_numpyro_0), units = "Mb")
 fit_numpyro_0[[6]]

###### Now compare the full model agaist Stan
 fit_numpyro = fit_block_plus_social_relations_model(
    data=model_dat,
    block_regression = ~ Merica,
    focal_regression = ~ Mass + Age + Love + Fire,
    target_regression = ~ Mass + Age + Love + Fire,
    dyad_regression = ~ Kinship*Dominant,
    return_predicted_network=TRUE,
    mode="numpyro",
    mcmc_parameters = list(
      seed = 1,  
      chains = 1,
      iter_warmup = 1000 ,
      iter_sampling = 1000 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.95,
      cores=2,
      chain_method = "vectorized")
  )

 res_numpyro = summarize_strand_results(fit_numpyro)
 format(object.size(res_numpyro), units = "Mb")
 fit_numpyro[[6]]

 fit_stan = fit_block_plus_social_relations_model(
    data=model_dat,
    block_regression = ~ Merica,
    focal_regression = ~ Mass + Age + Love + Fire,
    target_regression = ~ Mass + Age + Love + Fire,
    dyad_regression = ~ Kinship*Dominant,
    return_predicted_network=TRUE,
    mode="mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1000 ,
      iter_sampling = 1000 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )

 res_stan = summarize_strand_results(fit_stan)
 format(object.size(res_stan), units = "Mb")
 fit_stan[[6]]

#################################################################### Compare effects of covariates
 res_numpyro_0_tab  = strand_caterpillar_plot(res_numpyro_0, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=FALSE, site="NumPyro_0", export_as_table = TRUE)
 res_numpyro_tab  = strand_caterpillar_plot(res_numpyro, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=FALSE, site="NumPyro", export_as_table = TRUE)
 res_stan_tab = strand_caterpillar_plot(res_stan, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=FALSE, site="Stan", export_as_table = TRUE)

 Offset = res_numpyro_0_tab$Median[18]
 res_numpyro_0_tab$Median[18:27] = res_numpyro_0_tab$Median[18:27] + Offset
 res_numpyro_0_tab$LI[18:27] = res_numpyro_0_tab$LI[18:27] + Offset
 res_numpyro_0_tab$HI[18:27] = res_numpyro_0_tab$HI[18:27] + Offset

 Offset = res_numpyro_tab$Median[18]
 res_numpyro_tab$Median[18:27] = res_numpyro_tab$Median[18:27] + Offset
 res_numpyro_tab$LI[18:27] = res_numpyro_tab$LI[18:27] + Offset
 res_numpyro_tab$HI[18:27] = res_numpyro_tab$HI[18:27] + Offset

 Offset = res_stan_tab$Median[18]
 res_stan_tab$Median[18:27] = res_stan_tab$Median[18:27] + Offset
 res_stan_tab$LI[18:27] = res_stan_tab$LI[18:27] + Offset
 res_stan_tab$HI[18:27] = res_stan_tab$HI[18:27] + Offset

 True = c(sr_sigma[1], sr_effects_1[1], sr_effects_2[1], sr_effects_3[1], sr_effects_4[1],
         sr_sigma[2], sr_effects_1[2], sr_effects_2[2], sr_effects_3[2], sr_effects_4[2],
         dr_sigma, dr_effects_1, sr_rho, dr_rho, NA, c(t(B_1))+B_1[1,1], c(t(B_2))+B_1[1,1]
  )

 res_true = res_numpyro_tab

 res_true$Median = True
 res_true$Site = "True"
 res_true$LI = NA
 res_true$HI = NA

 vis_1 = rbind(res_numpyro_0_tab, res_numpyro_tab,res_stan_tab, res_true)
 vis_1 = vis_1[which(!vis_1$Variable %in% c("error sd", "offset, Any to Any")),]

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
        "lines")) + scale_color_manual(values=c("NumPyro_0"="chocolate4", "NumPyro"="chocolate1", "Stan"="deepskyblue2", "True" = "black")) + theme(legend.position="bottom")
 p

############################## Check Random effects
 par(mfrow=c(1,4))

################################### SR fits
 sr_res_numpyro = apply(res_numpyro$samples$srm_model_samples$focal_target_random_effects,2:3, mean)
 sr_res_stan = apply(res_stan$samples$srm_model_samples$focal_target_random_effects,2:3, mean)
 plot(sr_res_numpyro[,1] ~ sr_res_stan[,1])
 plot(sr_res_numpyro[,2] ~ sr_res_stan[,2])

################################### DR fits
 dr_res_numpyro = apply(res_numpyro$samples$srm_model_samples$dyadic_random_effects, 2:3, mean)
 dr_res_stan = apply(res_stan$samples$srm_model_samples$dyadic_random_effects, 2:3, mean)

 plot(c(dr_res_numpyro) ~ c(dr_res_stan))

################################### Predicted network fits
 p_res_numpyro = apply(res_numpyro$samples$predicted_network_sample, 2:3, mean)
 p_res_stan = apply(res_stan$samples$predicted_network_sample, 2:3, mean)

 plot(c(p_res_numpyro) ~ c(p_res_stan))

################################### VPCs
 strand_VPCs(fit_numpyro, include_reciprocity=TRUE)
 strand_VPCs(fit_stan, include_reciprocity=TRUE)

#################################### Run time
 fit_numpyro_0[[6]]
 fit_numpyro[[6]]
 fit_stan[[6]]

##################### Size
 format(object.size(res_numpyro_0), units = "Mb")
 format(object.size(res_numpyro), units = "Mb")
 format(object.size(res_stan), units = "Mb")

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



