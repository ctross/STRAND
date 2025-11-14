######################################################################################
#
#   Bernoulli Analyses with JAX  
#
######################################################################################

# Load libraries
library(STRAND)
library(ggplot2)

#Load package data
data(Colombia_Data)

# Create the STRAND data object
outcome = list(Friends = Colombia_Data$Friends)

dyad = list(Relatedness = standardize_strand(Colombia_Data$Relatedness), 
            Distance = standardize_strand(Colombia_Data$Distance)
            )

groups = data.frame(Ethnicity = as.factor(Colombia_Data$Individual$Ethnicity), 
                    Sex = as.factor(Colombia_Data$Individual$Sex)
                    )

indiv =  data.frame(Age = standardize_strand(Colombia_Data$Individual$Age), 
                    BMI = standardize_strand(Colombia_Data$Individual$BMI)
                     )

rownames(indiv) = rownames(Colombia_Data$Individual)
rownames(groups) = rownames(Colombia_Data$Individual)

dat = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       outcome_mode = "bernoulli",
                       link_mode = "logit")


# Model
fit_numpyro = fit_block_plus_social_relations_model(data=dat,
                                            block_regression = ~ Sex + Ethnicity,
                                            focal_regression = ~ Age + BMI,
                                            target_regression = ~ Age + BMI,
                                            dyad_regression = ~ Distance + Relatedness,
                                            mode="numpyro",
                                            mcmc_parameters = list(chains = 1, seed = 1, refresh = 1,
                                                                   iter_warmup = 1000, iter_sampling = 1000,
                                                                   max_treedepth = 12, adapt_delta = 0.95,
                                                                   cores=4, chain_method = "vectorized"))

# Summaries
res_numpyro = summarize_strand_results(fit_numpyro)

# Model
fit_stan = fit_block_plus_social_relations_model(data=dat,
                                            block_regression = ~ Sex + Ethnicity,
                                            focal_regression = ~ Age + BMI,
                                            target_regression = ~ Age + BMI,
                                            dyad_regression = ~ Distance + Relatedness,
                                            mode="mcmc",
                                            mcmc_parameters = list(chains = 1, refresh = 1,
                                                                   iter_warmup = 1000, iter_sampling = 1000,
                                                                   max_treedepth = 12, adapt_delta = 0.95))


# Summaries
res_numpyro = summarize_strand_results(fit_numpyro)
res_stan = summarize_strand_results(fit_stan)

# Plots
vis_numpyro = strand_caterpillar_plot(res_numpyro, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), 
                                      normalized=TRUE, only_slopes=TRUE, site = "NumPyro", export_as_table = TRUE)

vis_stan = strand_caterpillar_plot(res_stan, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), 
                                      normalized=TRUE, only_slopes=TRUE, site = "Stan", export_as_table = TRUE)


 vis_1 = rbind(vis_stan, vis_numpyro)


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


###############################################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.

################# Rhat and Effective Samples
# Check all the relevant parameters

# Stan
 fit_stan$fit$summary("focal_effects")
 fit_stan$fit$summary("target_effects")
 fit_stan$fit$summary("block_effects")

# NumPyro
 samples = fit_numpyro$fit$get_samples()
 jax_summary(samples$dyad_effects)
 jax_summary(samples$focal_effects)
 jax_summary(samples$target_effects)
 jax_summary(samples$block_effects)
