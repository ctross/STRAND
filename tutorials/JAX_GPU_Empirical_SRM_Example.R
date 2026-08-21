######################################################################################
#
#   Bernoulli Analyses with JAX using GPU  
#
######################################################################################

# Load libraries
library(STRAND)
library(ggplot2)

# ############################################################################################################################## GPU Setup 
# Setup a GPU back-end for NumPyro. We need these settings in the model call below.
 gpu_settings = list(jax_device_type = "cuda",
                     gpu_versions = c("numpy>=1.27","numpyro>=0.18", "jax[cuda12]>=0.7", "jaxlib>=0.7"), 
                     cusparse_path = NULL)

# If you get a cuSPARSE error, then give path to the correct file, e.g.: 
 gpu_settings$cusparse_path = "/home/cody_ross/.local/share/r-miniconda/pkgs/libcusparse-12.8.2.51-hf6fa245_0/lib/libcusparse.so.12"

 # The next two lines only need to be run once per machine, to build the GPU toolchain in a virtual environment
 #gpu_build_parameters = merge_mcmc_parameters(gpu_settings)
 #create_strand_venv(rebuild = TRUE, verbose = TRUE, mcmc_parameters = gpu_build_parameters) 

 # You may need to restart your session, and then comment out the previous two lines

# ########################################################################################################################################
# Then you can run your models as normal
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
                                            mcmc_parameters = list(
                                               seed = 1,
                                               chains = 2, 
                                               refresh = 1, 
                                               iter_warmup = 1000,                            
                                               iter_sampling = 1000, 
                                               max_treedepth = 12, 
                                               adapt_delta = 0.98,
                                               chain_method = "vectorized",
                                               jax_device_type = gpu_settings$jax_device_type,
                                               gpu_versions = gpu_settings$gpu_versions,
                                               cusparse_path = gpu_settings$cusparse_path)
                                            )

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
                                                                   max_treedepth = 12, adapt_delta = 0.98))


# Summaries
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

ggsave("res_srm.pdf", p, height=8, width=12)


###############################################################################################################
############################### Model fit diagnostics

#################################### Run time
fit_numpyro[[6]]
fit_stan[[6]]

# GPU is much faster

###############################################
# Dont forget to check rhat and ess
# These estimates come from the raw samples, prior to STRAND matching names to coefficients,
# so the parameter names are not helpful.

# Stan
 fit_stan$fit$summary("dyad_effects")
 fit_stan$fit$summary("focal_effects")
 fit_stan$fit$summary("target_effects")
 fit_stan$fit$summary("block_effects")

# JAX
 samples = fit_numpyro$fit$get_samples()
 # Due to some awkward differences in converting from Python to R, we need to reshape by hand
 jax_summary(samples$dyad_effects, n_iter = 1000, n_chains = 2)
 jax_summary(samples$focal_effects, n_iter = 1000, n_chains = 2)
 jax_summary(samples$target_effects, n_iter = 1000, n_chains = 2)
 jax_summary(samples$block_effects, n_iter = 1000, n_chains = 2)


 # Stan often has much nicer fit diagnostics, but much longer run times. As such, its good to make longer runs in NumPyro to compensate.


