#######################################################################################
#
#   Multiplex analyses with JAX using the GPU
#
#######################################################################################

 library(STRAND)
 library(stringr)
 library(ggplot2)
 library(psych)

# ############################################################################################################################## GPU Setup 
# Setup a GPU back-end for NumPyro. We need these settings in the model call below.
 gpu_settings = list(jax_device_type = "cuda",
                     gpu_versions = c("numpy>=1.27","numpyro>=0.18", "jax[cuda12]>=0.7", "jaxlib>=0.7"), 
                     cusparse_path = NULL)

# If you get a cuSPARSE error, then give path to the correct file, e.g.: 
 gpu_settings$cusparse_path = "/home/cody_ross/.local/share/r-miniconda/pkgs/libcusparse-12.8.2.51-hf6fa245_0/lib/libcusparse.so.12"

 # The next two lines only need to be run once per machine, to build the GPU toolchain in a virtual environment
   # gpu_build_parameters = merge_mcmc_parameters(gpu_settings)
   # create_strand_venv(rebuild = TRUE, verbose = TRUE, mcmc_parameters = gpu_build_parameters) 

# ########################################################################################################################################

# Load data
 data(RICH_Data)
 RICH = RICH_Data

# Outcomes stored as a labeled list
 outcome = list(
  Give = RICH$Give, 
  Take = RICH$Take, 
  Reduce = RICH$Reduce
 )

# Dyadic data as a labeled list
 dyad = list(
  Relatedness = standardize_strand(RICH$Relatedness), 
  Friends = RICH$Friends,
  Marriage = RICH$Marriage
 )

# Individual data in data-frame
 RICH$Individual$Age = standardize_strand(RICH$Individual$Age)
 RICH$Individual$Wealth = standardize_strand(RICH$Individual$Wealth)
 ind = RICH$Individual

# Individual blocking measures
 groups = data.frame(
  Ethnicity = as.factor(RICH$Individual$Ethnicity), 
  Sex = as.factor(RICH$Individual$Sex)
 )
 rownames(groups) = rownames(RICH$Individual)

# Merge data
dat = make_strand_data(
 outcome = outcome,
 block_covariates = groups, 
 individual_covariates = ind, 
 dyadic_covariates = dyad,
 outcome_mode="bernoulli",
 link_mode="logit",
 multiplex = TRUE
)

# Model - Full model, all controls
fit_numpyro = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="numpyro",
 mcmc_parameters = list(
   seed = 1,
   chains = 3, 
   refresh = 1, 
   iter_warmup = 2500, 
   iter_sampling = 1000, 
   max_treedepth = 13, 
   adapt_delta = 0.98,
   cores=16,
   chain_method = "vectorized",
   jax_device_type = gpu_settings$jax_device_type,
   gpu_versions = gpu_settings$gpu_versions,
   cusparse_path = gpu_settings$cusparse_path)
)

fit_stan = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 11, 
   adapt_delta = 0.95)
)


#########################################################################
# Model results
res_numpyro = summarize_multiplex_bsrm_results(fit_numpyro)
res_stan = summarize_multiplex_bsrm_results(fit_stan)

#########################################################################
# Merged plot
vis_1 = strand_caterpillar_plot(res_numpyro, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)
vis_2 = strand_caterpillar_plot(res_stan, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)

vis_1$Site = "NumPyro"
vis_2$Site = "Stan"

vis_1$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", vis_1$Variable)
vis_1$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", vis_1$Variable)
vis_1$Variable = gsub("dyadic effects coeffs, ", "", vis_1$Variable)
vis_1$Variable = gsub("dyadic effects ", "", vis_1$Variable)
vis_1$Variable = gsub("focal effects ", "", vis_1$Variable)
vis_1$Variable = gsub("target effects ", "", vis_1$Variable)


vis_2$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", vis_2$Variable)
vis_2$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", vis_2$Variable)
vis_2$Variable = gsub("dyadic effects coeffs, ", "", vis_2$Variable)
vis_2$Variable = gsub("dyadic effects ", "", vis_2$Variable)
vis_2$Variable = gsub("focal effects ", "", vis_2$Variable)
vis_2$Variable = gsub("target effects ", "", vis_2$Variable)

df = rbind(vis_1, vis_2)

df$Outcome = ifelse(str_detect(df$Variable, "Take"), "Take",
             ifelse(str_detect(df$Variable, "Give"), "Give",
             ifelse(str_detect(df$Variable, "Reduce"), "Reduce",
                    NA)))

df$Outcome = factor(df$Outcome)
df$Outcome = factor(df$Outcome, levels=c("Give", "Take", "Reduce"))

df$Variable = gsub("Give - ", "", df$Variable)
df$Variable = gsub("Take - ", "", df$Variable)
df$Variable = gsub("Reduce - ", "", df$Variable)

df$Variable = gsub(" - Give", "", df$Variable)
df$Variable = gsub(" - Take", "", df$Variable)
df$Variable = gsub(" - Reduce", "", df$Variable)

df$Variable = gsub("sd", "SD", df$Variable)

df$Variable = gsub("FoodInsecure", "Food Insecure", df$Variable)
df$Variable = gsub("Wealth", "Log Wealth", df$Variable)

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Site, color=Site,
        ymin = LI, ymax = HI)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(Submodel ~ 
        Outcome, scales = "free", space = "free") + ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("NumPyro"="chocolate1", "Stan"="deepskyblue2", "True" = "black")) + theme(legend.position="bottom")

ggsave("res.pdf",p, height=8, width=12)


###################################### Plots
####################################### VPCs
VPCs_1 = strand_VPCs(fit_numpyro, n_partitions = 5)
VPCs_2 = strand_VPCs(fit_stan, n_partitions = 5)

df1 = data.frame(do.call(rbind, VPCs_1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "NumPyro"
df1$Submodel = rep(c("Give","Take","Reduce"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic signal","Dyadic noise + Error"),3)

df2 = data.frame(do.call(rbind, VPCs_2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Stan"
df2$Submodel = rep(c("Give","Take","Reduce"),each=4)
df2$Variable2 = rep(c("Focal","Target","Dyadic signal","Dyadic noise + Error"),3)

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Give", "Take", "Reduce"))

df$Variable2 = factor(df$Variable2)


p = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Site, color=Site,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("NumPyro"="chocolate1", "Stan"="deepskyblue2", "True" = "black"))  + theme(legend.position="bottom")

ggsave("res2.pdf",p, height=8, width=12)


#################################### Run time
fit_numpyro[[6]]
fit_stan[[6]]

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
 jax_summary(samples$dyad_effects, n_iter = 1000, n_chains = 3)
 jax_summary(samples$focal_effects, n_iter = 1000, n_chains = 3)
 jax_summary(samples$target_effects, n_iter = 1000, n_chains = 3)
 jax_summary(samples$block_effects, n_iter = 1000, n_chains = 3)

