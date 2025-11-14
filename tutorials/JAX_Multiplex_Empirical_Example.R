#######################################################################################
#
#   Multiplex analyses with JAX 
#
#######################################################################################
# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors
 colors = plvs_vltra("dust_storm", rev=FALSE, elements=NULL, show=FALSE)

 library(STRAND)
 library(stringr)
 library(ggplot2)
 library(psych)

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
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 11, 
   adapt_delta = 0.95,
   cores=4,
   chain_method = "vectorized")
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

p

###################################################### SR correlations
sr_res_numpyro = apply(res_numpyro$samples$srm_model_samples$focal_target_random_effects,2:3, mean)
sr_res_stan = apply(res_stan$samples$srm_model_samples$focal_target_random_effects,2:3, mean)

par(mfrow=c(2,3))
for(q in 1:6)
plot(sr_res_numpyro[,q] ~ sr_res_stan[,q])



###################################### Plots
colors = plvs_vltra("dust_storm", rev=FALSE, elements=c(2,4))
colors = c(colors[1], "grey90", colors[2])

multiplex_plot(fit_numpyro, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = NULL, height=6, width=7, palette=colors)
multiplex_plot(fit_numpyro, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = NULL, height=6, width=7, palette=colors)

multiplex_plot(fit_stan, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE,  save_plot = NULL, height=6, width=7, palette=colors)
multiplex_plot(fit_stan, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE,  save_plot = NULL, height=6, width=7, palette=colors)


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

p


#################################### Run time
fit_numpyro[[6]]
fit_stan[[6]]

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
