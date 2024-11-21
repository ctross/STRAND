###########################################################################
#
#   Bernoulli Analyses: Probit versus Logit
#
#   Some have argued that SRM models require probit links. Here we show
#   that the choice of link function is irrelevant to anything other than
#   the scale of estimates and the run-time of the models.  
#
########################################

# Clear working space
rm(list = ls())

# Load libraries
library(PlvsVltra) # For colors: install_github('ctross/PlvsVltra')
colors = plvs_vltra("dust_storm", rev=FALSE, elements=NULL, show=FALSE)

library(STRAND)
library(ggplot2)

#Load package data
data(Colombia_Data)

###################################################################################### Logit model
# Create the STRAND data object
outcome = list(Friends = Colombia_Data$Friends)

dyad = list(Relatedness = Colombia_Data$Relatedness, 
            Distance = Colombia_Data$Distance
            )

groups = data.frame(Ethnicity = as.factor(Colombia_Data$Individual$Ethnicity), 
                    Sex = as.factor(Colombia_Data$Individual$Sex)
                    )

indiv =  data.frame(Age = Colombia_Data$Individual$Age, 
                    BMI = Colombia_Data$Individual$BMI
                     )

rownames(indiv) = rownames(Colombia_Data$Individual)
rownames(groups) = rownames(Colombia_Data$Individual)

dat_logit = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       outcome_mode = "bernoulli",
                       link_mode = "logit")


# Model
fit_logit = fit_block_plus_social_relations_model(data=dat_logit,
                                            block_regression = ~ Sex + Ethnicity,
                                            focal_regression = ~ Age + BMI,
                                            target_regression = ~ Age + BMI,
                                            dyad_regression = ~ Distance + Relatedness,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 500, iter_sampling = 500,
                                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Summaries
res_logit = summarize_strand_results(fit_logit)

###################################################################################### Probit model
set.seed(4)
# Create the STRAND data object
outcome = list(Friends = Colombia_Data$Friends)

dyad = list(Relatedness = Colombia_Data$Relatedness, 
            Distance = Colombia_Data$Distance
            )

groups = data.frame(Ethnicity = as.factor(Colombia_Data$Individual$Ethnicity), 
                    Sex = as.factor(Colombia_Data$Individual$Sex)
                    )

indiv =  data.frame(Age = Colombia_Data$Individual$Age, 
                    BMI = Colombia_Data$Individual$BMI
                     )

rownames(indiv) = rownames(Colombia_Data$Individual)
rownames(groups) = rownames(Colombia_Data$Individual)

dat_probit = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       outcome_mode = "bernoulli",
                       link_mode = "probit")


# Model. Probit models are more prone to errors of the form:
#   Rejecting initial value:
#   Log probability evaluates to log(0), i.e. negative infinity.
#   Stan can’t start sampling from this initial value.
# So we need to set init<2 in order to have it run.

fit_probit = fit_block_plus_social_relations_model(data=dat_probit,
                                            block_regression = ~ Sex + Ethnicity,
                                            focal_regression = ~ Age + BMI,
                                            target_regression = ~ Age + BMI,
                                            dyad_regression = ~ Distance + Relatedness,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 500, iter_sampling = 500,
                                                                          max_treedepth = NULL, adapt_delta = .98, init = 0.5)
)

# Summaries
res_probit = summarize_strand_results(fit_probit)


######################################################################## VPCs
VPCs_1 = strand_VPCs(fit_logit, n_partitions = 4)
VPCs_2 = strand_VPCs(fit_probit, n_partitions = 4)

df1 = data.frame(do.call(rbind, VPCs_1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Logit"
df1$Variable2 = rep(c("Focal","Target","Dyadic","Error"),1)

df2 = data.frame(do.call(rbind, VPCs_2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Probit"
df2$Variable2 = rep(c("Focal","Target","Dyadic","Error"),1)

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Variable2 = factor(df$Variable2)
df$Variable2 = factor(df$Variable2, levels=rev(c("Focal","Target","Dyadic","Error")))

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Site, color=Site,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + 
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("Logit" = colors[1], "Probit" = colors[3])) + theme(legend.position="bottom")

p

######################################################################## Reciprocity
VPCs_1 = strand_VPCs(fit_logit, n_partitions = 4, include_reciprocity = TRUE)
VPCs_2 = strand_VPCs(fit_probit, n_partitions = 4, include_reciprocity = TRUE)

df1 = data.frame(VPCs_1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Logit"

df2 = data.frame(VPCs_2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Probit"

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Link = df$Site 

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Link, color=Link,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + 
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("Logit" = colors[1], "Probit" = colors[3])) + theme(legend.position="bottom")

p

######################################################################## Block Effects
process_block_parameters(fit_logit, "AFROCOLOMBIAN to AFROCOLOMBIAN", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)
process_block_parameters(fit_logit, "EMBERA to EMBERA", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)
process_block_parameters(fit_logit, "EMBERA to AFROCOLOMBIAN", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)

process_block_parameters(fit_probit, "AFROCOLOMBIAN to AFROCOLOMBIAN", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)
process_block_parameters(fit_probit, "EMBERA to EMBERA", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)
process_block_parameters(fit_probit, "EMBERA to AFROCOLOMBIAN", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)


###################################################################### Plots of covariate effects
vis_1 = strand_caterpillar_plot(res_logit, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="Logit")
vis_2 = strand_caterpillar_plot(res_probit, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="Probit")


df = rbind(vis_1, vis_2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Link = df$Site 

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Link, color=Link,
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
        "lines")) + scale_color_manual(values=c("Logit" = colors[1], "Probit" = colors[3])) + theme(legend.position="bottom")

p

###############################################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
# Check all the relevant parameters
fit_logit$fit$summary()
fit_logit$fit$summary("focal_effects")
fit_logit$fit$summary("target_effects")
fit_logit$fit$summary("block_effects")

################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit_logit$fit$draws(), pars = c("focal_effects[1]","target_effects[1]","sr_L[2,1]","dr_L[2,1]"))

################# Rhat and Effective Samples
# Check all the relevant parameters
fit_probit$fit$summary()
fit_probit$fit$summary("focal_effects")
fit_probit$fit$summary("target_effects")
fit_probit$fit$summary("block_effects")

################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit_probit$fit$draws(), pars = c("focal_effects[1]","target_effects[1]","sr_L[2,1]","dr_L[2,1]"))

