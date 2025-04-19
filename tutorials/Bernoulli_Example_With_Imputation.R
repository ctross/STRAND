###########################################
#
#   Bernoulli Analyses - Imputation models 
#
###########################################

# Clear working space
rm(list = ls())

# Load libraries
library(STRAND)
library(ggplot2)

#Load package data
data(Colombia_Data)

# Create the STRAND data object
outcome = list(Friends = Colombia_Data$Friends)

dyad = list(Relatedness = standardize(Colombia_Data$Relatedness), 
            Distance = standardize(Colombia_Data$Distance)
            )

groups = data.frame(Ethnicity = as.factor(Colombia_Data$Individual$Ethnicity), 
                    Sex = as.factor(Colombia_Data$Individual$Sex)
                    )

indiv =  data.frame(Age = standardize(Colombia_Data$Individual$Age), 
                    BMI = standardize(Colombia_Data$Individual$BMI)
                     )

rownames(indiv) = rownames(Colombia_Data$Individual)
rownames(groups) = rownames(Colombia_Data$Individual)

### A clean data set
dat0 = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       outcome_mode = "bernoulli",
                       link_mode = "logit")

###################################################################
# Introduce some missing data
indiv$Age[c(7, 13, 30:35, 42, 46+2, 82:84, 91)] = NA
dyad$Relatedness[c(3, 7, 14, 59, 70), c(14:20)] = NA

outcome$Friends[c(1:5), c(1:15)] = NA
outcome$Friends[c(45:50), c(45:50)] = NA

groups$Sex[c(17, 39, 22:25, 42, 68, 70, 87)] = NA

############# By default STRAND will throw an error if there are NAs
dat = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       outcome_mode = "bernoulli",
                       link_mode = "logit")

############# Add: imputation = TRUE to turn off error checks for NAs
dat = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       outcome_mode = "bernoulli",
                       link_mode = "logit",
                       imputation = TRUE)

###################################################################
# Now fit new model which deals with missings
fit1 = fit_block_plus_social_relations_model_missings(data=dat,
                                            block_regression = ~ Sex + Ethnicity,
                                            focal_regression = ~ Age + BMI,
                                            target_regression = ~ Age + BMI,
                                            dyad_regression = ~ Distance + Relatedness,
                                            mode="mcmc",
                                            return_predicted_network=FALSE,
                                            priors=NULL,
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 500, iter_sampling = 500,
                                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Now fit new model on the true data without missings
fit2 = fit_block_plus_social_relations_model_missings(data=dat0,
                                            block_regression = ~ Sex + Ethnicity,
                                            focal_regression = ~ Age + BMI,
                                            target_regression = ~ Age + BMI,
                                            dyad_regression = ~ Distance + Relatedness,
                                            mode="mcmc",
                                            return_predicted_network=FALSE,
                                            priors=NULL,
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 500, iter_sampling = 500,
                                                                          max_treedepth = NULL, adapt_delta = .98)
)

# And the old basic model on the true data
fit0 = fit_block_plus_social_relations_model(data=dat0,
                                            block_regression = ~ Sex + Ethnicity,
                                            focal_regression = ~ Age + BMI,
                                            target_regression = ~ Age + BMI,
                                            dyad_regression = ~ Distance + Relatedness,
                                            mode="mcmc",
                                            return_predicted_network=FALSE,
                                            priors=NULL,
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 500, iter_sampling = 500,
                                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Summaries
res1 = summarize_strand_results(fit1)
res2 = summarize_strand_results(fit2)
res0 = summarize_strand_results(fit0)


# Plots
vis_1 = strand_caterpillar_plot(res1, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-Imputation")
vis_2 = strand_caterpillar_plot(res2, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-No-Imputation")
vis_0 = strand_caterpillar_plot(res0, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="Old-No-Imputation")


df = rbind(vis_1, vis_2, vis_0)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Method = df$Site 

p1 = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Method, color=Method,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(linewidth = 1, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("New-Imputation" = "darkred","Old-No-Imputation" = "black", "New-No-Imputation" = "goldenrod")) + theme(legend.position="bottom")

p1


################################################################# Reciprocity terms
# Simple correlations
cor1 = strand_VPCs(fit1, n_partitions = 4, include_reciprocity = TRUE, mode="cor")
cor2 = strand_VPCs(fit2, n_partitions = 4, include_reciprocity = TRUE, mode="cor")
cor0 = strand_VPCs(fit0, n_partitions = 4, include_reciprocity = TRUE, mode="cor")

df1 = data.frame(cor1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "New-Imputation"
df1$Submodel = rep(c("Generalized","Dyadic"),each=1)

df2 = data.frame(cor2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "New-No-Imputation"
df2$Submodel = rep(c("Generalized","Dyadic"),each=1)

df0 = data.frame(cor0[[3]])
colnames(df0) = c("Variable", "Median", "L", "H", "Mean", "SD")
df0$Site = "Old-No-Imputation"
df0$Submodel = rep(c("Generalized","Dyadic"),each=1)

# Adjusted correlations apply only for dyadic reciprocity. Adjusted correlations give the correlation in dyadic random effects+error. Unadjusted give correlation in dyadic random effects.

df = rbind(df1, df2, df0)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Method = df$Site 

p3 = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Method, color=Method,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(linewidth = 1, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("New-Imputation" = "darkred","Old-No-Imputation" = "black", "New-No-Imputation" = "goldenrod")) + theme(legend.position="bottom")

p3

######################################################################## VPCs
vpc1 = strand_VPCs(fit1, n_partitions = 4)
vpc2 = strand_VPCs(fit2, n_partitions = 4) 
vpc0 = strand_VPCs(fit0, n_partitions = 4) 

df1 = data.frame(do.call(rbind, vpc1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "New-Imputation"
df1$Submodel = rep(c("Friendship"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic","Error"),1)

df2 = data.frame(do.call(rbind, vpc2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "New-No-Imputation"
df2$Submodel = rep(c("Friendship"),each=4)
df2$Variable2 = rep(c("Focal","Target","Dyadic","Error"),1)

df0 = data.frame(do.call(rbind, vpc0[[2]]))
colnames(df0) = c("Variable", "Median", "L", "H", "Mean", "SD")
df0$Site = "Old-No-Imputation"
df0$Submodel = rep(c("Friendship"),each=4)
df0$Variable2 = rep(c("Focal","Target","Dyadic","Error"),1)

df = rbind(df1, df2, df0)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Friendship"))

df$Variable2 = factor(df$Variable2)
df$Variable2 = factor(df$Variable2, levels=rev(c("Focal","Target","Dyadic","Error","Dyadic+Error")))

df$Type = df$Site 

p4 = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Type, color=Type,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(linewidth = 1, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("New-Imputation" = "darkred","Old-No-Imputation" = "black", "New-No-Imputation" = "goldenrod")) + theme(legend.position="bottom")

p4

