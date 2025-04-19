###########################################
#
#   Binomial Analyses with sampling biases  
#
########################################

# Clear working space
rm(list = ls())

# Load libraries
library(STRAND)
library(ggplot2)
library(igraph)

# Create data
set.seed(420)

V = 1            # One blocking variable
G = 3            # Three categories in this variable
N_id = 85        # Number of bonobos

Group = sample(1:3, N_id, replace=TRUE)
Sex = sample(1:2, N_id, replace=TRUE)
B = matrix(-13, nrow=G, ncol=G)
diag(B) = -9.2  # Block matrix

B[1,3] = -10.1
B[3,2] = -11.9

Coloration = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
Fluffiness = matrix(rnorm(N_id, 0, 1), nrow=N_id, ncol=1)
SizeDiff = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
Relatedness = array(rnorm(N_id*N_id, 0, 1), c(N_id, N_id, 1))
                                                               
A = simulate_sbm_plus_srm_network_with_measurement_bias(N_id = N_id, 
                                                   B=list(B=B),
                                                   V=V,
                                                   groups=data.frame(Group=factor(Group)),
                                                   individual_predictors=Coloration,
                                                   individual_effects=matrix(c(1.9, -1.5),ncol=1, nrow=2),
                                                   dyadic_predictors = SizeDiff,
                                                   dyadic_effects = c(-1.5),
                                                   sr_mu = c(0,0),
                                                   sr_sigma = c(1.4, 0.8),
                                                   sr_rho = 0.6,
                                                   dr_mu = 0,
                                                   dr_sigma = 1.0,
                                                   dr_rho = 0.75,
                                                   exposure_mu = 4.5,
                                                   exposure_sigma = 1.9,
                                                   exposure_max = 20,
                                                   censoring_mu = -4.5,
                                                   censoring_sigma = 0.5,
                                                   exposure_predictors = Coloration,
                                                   censoring_predictors = Coloration,
                                                   exposure_effects = c(-2.7),
                                                   censoring_effects = c(2.1),
                                                   outcome_mode = "binomial",
                                                   link_mode = "logit"
                                                   )

# Plot data
Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
V(Net)$color = c("turquoise4","brown4", "goldenrod3")[A$group_ids$Group]
E(Net)$color = c("grey60","black")[is.mutual(Net)+1]
plot(Net, edge.arrow.size =0.3, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)

# Prep dyadic data
# Add colnames and rownames
animal_names = paste("Bonobo", 1:N_id)
SizeDiff = as.matrix(SizeDiff[,,1])
Relatedness = as.matrix(Relatedness[,,1])
colnames(SizeDiff) = rownames(SizeDiff) = animal_names
colnames(Relatedness) = rownames(Relatedness) = animal_names
colnames(A$net) = rownames(A$net) = animal_names
colnames(A$true_samps) = rownames(A$true_samps) = animal_names

# Make lists
grooming = list(Grooming = A$net)
exposure = list(Grooming = A$true_samps)
dyad = list(SizeDiff = SizeDiff, Relatedness=Relatedness)

# Prep individual data
# Make data frames
block = data.frame(Group = as.factor(Group), Sex = as.factor(ifelse(Sex==1,"M","F")))
indiv =  data.frame(Coloration = Coloration, Fluffiness=Fluffiness)

# Add colnames and rownames
rownames(block) = animal_names
rownames(indiv) = animal_names

# Add in data needed for measurement model supplements
sampled = A$true_exposure
sampled_exposure = rep(20, N_id)
sampled_mask = rep(0, N_id)

detected = A$detected
detected_exposure = A$trials 
detected_mask = rep(0, N_id)


########################################################## Build clean data object
model_dat = make_strand_data(outcome = grooming,
                             individual_covariates = indiv, 
                             block_covariates = block,
                             dyadic_covariates = dyad,
                             m_e_data = list(sampled=sampled, 
                                             sampled_exposure=sampled_exposure, 
                                             sampled_mask=sampled_mask, 
                                             detected=detected, 
                                             detected_exposure=detected_exposure, 
                                             detected_mask=detected_mask),
                             outcome_mode = "binomial",
                             link_mode = "logit",
                             exposure = exposure,
                             imputation=TRUE
                             )

###################################################################
# Introduce some missing data
indiv$Coloration[c(7, 13, 30:35, 42, 46+2, 82:83)] = NA

dyad$Relatedness[c(3, 7, 14, 59, 70), c(14:20)] = NA

grooming$Grooming[c(1:5), c(1:15)] = NA
grooming$Grooming[c(45:50), c(45:50)] = NA

block$Sex[c(17, 39, 22:25, 42, 68, 70, 84)] = NA

sampled[c(16,25, 30:34, 41)] = NA
sampled_exposure[c(16,25, 30:34, 50)] = NA
sampled_mask[c(16,25, 30:34, 50)] = NA

detected[c(17,27, 33:34)] = NA
detected_exposure[c(15, 25, 30:34)] = NA

########################################################## Build data object with NAs
model_dat_nas = make_strand_data(outcome = grooming,
                             individual_covariates = indiv, 
                             block_covariates = block,
                             dyadic_covariates = dyad,
                             m_e_data = list(sampled=sampled, 
                                             sampled_exposure=sampled_exposure, 
                                             sampled_mask=sampled_mask, 
                                             detected=detected, 
                                             detected_exposure=detected_exposure, 
                                             detected_mask=detected_mask),
                             outcome_mode = "binomial",
                             link_mode = "logit",
                             exposure = exposure,
                             imputation=TRUE
                             )

# Imputation model on clean data.                                                         
fit1 =  fit_block_plus_social_relations_model_with_measurement_bias_missings(data=model_dat,
                              block_regression = ~ Group + Sex,
                              focal_regression = ~ Coloration + Fluffiness,
                              target_regression = ~ Coloration + Fluffiness,
                              sampling_regression = ~ Coloration,
                              censoring_regression = ~ Coloration,
                              dyad_regression = ~  SizeDiff + Relatedness,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 500, iter_sampling = 500,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Imputation model on NA affected data. 
fit2 =  fit_block_plus_social_relations_model_with_measurement_bias_missings(data=model_dat_nas,
                              block_regression = ~ Group + Sex,
                              focal_regression = ~ Coloration + Fluffiness,
                              target_regression = ~ Coloration + Fluffiness,
                              sampling_regression = ~ Coloration,
                              censoring_regression = ~ Coloration,
                              dyad_regression = ~  SizeDiff + Relatedness,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 500, iter_sampling = 500,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Base model on clean data.
fit0 =  fit_block_plus_social_relations_model_with_measurement_bias(data=model_dat,
                              block_regression = ~ Group + Sex,
                              focal_regression = ~ Coloration + Fluffiness,
                              target_regression = ~ Coloration + Fluffiness,
                              sampling_regression = ~ Coloration,
                              censoring_regression = ~ Coloration,
                              dyad_regression = ~  SizeDiff + Relatedness,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 500, iter_sampling = 500,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

######################################## Estimates
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
