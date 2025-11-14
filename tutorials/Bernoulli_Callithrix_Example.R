##############################################################################################
#
#   Understanding the predicted network
#
##############################################################################################
 
# When users fit a block_plus_social_relations_model, or similar models, they can export a 
# predicted network of ties weights. These tie weight are returned, even for dyads that are
# masked out of the regression. The predictions use the supplied covariate data. Sometimes,
# users will want the full set of predictions, other times user will want to post-process the
# predicted network before using it. Lets look at two examples.

# Load libraries
library(STRAND)
library(ggplot2)

#############################################################################################################
####### Case 1: Here we have a case when different groups of animals are kept in different cages
# Between-cage ties are impossible, and can be treated as strucutural zeros and masked out of the regression

# Load package data
data(Callithrix_Data)

# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[1]]$Aggressed)
mask = list(Aggressed = Callithrix_Data[[1]]$NoOpportunity)               # 0s are possible ties, 1s are censored/impossible. The mask will remove such terms from the likelihood of the model.
dyad = list(RankDiff = standardize_strand(Callithrix_Data[[1]]$RankDiff)) # rank distances can be computed, even between individuals in different cages, but arent necessarily meaningful
indiv =  Callithrix_Data[[1]]$Individual

model_dat_cm = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             mask=mask,                    # if mask != NULL, then it will be used by any model that reads this data.
                             outcome_mode = "bernoulli",
                             link_mode = "logit"
                             )
# Model
fit_cm_1 =  fit_block_plus_social_relations_model(data=model_dat_cm,
                              block_regression = ~ 1,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff,
                              return_predicted_network =TRUE,
                              mode="mcmc",
                              mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = 11, adapt_delta = 0.95)
)

# Get results. Data behind mask dont affect parameter estimates.
res_cm_mask = summarize_strand_results(fit_cm_1)

# Note that the data behind the mask *do* influence the predicted ties! Effects of covariates behind the mask are extrapolated to predict for the masked set.
pred_network_1a = apply(res_cm_mask$samples$predicted_network_sample, c(2,3), mean)
pred_network_1b = pred_network_1a*(1-mask[[1]])
par(mfrow=c(2,2))
image(pred_network_1a)  # The off-diagonal blocks are meaningless in this case
image(pred_network_1b)  # Remove those predictions using the mask

# Now we can plot binarized network graphs if we want
binarized_network_1a = graph_from_adjacency_matrix(ifelse(pred_network_1a>0.5,1,0), mode="directed") 
binarized_network_1b = graph_from_adjacency_matrix(ifelse(pred_network_1b>0.5,1,0), mode="directed")  

V(binarized_network_1a)$color = c(rep("darkred",5), rep("slateblue",6), rep("grey50",5), rep("orange4",4))
V(binarized_network_1b)$color = c(rep("darkred",5), rep("slateblue",6), rep("grey50",5), rep("orange4",4))

plot(binarized_network_1a, edge.arrow.size = 0.5)    # Not realistic
plot(binarized_network_1b, edge.arrow.size = 0.5)    # Realistic

################################################################ Condition, S++, high food
# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[2]]$Aggressed)
mask = list(Aggressed = Callithrix_Data[[2]]$NoOpportunity) # 0s are possible ties, 1s are censored/impossible. The mask will remove such terms from the model.
dyad = list(RankDiff = standardize_strand(Callithrix_Data[[2]]$RankDiff))
indiv =  Callithrix_Data[[2]]$Individual

model_dat_spp = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             mask=mask,
                             outcome_mode = "bernoulli",
                             link_mode = "logit"
                             )
# Model
fit_spp_1 =  fit_block_plus_social_relations_model(data=model_dat_spp,
                              block_regression = ~ 1,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff,
                              mode="mcmc",
                              mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = 11, adapt_delta = 0.95)
)

# Get results.
res_spp_mask = summarize_strand_results(fit_spp_1)

############################################################################################## 
## Interaction method tests
##############################################################################################

################################################################ Condition, C-, low food
# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[1]]$Aggressed)

dyad = list(RankDiff = standardize_strand(Callithrix_Data[[1]]$RankDiff),
            NoOpportunity = Callithrix_Data[[1]]$NoOpportunity)

indiv = Callithrix_Data[[1]]$Individual

model_dat_cm = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             outcome_mode = "bernoulli",
                             link_mode = "logit"
                             )
# Model
fit_cm_2 =  fit_block_plus_social_relations_model(data=model_dat_cm,
                              block_regression = ~ 1,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff * NoOpportunity,
                              mode="mcmc",
                              mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = 11, adapt_delta = 0.95)
)

# Get results.
res_cm_interaction = summarize_strand_results(fit_cm_2)


################################################################ Condition, S++, high food
# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[2]]$Aggressed)

dyad = list(RankDiff = standardize_strand(Callithrix_Data[[2]]$RankDiff),
            NoOpportunity = Callithrix_Data[[2]]$NoOpportunity)

indiv = Callithrix_Data[[2]]$Individual

model_dat_spp = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             outcome_mode = "bernoulli",
                             link_mode = "logit"
                             )
# Model
fit_spp_2 =  fit_block_plus_social_relations_model(data=model_dat_spp,
                              block_regression = ~ 1,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff * NoOpportunity,
                              mode="mcmc",
                              mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = 11, adapt_delta = 0.95)
)

# Get results.
res_spp_interaction = summarize_strand_results(fit_spp_2)

##############################################################################################
#################################################################### Compare effects of covariates
res_cm_mask_tab  = strand_caterpillar_plot(res_cm_mask, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="C-", export_as_table = TRUE)
res_spp_mask_tab = strand_caterpillar_plot(res_spp_mask, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="S++", export_as_table = TRUE)

res_cm_interaction_tab  = strand_caterpillar_plot(res_cm_interaction, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="C-", export_as_table = TRUE)
res_spp_interaction_tab = strand_caterpillar_plot(res_spp_interaction, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="S++", export_as_table = TRUE)

vis_1 = rbind(res_cm_mask_tab, res_spp_mask_tab, res_cm_interaction_tab, res_spp_interaction_tab)
vis_1 = vis_1[which(!vis_1$Variable %in% c("error sd")),]

vis_1$Variable =
c("Female", "Female", 
"RankDiff", "Intercept", 
"Female", "Female", 
"RankDiff", "Intercept",
"Female", "Female", 
"RankDiff", "NoOpportunity", 
"RankDiff:NoOpportunity", "Intercept", 
"Female", "Female", 
"RankDiff", "NoOpportunity", 
"RankDiff:NoOpportunity", "Intercept"
)

vis_1$Type =
c("Mask", "Mask", 
"Mask", "Mask", 
"Mask", "Mask", 
"Mask", "Mask",
"Interaction", "Interaction", 
"Interaction", "Interaction", 
"Interaction", "Interaction", 
"Interaction", "Interaction", 
"Interaction", "Interaction", 
"Interaction", "Interaction"
)

vis_1$Condition = vis_1$Site

vis_1$Type2 = factor(paste(vis_1$Condition, vis_1$Type))

p = ggplot(vis_1, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Type2, color=Condition, linetype=Type)) + geom_linerange(linewidth = 1, position=position_dodge(width=0.32)) + 
        geom_point(size = 2, position=position_dodge(width=0.32)) + facet_grid(Submodel ~ ., scales = "free", space = "free") + geom_hline(aes(yintercept = 0), 
        color = "black", linetype = "dashed") + labs(y = "Regression parameters", x = "") + 
        theme(strip.text.x = element_text(size = 14, face = "bold"), 
          strip.text.y = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14), 
        axis.title.y = element_text(size = 14, face = "bold"), 
                legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.title.x = element_blank()) + theme(strip.text.y = element_text(angle = 360)) + 
        coord_flip() + theme(panel.spacing = unit(1, 
        "lines")) + scale_color_manual(values=c("C-"="#697f4f", "S++"="black")) + theme(legend.position="bottom")
p

##############################################################################################
################################################################# Compare random effect structure
vis_2_cm_mask = strand_caterpillar_plot(res_cm_mask, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE,  only_technicals=TRUE, only_slopes=FALSE, site="C-", export_as_table = TRUE)
vis_2_spp_mask = strand_caterpillar_plot(res_spp_mask, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE,  only_technicals=TRUE, only_slopes=FALSE, site="S++", export_as_table = TRUE)

vis_2_cm_interaction = strand_caterpillar_plot(res_cm_interaction, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE,  only_technicals=TRUE, only_slopes=FALSE, site="C-", export_as_table = TRUE)
vis_2_spp_interaction = strand_caterpillar_plot(res_spp_interaction, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE,  only_technicals=TRUE, only_slopes=FALSE, site="S++", export_as_table = TRUE)

vis_2 = rbind(vis_2_cm_mask, vis_2_spp_mask, vis_2_cm_interaction, vis_2_spp_interaction)
vis_2$Condition = vis_2$Site

vis_2$Type =
c("Mask", "Mask", 
"Mask", "Mask", 
"Mask", "Mask", 
"Mask", "Mask",
"Mask", "Mask", 
"Interaction", "Interaction", 
"Interaction", "Interaction", 
"Interaction", "Interaction", 
"Interaction", "Interaction", 
"Interaction", "Interaction"
)

vis_2$Type2 = factor(paste(vis_2$Condition, vis_2$Type))

p2 = ggplot(vis_2, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Type2, color=Condition, linetype=Type)) + geom_linerange(linewidth = 1, position=position_dodge(width=0.32)) + 
        geom_point(size = 2, position=position_dodge(width=0.32)) + facet_grid(SubModel2 ~ ., scales = "free", space = "free") + geom_hline(aes(yintercept = 0), 
        color = "black", linetype = "dashed") + labs(y = "Regression parameters", x = "") + 
        theme(strip.text.x = element_text(size = 14, face = "bold"), 
          strip.text.y = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14), 
        axis.title.y = element_text(size = 14, face = "bold"), 
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.title.x = element_blank()) + theme(strip.text.y = element_text(angle = 360)) + 
        coord_flip() + theme(panel.spacing = unit(1, 
        "lines")) + scale_color_manual(values=c("C-"="#697f4f", "S++"="black")) + theme(legend.position="bottom")
p2

##############################################################################################
################################################################# Reciprocity terms
# Simple correlations
vis_cm_mask = strand_VPCs(fit_cm_1, n_partitions = 5, merge_noise_terms = TRUE, include_reciprocity = TRUE, mode="cor")
vis_spp_mask = strand_VPCs(fit_spp_1, n_partitions = 5, merge_noise_terms = TRUE, include_reciprocity = TRUE, mode="cor")

df1 = data.frame(vis_cm_mask[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "C-"
df1$Submodel = rep(c("Generalized","Dyadic"),each=1)

df2 = data.frame(vis_spp_mask[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "S++"
df2$Submodel = rep(c("Generalized","Dyadic"),each=1)


df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Context = df$Site 

p3 = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Context, color=Context,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(linewidth = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("C-" = "#697f4f", "S++" = "black")) + theme(legend.position="bottom")

p3


# AMEN-style adjusted correlations
vis_cm_mask = strand_VPCs(fit_cm_1, n_partitions = 4, include_reciprocity = TRUE, mode="adj")
vis_spp_mask = strand_VPCs(fit_spp_1, n_partitions = 4, include_reciprocity = TRUE, mode="adj")

df1 = data.frame(vis_cm_mask[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "C-"
df1$Submodel = rep(c("Generalized","Dyadic"),each=1)

df2 = data.frame(vis_spp_mask[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "S++"
df2$Submodel = rep(c("Generalized","Dyadic"),each=1)


df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Context = df$Site 

p4 = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Context, color=Context,
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
        "lines")) + scale_color_manual(values=c("C-" = "#697f4f", "S++" = "black")) + theme(legend.position="bottom")

p4

##############################################################################################
######################################################################## VPCs
# STRAND 4-way variance partition
vis_cm_mask = strand_VPCs(fit_cm_1, n_partitions = 5, merge_noise_terms = TRUE)
vis_spp_mask = strand_VPCs(fit_spp_1, n_partitions = 5, merge_noise_terms = TRUE)

df1 = data.frame(do.call(rbind, vis_cm_mask[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "C-"
df1$Submodel = rep(c("Aggressed"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic","Error"),1)

df2 = data.frame(do.call(rbind, vis_spp_mask[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "S++"
df2$Submodel = rep(c("Aggressed"),each=4)
df2$Variable2 = rep(c("Focal","Target","Dyadic","Error"),1)

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Aggressed"))

df$Variable2 = factor(df$Variable2)
df$Variable2 = factor(df$Variable2, levels=rev(c("Focal","Target","Dyadic","Error")))

p5 = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Site, color=Site,
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
        "lines")) + scale_color_manual(values=c("C-" = "#697f4f", "S++" = "black"))  + theme(legend.position="bottom")

p5


# AMEN-style 3-way variance partition
vis_cm_mask = strand_VPCs(fit_cm_1, n_partitions = 3)
vis_spp_mask = strand_VPCs(fit_spp_1, n_partitions = 3)

df1 = data.frame(do.call(rbind, vis_cm_mask[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "C-"
df1$Submodel = rep(c("Aggressed"),each=3)
df1$Variable2 = rep(c("Focal","Target","Dyadic+Error"),1)

df2 = data.frame(do.call(rbind, vis_spp_mask[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "S++"
df2$Submodel = rep(c("Aggressed"),each=3)
df2$Variable2 = rep(c("Focal","Target","Dyadic+Error"),1)

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Aggressed"))

df$Variable2 = factor(df$Variable2)
df$Variable2 = factor(df$Variable2, levels=rev(c("Focal","Target","Dyadic+Error")))

p6 = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Site, color=Site,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(linewidth=1, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("C-" = "#697f4f", "S++" = "black"))  + theme(legend.position="bottom")

p6

##############################################################################################
# Masking approach and interaction approach work similarly

##############################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
# Check all the relevant parameters
fit_spp_1$fit$summary()
fit_spp_1$fit$summary("focal_effects")
fit_spp_1$fit$summary("target_effects")


################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit_spp_1$fit$draws(), pars = c("focal_effects[1]","target_effects[1]","sr_L[2,1]","dr_L[2,1]"))


