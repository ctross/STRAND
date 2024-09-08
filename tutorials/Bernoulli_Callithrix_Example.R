########################################
#
#   Bernoulli Analyses  
#
########################################

# Clear working space
rm(list = ls())

# Load libraries
library(STRAND)
library(ggplot2)

#Load package data
data(Callithrix_Data)

############################################################################################## 
## Two ways to deal with censoring/masking/dyads with no possiblity of interacting
##############################################################################################

# New Method: Use a masking layer.
# Old Method: Use interaction.

# We show both here. First the masking approach. Then the interaction method. The masking method is better,
# as it yields more accurate estimation of the random effects structure.

# Either way, start by creating a matrix of 0s and 1s for each layer of the outcome list. Use the same names as in the outcome list.
# 0s represent non-censored ties---i.e., dyadic outcomes for which there is a possibility of an interaction/tie/link/flow/etc.
# 1s represent censored ties---i.e., dyadic outcomes for which there is no possibility of an interaction/tie/link/flow/etc.

# In this example, the variable: Callithrix_Data[[1]]$NoOpportunity is a censoring mask. 
# There are 4 groups of Callithrix, and ties can only occur within groups. Individuals are sorted by group id. Thus, running:
Callithrix_Data[[1]]$NoOpportunity
# shows that there are 4 blocks of 0s along the diagonal. These 0s represent "real data" appearing in the outcome matrix.
# There are also 1s appearing in the rest of the matrix. These 1s represent between-group dyad interactions that are impossible (e.g., because the groups live in different areas). 

# By running:
# mask = list(Aggressed = Callithrix_Data[[1]]$NoOpportunity)
# and then adding "mask=mask" to the make_strand_data() function, STRAND will only model outcome[i,j,m] if mask[i,j,m]==0.

# If make_strand_data() is called, and mask != NULL, then the mask layer will be used by any model that reads this data.

############################################################################################## 
## Censoring mask
##############################################################################################
################################################################ Condition, C-, low food
# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[1]]$Aggressed)
mask = list(Aggressed = Callithrix_Data[[1]]$NoOpportunity) # 0s are possible ties, 1s are censored/impossible. The mask will remove such terms from the likelihood of the model.
dyad = list(RankDiff = Callithrix_Data[[1]]$RankDiff)
indiv =  data.frame(Female = Callithrix_Data[[1]]$Female)

model_dat_cm = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             mask=mask,                    # if mask != NULL, then it will be used by any model that reads this data.
                             outcome_mode = "bernoulli"
                             )
# Model
fit_cm =  fit_social_relations_model(data=model_dat_cm,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Get results.
res_cm_mask = summarize_strand_results(fit_cm)


################################################################ Condition, S++, high food
# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[2]]$Aggressed)
mask = list(Aggressed = Callithrix_Data[[2]]$NoOpportunity) # 0s are possible ties, 1s are censored/impossible. The mask will remove such terms from the model.
dyad = list(RankDiff = Callithrix_Data[[2]]$RankDiff)
indiv =  data.frame(Female = Callithrix_Data[[2]]$Female)

model_dat_spp = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             mask=mask,
                             outcome_mode = "bernoulli"
                             )
# Model
fit_spp =  fit_social_relations_model(data=model_dat_spp,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Get results.
res_spp_mask = summarize_strand_results(fit_spp)

############################################################################################## 
## Interaction method
##############################################################################################
################################################################ Condition, C-, low food
# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[1]]$Aggressed)

dyad = list(RankDiff = Callithrix_Data[[1]]$RankDiff,
            NoOpportunity = Callithrix_Data[[1]]$NoOpportunity)

indiv =  data.frame(Female = Callithrix_Data[[1]]$Female)

model_dat_cm = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             outcome_mode = "bernoulli"
                             )
# Model
fit_cm =  fit_social_relations_model(data=model_dat_cm,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff * NoOpportunity,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Get results.
res_cm_interaction = summarize_strand_results(fit_cm)


################################################################ Condition, S++, high food
# Create the STRAND data object
outcome = list(Aggressed = Callithrix_Data[[2]]$Aggressed)

dyad = list(RankDiff = Callithrix_Data[[2]]$RankDiff,
            NoOpportunity = Callithrix_Data[[2]]$NoOpportunity)

indiv =  data.frame(Female = Callithrix_Data[[2]]$Female)

model_dat_spp = make_strand_data(outcome = outcome,
                             individual_covariates = indiv, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             outcome_mode = "bernoulli"
                             )
# Model
fit_spp =  fit_social_relations_model(data=model_dat_spp,
                              focal_regression = ~ Female,
                              target_regression = ~ Female,
                              dyad_regression = ~  RankDiff * NoOpportunity,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 100,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

# Get results.
res_spp_interaction = summarize_strand_results(fit_spp)

#################################################################### Compare effects of covariates
res_cm_mask_tab  = strand_caterpillar_plot(res_cm_mask, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="C-", export_as_table = TRUE)
res_spp_mask_tab = strand_caterpillar_plot(res_spp_mask, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="S++", export_as_table = TRUE)

res_cm_interaction_tab  = strand_caterpillar_plot(res_cm_interaction, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="C-", export_as_table = TRUE)
res_spp_interaction_tab = strand_caterpillar_plot(res_spp_interaction, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="S++", export_as_table = TRUE)

vis_1 = rbind(res_cm_mask_tab, res_spp_mask_tab, res_cm_interaction_tab, res_spp_interaction_tab)

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

p = ggplot(vis_1, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Type2, color=Condition, linetype=Type)) + geom_linerange(size = 1, position=position_dodge(width=0.32)) + 
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

#ggsave("Callithrix_slopes.pdf", p, width=10, height=5.5)

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

p2 = ggplot(vis_2, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Type2, color=Condition, linetype=Type)) + geom_linerange(size = 1, position=position_dodge(width=0.32)) + 
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

#ggsave("Callithrix_corr.pdf", p2, width=8, height=3.5)

