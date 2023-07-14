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

res_cm = summarize_strand_results(fit_cm)


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

res_spp = summarize_strand_results(fit_spp)

########################################################

vis_cm  = strand_caterpillar_plot(res_cm, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="C-", export_as_table = TRUE)
vis_spp = strand_caterpillar_plot(res_spp, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, only_slopes=TRUE, site="S++", export_as_table = TRUE)

vis_1 = rbind(vis_cm, vis_spp)

vis_1$Variable =
c("Female", "Female", 
"RankDiff", "NoOpportunity", 
"RankDiff:NoOpportunity", "Intercept", 
"Female", "Female", 
"RankDiff", "NoOpportunity", 
"RankDiff:NoOpportunity", "Intercept"
)

vis_1$Condition = vis_1$Site

p = ggplot(vis_1, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Condition, color=Condition)) + geom_linerange(size = 1, position=position_dodge(width=0.32)) + 
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

ggsave("Callithrix_slopes.pdf", p, width=10, height=5.5)

vis_2_cm = strand_caterpillar_plot(res_cm, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE,  only_technicals=TRUE, only_slopes=FALSE, site="C-", export_as_table = TRUE)
vis_2_spp = strand_caterpillar_plot(res_spp, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE,  only_technicals=TRUE, only_slopes=FALSE, site="S++", export_as_table = TRUE)
vis_2 = rbind(vis_2_cm, vis_2_spp)
vis_2$Condition = vis_2$Site

p2 = ggplot(vis_2, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Condition, color=Condition)) + geom_linerange(size = 1, position=position_dodge(width=0.32)) + 
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

ggsave("Callithrix_corr.pdf", p2, width=8, height=3.5)

