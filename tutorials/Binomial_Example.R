########################################
#
#   Binomial Analyses  
#
########################################

# Clear working space
rm(list = ls())

# Load libraries
library(STRAND)

# Import data
data(Baboon_Data)

# Number of grooming event and a sample-size measure
nets = list(Grooming = Baboon_Data$Grooming)
exposure_nets = list(Exposure = Baboon_Data$Exposure)

# Dyadic variable: transpose of Presenting
dyad = list(Presenting = t(Baboon_Data$Presenting))


model_dat = make_strand_data(self_report = nets,
                             individual_covariates = NULL, 
                             block_covariates = NULL,
                             dyadic_covariates = dyad,
                             outcome_mode = "binomial",
                             exposure = exposure_nets
                             )

#model
fit =  fit_social_relations_model(data=model_dat,
                              focal_regression = ~ 1,
                              target_regression = ~ 1,
                              dyad_regression = ~  Presenting,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

res = summarize_strand_results(fit)


############################### Plots
vis_1 = strand_caterpillar_plot(res, submodels=c("Focal efffects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, site="HP", only_slopes=TRUE)
vis_1
#ggsave("Baboon_slopes.pdf", vis_1, width=7, height=2.5)

vis_2 = strand_caterpillar_plot(res, submodels=c("Focal efffects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, site="HP", only_technicals=TRUE, only_slopes=FALSE)
vis_2
#ggsave("Baboon_corr.pdf", vis_2, width=6, height=2.5)


