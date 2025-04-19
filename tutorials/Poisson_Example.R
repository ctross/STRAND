########################################
#
#   Poisson Analyses  
#
# Here we show how to deal with censored data. 
# i.e., cases where no ties were observable.
# Using a mask layer is the prefered method.
#
########################################

# Clear working space
rm(list = ls()) 

# Load libraries
library(STRAND)
library(ggplot2)

# Load data
data(Bat_Data)

# Number of minutes of blood lickings
nets = list(Lick = round(Bat_Data$Lick/60,0))

# Dyadic variables
dyad = list(Relatedness = standardize(Bat_Data$Relatedness), 
            NoOpportunity = Bat_Data$NoOpportunity
              )

# Block variables
group_ids = data.frame(Sex = as.factor(Bat_Data$Individual$Sex))
rownames(group_ids) = rownames(Bat_Data$Lick)

############################################################################### Old style, deal with censoring via interaction at dyad level
model_dat_0 = make_strand_data(outcome = nets,
                              block_covariates = group_ids, 
                              individual_covariates = NULL, 
                              dyadic_covariates = dyad,
                              outcome_mode = "poisson",
                              link_mode = "log"
                              )


# Model 0
fit0 = fit_block_plus_social_relations_model(data=model_dat_0,
                                            block_regression = ~ Sex ,
                                            focal_regression = ~ 1,
                                            target_regression = ~ 1,
                                            dyad_regression = ~ NoOpportunity * Relatedness,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                                        max_treedepth = NULL, adapt_delta = 0.98)
)

# Summary 0
res0 = summarize_strand_results(fit0)

# Simple forest plots
vis_1 = strand_caterpillar_plot(res0, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE,  only_slopes=TRUE)
vis_1
#ggsave("Bat_slopes.pdf", vis_1, width=7, height=2.5)

vis_2 = strand_caterpillar_plot(res0, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, site="XX", only_technicals=TRUE, only_slopes=FALSE)
vis_2
#ggsave("Bat_corr.pdf", vis_2, width=6, height=2.5)

############################################################### To compute contrasts with new tools, do this:
process_block_parameters(input=fit0, focal="Female to Female", base="Male to Female", HPDI=0.9)
process_block_parameters(input=fit0, focal="Female to Male", base="Male to Female", HPDI=0.9)
process_block_parameters(input=fit0, focal="Male to Male", base="Male to Female", HPDI=0.9)

############################################################################### New style, deal with censoring via masking layer
mask = list(Lick = Bat_Data$NoOpportunity) # 0s are possible ties, 1s are censored/impossible. The mask will remove such terms from the likelihood of the model.

model_dat = make_strand_data(outcome = nets,
                              block_covariates = group_ids, 
                              individual_covariates = NULL, 
                              dyadic_covariates = dyad,
                              outcome_mode = "poisson",
                              link_mode = "log",
                              mask = mask
                              )

# Model 1
fit = fit_block_plus_social_relations_model(data=model_dat,
                                            block_regression = ~ Sex ,
                                            focal_regression = ~ 1,
                                            target_regression = ~ 1,
                                            dyad_regression = ~ Relatedness,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                                        max_treedepth = NULL, adapt_delta = 0.98)
)

# Summary 1
res = summarize_strand_results(fit)

# Simple forest plots
vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE,  only_slopes=TRUE)
vis_1
#ggsave("Bat_slopes.pdf", vis_1, width=7, height=2.5)

vis_2 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, site="XX", only_technicals=TRUE, only_slopes=FALSE)
vis_2
#ggsave("Bat_corr.pdf", vis_2, width=6, height=2.5)

# Same results

############################################################### To compute contrasts with new tools, do this:
process_block_parameters(input=fit, focal="Female to Female", base="Male to Female", HPDI=0.9)
process_block_parameters(input=fit, focal="Female to Male", base="Male to Female", HPDI=0.9)
process_block_parameters(input=fit, focal="Male to Male", base="Male to Female", HPDI=0.9)

###############################################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
# Check all the relevant parameters
fit$fit$summary()
fit$fit$summary("block_effects")

################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit$fit$draws(), pars = c("block_effects[1]","block_effects[2]","sr_L[2,1]","dr_L[2,1]"))

