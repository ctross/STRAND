#######################################
#
#   Binomial Analyses  
#
########################################

# Clear working space
rm(list = ls())

# Load libraries
library(STRAND)
library(ggplot2)

# Import data
data(Baboon_Data)

# Number of grooming event and a sample-size measure
# Here, the term "exposure" relates to the number of trials for a binomial distribution
nets = list(Grooming = Baboon_Data$Grooming)
exposure = list(Grooming = Baboon_Data$Exposure)

# Dyadic variable: transpose of Presenting
dyad = list(Presenting = t(Baboon_Data$Presenting),
            Threatening = t(Baboon_Data$Threatening))

indiv = Baboon_Data$Individual

block = data.frame(Sex = as.factor(indiv$Sex))
rownames(block) = rownames(indiv)

model_dat = make_strand_data(outcome = nets,
                             individual_covariates = indiv, 
                             block_covariates = block,
                             dyadic_covariates = dyad,
                             outcome_mode = "binomial",
                             exposure = exposure
                             )

# model
fit =  fit_block_plus_social_relations_model(data=model_dat,
                              block_regression = ~ Sex,
                              focal_regression = ~ Age,
                              target_regression = ~ Age,
                              dyad_regression = ~  Presenting + Threatening,
                              mode="mcmc",
                              stan_mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 500, iter_sampling = 500,
                                                          max_treedepth = NULL, adapt_delta = .98)
)

res = summarize_strand_results(fit)


############################### Plots
vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, site="HP", only_slopes=TRUE)
vis_1
#ggsave("Baboon_slopes.pdf", vis_1, width=7, height=2.5)

vis_2 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, site="HP", only_technicals=TRUE, only_slopes=FALSE)
vis_2
#ggsave("Baboon_corr.pdf", vis_2, width=6, height=2.5)

############################################################### To compute contrasts with new tools, do this:
process_block_parameters(input=fit, focal="female to male", base="male to female", HPDI=0.9)

############################################################### To compute contrasts by hand do this:
############################### Posterior contrasts
## Get contrast for block effects
sex_samps = res$sample$srm_model_samples$block_parameters[[2]]
# The above line of code pulls out the MCMC samples for the within- and between-block tie parameters
# The first slot, [[1]], is for the intercept, and the next slot, [[2]], is for the Sex variable.

# Check the order of the factor levels
attr(model_dat,"group_ids_levels")$Sex
# "female" is in slot 1, and "male" is in slot 2

# Now, compute the contrast of female-to-male versus male-to-female ties
mean(sex_samps[,1,2] - sex_samps[,2,1])
HPDI(sex_samps[,1,2] - sex_samps[,2,1], prob=0.89)

############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
fit$fit$summary()
fit$fit$summary("focal_effects")
fit$fit$summary("target_effects")
fit$fit$summary("block_effects")

################# Traceplots
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit$fit$draws(), pars = c("focal_effects[1]","target_effects[1]","sr_L[2,1]","dr_L[2,1]"))

