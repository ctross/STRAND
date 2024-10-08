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
data(Colombia_Data)

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

dat = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad)


#model
fit = fit_block_plus_social_relations_model(data=dat,
                                            block_regression = ~ Sex + Ethnicity,
                                            focal_regression = ~ Age + BMI,
                                            target_regression = ~ Age + BMI,
                                            dyad_regression = ~ Distance + Relatedness,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 300, iter_sampling = 300,
                                                                          max_treedepth = NULL, adapt_delta = .98)
)

res = summarize_strand_results(fit)

vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE)
vis_1
#ggsave("Colombia_slopes.pdf", vis_1, width=10, height=5.5)

vis_2 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE,  only_technicals=TRUE, only_slopes=FALSE)
vis_2
#ggsave("Colombia_corr.pdf", vis_2, width=6, height=2.5)

