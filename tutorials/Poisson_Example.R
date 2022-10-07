########################################
#
#   Poisson Analyses  
#
########################################

# Clear working space
rm(list = ls()) 

# Load libraries
library(STRAND)

# Load data
data(Bat_Data)

# Number of minutes of blood lickings
nets = list(Lick = round(Bat_Data$Lick/60,0))

# Dyadic variables
dyad = list(Relatedness = Bat_Data$Relatedness, 
            NoOppertunity = Bat_Data$NoOppertunity
              )

# Block variables
group_ids = data.frame(Sex = as.factor(Bat_Data$Sex))

model_dat = make_strand_data(self_report = nets,
                              block_covariates = group_ids, 
                              individual_covariates = NULL, 
                              dyadic_covariates = dyad,
                              outcome_mode = "poisson"
                              )

# Model
fit = fit_block_plus_social_relations_model(data=model_dat,
                                            block_regression = ~ Sex ,
                                            focal_regression = ~ 1,
                                            target_regression = ~ 1,
                                            dyad_regression = ~ NoOppertunity + Relatedness,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                                        max_treedepth = NULL, adapt_delta = .98)
)

res = summarize_strand_results(fit)

vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE,  only_slopes=TRUE)
vis_1
#ggsave("Bat_slopes.pdf", vis_1, width=7, height=2.5)

vis_2 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, site="XX", only_technicals=TRUE, only_slopes=FALSE)
vis_2
#ggsave("Bat_corr.pdf", vis_2, width=6, height=2.5)

