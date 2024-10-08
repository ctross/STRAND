########################################
#
#   Latent Network Analyses  
#
########################################

# Clear working space
rm(list = ls())

# Load libraries
library(STRAND)
library(ggplot2)

#Load package data
data(FoodSharing_Data)

# Create the STRAND data object
outcome = list(TransferOut = FoodSharing_Data$TransferOut, 
               TransferIn = FoodSharing_Data$TransferIn)

dyad = list(Relatedness = FoodSharing_Data$Relatedness, 
            Friends = FoodSharing_Data$Friends
            )

indiv = FoodSharing_Data$Individual 

groups = data.frame(Ethnicity = as.factor(FoodSharing_Data$Individual$Ethnicity), 
                    Sex = as.factor(FoodSharing_Data$Individual$Sex)
                    )

rownames(groups) = rownames(indiv)  

dat = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad)

# Run model
fit = fit_latent_network_model(data=dat,
                                block_regression = ~ Ethnicity + Sex,
                                focal_regression = ~ Age + Wealth + FoodInsecure + CantWork + GripStrength + Depressed,
                                target_regression = ~ Age + Wealth + FoodInsecure + CantWork + GripStrength + Depressed,
                                dyad_regression = ~ Relatedness + Friends,
                                fpr_regression = ~ Age + Wealth + Depressed,
                                rtt_regression = ~ Age + Wealth + Depressed,
                                theta_regression = ~ 1,
                                mode="mcmc",
                                return_predicted_network = FALSE,
                                stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500,
                                iter_sampling = 500, max_treedepth = NULL, adapt_delta = NULL)
                                              )

res = summarize_strand_results(fit)

vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), normalized=TRUE)
vis_1
 # ggsave("Colombia_slopes_latent.pdf", vis_1, width=8, height=8)

vis_2 = strand_caterpillar_plot(res, submodels=c("False positive rate", "Recall of true ties","Theta: question-order effects"), normalized=TRUE)
vis_2
 # ggsave("Colombia_slopes_measurement.pdf", vis_2, width=8, height=8)

vis_3 = strand_caterpillar_plot(res, submodels=c("False positive rate", "Recall of true ties","Theta: question-order effects"), only_slopes=FALSE, normalized=FALSE, only_technicals=TRUE)
vis_3
 # ggsave("Colombia_intercepts_measurement.pdf", vis_3, width=8, height=4)

            
