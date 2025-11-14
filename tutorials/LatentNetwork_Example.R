###########################################################################################################
#
#   Latent Network Analyses  
#
###########################################################################################################

# Load libraries
library(STRAND)
library(ggplot2)

#Load package data
data(FoodSharing_Data)

# Create the STRAND data object
outcome = list(TransferOut = FoodSharing_Data$TransferOut, 
               TransferIn = FoodSharing_Data$TransferIn)

dyad = list(Relatedness = standardize_strand(FoodSharing_Data$Relatedness), 
            Friends = FoodSharing_Data$Friends
            )

FoodSharing_Data$Individual$Age = standardize_strand(FoodSharing_Data$Individual$Age)
FoodSharing_Data$Individual$Wealth = standardize_strand(FoodSharing_Data$Individual$Wealth)
FoodSharing_Data$Individual$GripStrength = standardize_strand(FoodSharing_Data$Individual$GripStrength)
FoodSharing_Data$Individual$Education = standardize_strand(FoodSharing_Data$Individual$Education)
indiv = FoodSharing_Data$Individual 

groups = data.frame(Ethnicity = as.factor(FoodSharing_Data$Individual$Ethnicity), 
                    Sex = as.factor(FoodSharing_Data$Individual$Sex)
                    )

rownames(groups) = rownames(indiv)  

dat = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       outcome_mode = "bernoulli",
                       link_mode = "logit")

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
                                mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500,
                                iter_sampling = 500, max_treedepth = 11, adapt_delta = 0.95)
                                              )

# Summary
res = summarize_strand_results(fit)

# Simple visualizations
vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), normalized=TRUE)
vis_1

vis_2 = strand_caterpillar_plot(res, submodels=c("False positive rate", "Recall of true ties","Theta: question-order effects"), normalized=TRUE)
vis_2

vis_3 = strand_caterpillar_plot(res, submodels=c("False positive rate", "Recall of true ties","Theta: question-order effects"), only_slopes=FALSE, normalized=FALSE, only_technicals=TRUE)
vis_3

            
############################################################## To compute contrasts with new tools, do this:
process_block_parameters(fit, "AFROCOLOMBIAN to AFROCOLOMBIAN", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)
process_block_parameters(fit, "EMBERA to EMBERA", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)
process_block_parameters(fit, "EMBERA to AFROCOLOMBIAN", "AFROCOLOMBIAN to EMBERA", HPDI=0.9)

################################################################# Reciprocity terms
# Simple correlations
cor1 = strand_VPCs(fit, n_partitions = 5, include_reciprocity = TRUE, mode="cor")
cor2 = strand_VPCs(fit, n_partitions = 5, include_reciprocity = TRUE, mode="adj")

df1 = data.frame(cor1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Unadjusted"
df1$Submodel = rep(c("Generalized","Dyadic"),each=1)

df2 = data.frame(cor2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Adjusted"
df2$Submodel = rep(c("Generalized","Dyadic"),each=1)

# Adjusted correlations apply only for dyadic reciprocity. Adjusted correlations give the correlation in random effects+error. Unadjusted give correlation in random effects.

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Method = df$Site 

p3 = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Method, color=Method,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("Unadjusted" = "darkred", "Adjusted" = "black")) + theme(legend.position="bottom")

p3

######################################################################## VPCs
# STRAND variance partition
vpc1 = strand_VPCs(fit, n_partitions = 5)
vpc2 = strand_VPCs(fit, n_partitions = 3)

df1 = data.frame(do.call(rbind, vpc1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "4-way VPC"
df1$Submodel = rep(c("Transfers"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic signal","Dyadic noise + Error"),1)

df2 = data.frame(do.call(rbind, vpc2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "3-way VPC"
df2$Submodel = rep(c("Transfers"),each=3)
df2$Variable2 = rep(c("Focal","Target","Dyadic+Error"),1)

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Transfers"))

df$Variable2 = factor(df$Variable2)

df$Type = df$Site 

p4 = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Type, color=Type,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("4-way VPC" = "darkred", "3-way VPC" = "black"))  + theme(legend.position="bottom")

p4

###############################################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
# Check all the relevant parameters
fit$fit$summary()
fit$fit$summary("focal_effects")
fit$fit$summary("target_effects")
fit$fit$summary("block_effects")

################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit$fit$draws(), pars = c("focal_effects[1]","target_effects[1]","sr_L[2,1]","dr_L[2,1]"))

