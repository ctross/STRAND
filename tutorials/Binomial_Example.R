##############################################################################################
#
#   Binomial Analyses  
#
##############################################################################################

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
            Threatening = standardize_strand(t(Baboon_Data$Threatening),center=FALSE)
            )

Baboon_Data$Individual$Age = standardize_strand(Baboon_Data$Individual$Age)
indiv = Baboon_Data$Individual

block = data.frame(Sex = as.factor(indiv$Sex))
rownames(block) = rownames(indiv)

model_dat = make_strand_data(outcome = nets,
                             individual_covariates = indiv, 
                             block_covariates = block,
                             dyadic_covariates = dyad,
                             outcome_mode = "binomial",
                             link_mode = "logit",
                             exposure = exposure
                             )

# model
fit =  fit_block_plus_social_relations_model(data=model_dat,
                              block_regression = ~ Sex,
                              focal_regression = ~ Age,
                              target_regression = ~ Age,
                              dyad_regression = ~  Presenting + Threatening,
                              mode="mcmc",
                              mcmc_parameters = list(chains = 1, refresh = 1,
                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                          max_treedepth = 11, adapt_delta = 0.95)
)

res = summarize_strand_results(fit)

##############################################################################################
############################### Plots
vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, site="HP", only_slopes=TRUE)
vis_1

vis_2 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, site="HP", only_technicals=TRUE, only_slopes=FALSE)
vis_2


############################################################### To compute contrasts with new tools, do this:
process_block_parameters(input=fit, focal="female to male", base="male to female", HPDI=0.9)
process_block_parameters(input=fit, focal="female to female", base="male to female", HPDI=0.9)
process_block_parameters(input=fit, focal="male to male", base="male to female", HPDI=0.9)


################################################################# Reciprocity terms
# Simple correlations
cor1 = strand_VPCs(fit, n_partitions = 5, merge_noise_terms = TRUE, include_reciprocity = TRUE, mode="cor")
cor2 = strand_VPCs(fit, n_partitions = 5, merge_noise_terms = TRUE, include_reciprocity = TRUE, mode="adj")

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
        "lines")) + scale_color_manual(values=c("Unadjusted" = "darkred", "Adjusted" = "black")) + theme(legend.position="bottom")

p3

######################################################################## VPCs
# STRAND variance partitions
vpc1 = strand_VPCs(fit, n_partitions = 4)
vpc2 = strand_VPCs(fit, n_partitions = 3)
vpc3 = strand_VPCs(fit, n_partitions = 5, merge_noise_terms = TRUE)

df1 = data.frame(do.call(rbind, vpc1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "4-way VPC"
df1$Submodel = rep(c("Grooming"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic","Error"),1)

df2 = data.frame(do.call(rbind, vpc2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "3-way VPC"
df2$Submodel = rep(c("Grooming"),each=3)
df2$Variable2 = rep(c("Focal","Target","Dyadic + Error"),1)

df3 = data.frame(do.call(rbind, vpc3[[2]]))
colnames(df3) = c("Variable", "Median", "L", "H", "Mean", "SD")
df3$Site = "4-way VPC (alternative)"
df3$Submodel = rep(c("Grooming"),each=4)
df3$Variable2 = rep(c("Focal","Target","Dyadic Signal","Dyadic Noise + Error"),1)

df = rbind(df1, df2, df3)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Grooming"))

df$Variable2 = factor(df$Variable2)
df$Variable2 = factor(df$Variable2, levels=rev(c("Focal","Target","Dyadic","Error","Dyadic + Error", "Dyadic Signal","Dyadic Noise + Error")))

df$Type = df$Site 

p4 = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Type, color=Type,
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
        "lines")) + scale_color_manual(values=c("4-way VPC" = "darkred", "3-way VPC" = "black", "4-way VPC (alternative)" = "orange1"))  + theme(legend.position="bottom")

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

