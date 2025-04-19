########################################
#
#   Multiplex Bernoulli Analyses  
#
########################################

# Clear working space
rm(list = ls())

# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors
 colors = plvs_vltra("dust_storm", rev=FALSE, elements=NULL, show=FALSE)

library(STRAND)
library(stringr)
library(ggplot2)
library(psych)

# Load data
 data(RICH_Data)
 RICH = RICH_Data

# Outcomes stored as a labeled list
outcome = list(
 Give = RICH$Give, 
 Take = RICH$Take, 
 Reduce = RICH$Reduce
)

# Dyadic data as a labeled list
dyad = list(
 Relatedness = standardize(RICH$Relatedness), 
 Friends = RICH$Friends,
 Marriage = RICH$Marriage
)

# Individual data in data-frame
RICH$Individual$Age = standardize(RICH$Individual$Age)
RICH$Individual$Wealth = standardize(RICH$Individual$Wealth)
ind = RICH$Individual

# Individual blocking measures
groups = data.frame(
 Ethnicity = as.factor(RICH$Individual$Ethnicity), 
 Sex = as.factor(RICH$Individual$Sex)
)
rownames(groups) = rownames(RICH$Individual)

# Merge data
dat = make_strand_data(
 outcome = outcome,
 block_covariates = groups, 
 individual_covariates = ind, 
 dyadic_covariates = dyad,
 outcome_mode="bernoulli",
 link_mode="logit",
 multiplex = TRUE
)

# Model 1
fit_1 = fit_multiplex_model(
 data=dat,
 block_regression = ~ 1,
 focal_regression = ~ 1,
 target_regression = ~ 1,
 dyad_regression = ~ 1,
 mode="mcmc",
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 500, 
   iter_sampling = 500, 
   max_treedepth = NULL, 
   adapt_delta = 0.98)
)

# Model 2 - Full model, all controls
fit_2 = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 500, 
   iter_sampling = 500, 
   max_treedepth = NULL, 
   adapt_delta = 0.98)
)

# Model 1 results
res_1 = summarize_strand_results(fit_1)

colors = plvs_vltra("dust_storm", rev=FALSE, elements=c(2,4))
colors = c(colors[1], "grey90", colors[2])

multiplex_plot(fit_1, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = "Colombia_1_dyadic.pdf", height=6, width=7, palette=colors)
multiplex_plot(fit_1, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = "Colombia_1_generalized.pdf", height=6, width=7, palette=colors)

# Model 2 results
res_2 = summarize_strand_results(fit_2)

multiplex_plot(fit_2, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE,  save_plot = "Colombia_2_dyadic.pdf", height=6, width=7, palette=colors)
multiplex_plot(fit_2, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE,  save_plot = "Colombia_2_generalized.pdf", height=6, width=7, palette=colors)

# Example Contrasts, new way
process_block_parameters(fit_2, "Afrocolombian to Afrocolombian", "Afrocolombian to Embera", HPDI=0.9)

# Merged plot
vis_1 = strand_caterpillar_plot(res_1, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)
vis_2 = strand_caterpillar_plot(res_2, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)

vis_1$Site = "Base model"
vis_2$Site = "Full model"

vis_1$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", vis_1$Variable)
vis_1$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", vis_1$Variable)
vis_1$Variable = gsub("dyadic effects coeffs, ", "", vis_1$Variable)
vis_1$Variable = gsub("dyadic effects ", "", vis_1$Variable)
vis_1$Variable = gsub("focal effects ", "", vis_1$Variable)
vis_1$Variable = gsub("target effects ", "", vis_1$Variable)


vis_2$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", vis_2$Variable)
vis_2$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", vis_2$Variable)
vis_2$Variable = gsub("dyadic effects coeffs, ", "", vis_2$Variable)
vis_2$Variable = gsub("dyadic effects ", "", vis_2$Variable)
vis_2$Variable = gsub("focal effects ", "", vis_2$Variable)
vis_2$Variable = gsub("target effects ", "", vis_2$Variable)

df = rbind(vis_1, vis_2)

df$Outcome = ifelse(str_detect(df$Variable, "Take"), "Take",
             ifelse(str_detect(df$Variable, "Give"), "Give",
             ifelse(str_detect(df$Variable, "Reduce"), "Reduce",
                    NA)))

df$Outcome = factor(df$Outcome)
df$Outcome = factor(df$Outcome, levels=c("Give", "Take", "Reduce"))

df$Variable = gsub("Give - ", "", df$Variable)
df$Variable = gsub("Take - ", "", df$Variable)
df$Variable = gsub("Reduce - ", "", df$Variable)

df$Variable = gsub(" - Give", "", df$Variable)
df$Variable = gsub(" - Take", "", df$Variable)
df$Variable = gsub(" - Reduce", "", df$Variable)

df$Variable = gsub("sd", "SD", df$Variable)

df$Variable = gsub("FoodInsecure", "Food Insecure", df$Variable)
df$Variable = gsub("Wealth", "Log Wealth", df$Variable)

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Site, color=Site,
        ymin = LI, ymax = HI)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(Submodel ~ 
        Outcome, scales = "free", space = "free") + ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("Base model" = colors[1], "Full model" = colors[3])) + theme(legend.position="bottom")

p

# ggsave("rich_res.pdf",p, width=9, height=4.5)



VPCs_1 = strand_VPCs(fit_1, n_partitions = 4)
VPCs_2 = strand_VPCs(fit_2, n_partitions = 4)

df1 = data.frame(do.call(rbind, VPCs_1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Base model"
df1$Submodel = rep(c("Give","Take","Reduce"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic","Error"),3)

df2 = data.frame(do.call(rbind, VPCs_2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Full model"
df2$Submodel = rep(c("Give","Take","Reduce"),each=4)
df2$Variable2 = rep(c("Focal","Target","Dyadic","Error"),3)

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Give", "Take", "Reduce"))

df$Variable2 = factor(df$Variable2)
df$Variable2 = factor(df$Variable2, levels=rev(c("Focal","Target","Dyadic","Error")))

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Site, color=Site,
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
        "lines")) + scale_color_manual(values=c("Base model" = colors[1], "Full model" = colors[3])) + theme(legend.position="bottom")

p

# ggsave("rich_res_vpc.pdf",p, width=9, height=3)


# Compare with AMEN probit model
Q = function(x){
  return(c(median(x),HPDI(x)))
}

################################################################################################ Fit AMEN
library(amen)
fit_SRM = ame(RICH$Give, family = "bin")
summary(fit_SRM)

## The variance parameters for nodes are reported as variances and covariances
## Converting to correlation (as in the simulation) is effected
## with some straightforward calculations for the generalized reciprocity correlation
ame_sender_var = fit_SRM$VC[,1]
ame_covariance = fit_SRM$VC[,2]
ame_target_var = fit_SRM$VC[,3]
ame_dyadic_cor = fit_SRM$VC[,4] 


ame_generalized_reciprocity = ame_covariance / (sqrt (ame_sender_var * ame_target_var))

## One can also calculate the Variance Partition Coefficients
ame_dyadic_var = 1 
ame_total_var = ame_sender_var + ame_target_var + ame_dyadic_var + 2*ame_covariance

ame_sender_VPC = ame_sender_var / ame_total_var
ame_target_VPC = ame_target_var / ame_total_var
ame_dyadic_VPC = ame_dyadic_var / ame_total_var

AME_Measures_Give = rbind(Q(ame_generalized_reciprocity), Q(ame_dyadic_cor), Q(ame_sender_VPC), Q(ame_target_VPC), rep(NA,3), rep(NA,3), Q(ame_dyadic_VPC))

colnames(AME_Measures_Give) = c("M","L","H")
rownames(AME_Measures_Give) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC")

fit_SRM = ame(RICH$Take, family = "bin")
summary(fit_SRM)

## The variance parameters for nodes are reported as variances and covariances
## Converting to correlation (as in the simulation) is effected
## with some straightforward calculations for the generalized reciprocity correlation
ame_sender_var = fit_SRM$VC[,1]
ame_covariance = fit_SRM$VC[,2]
ame_target_var = fit_SRM$VC[,3]
ame_dyadic_cor = fit_SRM$VC[,4] 


ame_generalized_reciprocity = ame_covariance / (sqrt (ame_sender_var * ame_target_var))

## One can also calculate the Variance Partition Coefficients
ame_dyadic_var = 1 
ame_total_var = ame_sender_var + ame_target_var + ame_dyadic_var + 2*ame_covariance

ame_sender_VPC = ame_sender_var / ame_total_var
ame_target_VPC = ame_target_var / ame_total_var
ame_dyadic_VPC = ame_dyadic_var / ame_total_var

AME_Measures_Take = rbind(Q(ame_generalized_reciprocity), Q(ame_dyadic_cor), Q(ame_sender_VPC), Q(ame_target_VPC), rep(NA,3), rep(NA,3), Q(ame_dyadic_VPC))

colnames(AME_Measures_Take) = c("M","L","H")
rownames(AME_Measures_Take) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC")

fit_SRM = ame(RICH$Reduce, family = "bin")
summary(fit_SRM)

## The variance parameters for nodes are reported as variances and covariances
## Converting to correlation (as in the simulation) is effected
## with some straightforward calculations for the generalized reciprocity correlation
ame_sender_var = fit_SRM$VC[,1]
ame_covariance = fit_SRM$VC[,2]
ame_target_var = fit_SRM$VC[,3]
ame_dyadic_cor = fit_SRM$VC[,4] 


ame_generalized_reciprocity = ame_covariance / (sqrt (ame_sender_var * ame_target_var))

## One can also calculate the Variance Partition Coefficients
ame_dyadic_var = 1 
ame_total_var = ame_sender_var + ame_target_var + ame_dyadic_var + 2*ame_covariance

ame_sender_VPC = ame_sender_var / ame_total_var
ame_target_VPC = ame_target_var / ame_total_var
ame_dyadic_VPC = ame_dyadic_var / ame_total_var

AME_Measures_Reduce = rbind(Q(ame_generalized_reciprocity), Q(ame_dyadic_cor), Q(ame_sender_VPC), Q(ame_target_VPC), rep(NA,3), rep(NA,3), Q(ame_dyadic_VPC))

colnames(AME_Measures_Reduce) = c("M","L","H")
rownames(AME_Measures_Reduce) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC")

VPCs_1 = strand_VPCs(fit_1, n_partitions = 3)


df1 = data.frame(do.call(rbind, VPCs_1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Base model"
df1$Submodel = rep(c("Give","Take","Reduce"),each=3)
df1$Variable2 = rep(c("Focal","Target","Dyadic+Error"),3)
df1$Software = "STRAND"

df2 = rbind(AME_Measures_Give, AME_Measures_Take, AME_Measures_Reduce)
df2 = df2[which(rownames(df2) %in% c("Focal VPC", "Target VPC", "DyadicError VPC")),]

df2b = df1
df2b[,2:4] = df2
df2b$Software = "amen"

df = rbind(df1, df2b)

df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Give", "Take", "Reduce"))

df$Variable2 = factor(df$Variable2)
df$Variable2 = factor(df$Variable2, levels=rev(c("Focal","Target","Dyadic+Error")))

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Software, color=Software,
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
        "lines")) + scale_color_manual(values=c("STRAND" = colors[1], "amen" = colors[3])) + theme(legend.position="bottom")

p

# ggsave("rich_res_vpc_3.pdf",p, width=9, height=3)
# Basically the same as AMEN layer-wise

###############################################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
# Check all the relevant parameters
fit_1$fit$summary()
fit_1$fit$summary("block_effects")

################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit_1$fit$draws(), pars = c("block_effects[1,1]","block_effects[2,1]","sr_L[2,1]","dr_L[2,1]"))

###############################################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
# Check all the relevant parameters
fit_2$fit$summary()
fit_2$fit$summary("focal_effects")
fit_2$fit$summary("target_effects")
fit_2$fit$summary("block_effects")

################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit_2$fit$draws(), pars = c("focal_effects[1,1]","target_effects[1,2]","sr_L[2,1]","dr_L[2,1]"))

