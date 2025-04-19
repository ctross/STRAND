########################################
#
#   Multiplex Binomial Analyses  
#
########################################

# Clear working space
rm(list = ls())

# install_github('ctross/PlvsVltra')
library(PlvsVltra) # For colors

library(STRAND)
library(stringr)
library(ggplot2)
library(psych)

# Load data
data(Baboon_Data)
Baboon = Baboon_Data

# Outcomes stored as a labeled list
outcome = list(
 Groom = Baboon$Grooming, 
 Present = Baboon$Presenting, 
 Threat = Baboon$Threatening
)

# Exposure stored as a labeled list
exposure = list(
 Groom = Baboon$Exposure, 
 Present = Baboon$Exposure, 
 Threat = Baboon$Exposure
)

# Individual data in data-frame
Baboon$Individual$Age = standardize(Baboon$Individual$Age)
ind = Baboon$Individual 

# Merge data
dat = make_strand_data(
 outcome = outcome,
 block_covariates = NULL, 
 individual_covariates = ind, 
 dyadic_covariates = NULL,
 exposure = exposure,
 outcome_mode="binomial",
 link_mode="logit",
 multiplex = TRUE
)

# Model 3
fit3 = fit_multiplex_model(
 data=dat,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ 1,
 mode="mcmc",
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = NULL, 
   adapt_delta = 0.98)
)

res3 = summarize_strand_results(fit3)

#################################################################
# Correlation matrix plots
colors = plvs_vltra("dust_storm", rev=FALSE, elements=c(2,4))
colors = c(colors[1], "grey90", colors[2])
multiplex_plot(fit3, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE,  height=6, width=7, palette=colors)
multiplex_plot(fit3, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE,  height=6, width=7, palette=colors)

# Slopes plot
vis3 = strand_caterpillar_plot(res3, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=TRUE)

vis3$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", vis3$Variable)
vis3$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", vis3$Variable)
vis3$Variable = gsub("dyadic effects coeffs, ", "", vis3$Variable)
vis3$Variable = gsub("dyadic effects ", "", vis3$Variable)
vis3$Variable = gsub("focal effects ", "", vis3$Variable)
vis3$Variable = gsub("target effects ", "", vis3$Variable)
vis3$Variable = gsub("Sexmale", "Male", vis3$Variable)
vis3$Variable = gsub("sd", "SD", vis3$Variable)

df = rbind(vis3)

df$Outcome = ifelse(str_detect(df$Variable, "Groom"), "Groom",
             ifelse(str_detect(df$Variable, "Present"), "Present",
             ifelse(str_detect(df$Variable, "Threat"), "Threat",
                    NA)))

df$Variable = gsub("Groom - ", "", df$Variable)
df$Variable = gsub("Present - ", "", df$Variable)
df$Variable = gsub("Threat - ", "", df$Variable)

df$Variable = gsub(" - Groom", "", df$Variable)
df$Variable = gsub(" - Present", "", df$Variable)
df$Variable = gsub(" - Threat", "", df$Variable)

p = ggplot(df, aes(x = Variable, y = Median, 
        ymin = LI, ymax = HI)) + geom_hline(aes(yintercept = 0), 
        color = "black", linetype = "dashed") + geom_linerange(size = 1, color=colors[3]) + 
        geom_point(size = 2, color=colors[3]) + facet_grid(Submodel ~ 
        Outcome, scales = "free", space="free_y") + labs(y = "Regression parameters", 
        x = "") + theme(strip.text.x = element_text(size = 12, 
        face = "bold"), strip.text.y = element_text(size = 12, 
        face = "bold"), axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 14, face = "bold"), 
        axis.title.x = element_blank()) + theme(strip.text.y = element_text(angle = 360)) + 
        coord_flip() + theme(panel.spacing = grid::unit(1, 
        "lines")) + theme(legend.position="bottom")

p


######################## VPCs
VPCs_1 = strand_VPCs(fit3, n_partitions = 4)

df1 = data.frame(do.call(rbind, VPCs_1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Base model"
df1$Submodel = rep(c("Groom","Present","Threat"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic","Error"),3)


df = rbind(df1)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Groom", "Present", "Threat"))

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


######################################################################## Reciprocity
VPCs_1 = strand_VPCs(fit3, n_partitions = 4, include_reciprocity = TRUE)

df1 = data.frame(VPCs_1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Base Cor"
df1$Submodel = rep(c("Generalized","Dyadic"),each=15)

df1$Median = as.numeric(df1$Median)
df1$L = as.numeric(df1$L)
df1$H = as.numeric(df1$H)

df1$Submodel = factor(df1$Submodel)

df1$Model = df1$Site 

## Adjusted
VPCs_2 = strand_VPCs(fit3, n_partitions = 4, include_reciprocity = TRUE, mode="adj")

df2 = data.frame(VPCs_2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Adjusted Cor"
df2$Submodel = rep(c("Generalized","Dyadic"),each=15)

df2$Median = as.numeric(df2$Median)
df2$L = as.numeric(df2$L)
df2$H = as.numeric(df2$H)

df2$Submodel = factor(df2$Submodel)

df2$Model = df2$Site 

df = rbind(df1,df2)

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Model, color=Model,
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
        "lines")) + scale_color_manual(values=c("Base Cor" = colors[3], "Adjusted Cor" = colors[5])) + theme(legend.position="bottom")

p


###############################################################################################################
############################### Model fit diagnostics
# Note that these functions are run on the raw stan object, so variables are not mapped yet to parameter names.
# See Supplementary Appendix for a list of parameter names
################# Rhat and Effective Samples
# Check all the relevant parameters
fit3$fit$summary()
fit3$fit$summary("focal_effects")
fit3$fit$summary("target_effects")
fit3$fit$summary("block_effects")

################# Traceplots
# Check all the relevant parameters
bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(fit3$fit$draws(), pars = c("focal_effects[1,1]","target_effects[1,2]","sr_L[2,1]","dr_L[2,1]"))

