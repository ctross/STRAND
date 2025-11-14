#####################################################################################
#
#   Multiplex Bernoulli Analyses: Probit versus Logit
#
#   Some have argued that Multiplex SRM models require probit links. Here we show
#   that the choice of link function is irrelevant to anything other than
#   the scale of estimates and the run-time of the models.  
#
#####################################################################################

# Load libraries
 library(PlvsVltra) # For colors: install_github('ctross/PlvsVltra')
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
 Relatedness = standardize_strand(RICH$Relatedness), 
 Friends = RICH$Friends,
 Marriage = RICH$Marriage
)

# Individual data in data-frame
ind = data.frame(
 Age = standardize_strand(RICH$Individual$Age), 
 FoodInsecure = RICH$Individual$FoodInsecure,
 Wealth = standardize_strand(RICH$Individual$Wealth),
 Depressed = RICH$Individual$Depressed
)

# Individual blocking measures
groups = data.frame(
 Ethnicity = as.factor(RICH$Individual$Ethnicity), 
 Sex = as.factor(RICH$Individual$Sex)
)

# Merge data
dat = make_strand_data(
 outcome = outcome,
 block_covariates = groups, 
 individual_covariates = ind, 
 dyadic_covariates = dyad,
 outcome_mode="bernoulli",
 link_mode = "logit",
 multiplex = TRUE
)

# Merge data probit
dat2 = make_strand_data(
 outcome = outcome,
 block_covariates = groups, 
 individual_covariates = ind, 
 dyadic_covariates = dyad,
 outcome_mode="bernoulli",
 link_mode = "probit",
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
 mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.96)
)

# Model 2
fit_2 = fit_multiplex_model(
 data=dat2,
 block_regression = ~ 1,
 focal_regression = ~ 1,
 target_regression = ~ 1,
 dyad_regression = ~ 1,
 mode="mcmc",
 mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.96)
)

############################################# Results
res_1 = summarize_strand_results(fit_1)
res_2 = summarize_strand_results(fit_2)

colors = plvs_vltra("dust_storm", rev=FALSE, elements=c(2,4))
colors = c(colors[1], "grey90", colors[2])

########################################## Corr plots
multiplex_plot(fit_1, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE, height=6, width=7, palette=colors)
multiplex_plot(fit_2, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE, height=6, width=7, palette=colors)

multiplex_plot(fit_1, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE, height=6, width=7, palette=colors)
multiplex_plot(fit_2, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE, height=6, width=7, palette=colors)


######################################################################## VPCs
VPCs_1 = strand_VPCs(fit_1, n_partitions = 5)
VPCs_2 = strand_VPCs(fit_2, n_partitions = 5)

df1 = data.frame(do.call(rbind, VPCs_1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Logit"
df1$Submodel = rep(c("Give","Take","Reduce"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic signal","Dyadic noise + Error"),3)

df2 = data.frame(do.call(rbind, VPCs_2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Probit"
df2$Submodel = rep(c("Give","Take","Reduce"),each=4)
df2$Variable2 = rep(c("Focal","Target","Dyadic signal","Dyadic noise + Error"),3)

df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Give", "Take", "Reduce"))

df$Variable2 = factor(df$Variable2)

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
        "lines")) + scale_color_manual(values=c("Logit" = colors[1], "Probit" = colors[3])) + theme(legend.position="bottom")

p


######################################################################## Reciprocity
VPCs_1 = strand_VPCs(fit_1, n_partitions = 5, include_reciprocity = TRUE)
VPCs_2 = strand_VPCs(fit_2, n_partitions = 5, include_reciprocity = TRUE)

df1 = data.frame(VPCs_1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Logit"
df1$Submodel = rep(c("Generalized","Dyadic"),each=15)

df2 = data.frame(VPCs_2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Probit"
df2$Submodel = rep(c("Generalized","Dyadic"),each=15)


df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Link = df$Site 

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Link, color=Link,
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
        "lines")) + scale_color_manual(values=c("Logit" = colors[1], "Probit" = colors[3])) + theme(legend.position="bottom")

p


######################################################################## Reciprocity, adjusted
# Also possible to make amen-style dyadic reciprocity estimates by rescaling by the ratio of "dyadic variance" to "dyadic variance plus error variance" with mode="adj"
VPCs_1 = strand_VPCs(fit_1, n_partitions = 4, include_reciprocity = TRUE, mode="adj")
VPCs_2 = strand_VPCs(fit_2, n_partitions = 4, include_reciprocity = TRUE, mode="adj")

df1 = data.frame(VPCs_1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "Logit"
df1$Submodel = rep(c("Generalized","Dyadic"),each=15)

df2 = data.frame(VPCs_2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Probit"
df2$Submodel = rep(c("Generalized","Dyadic"),each=15)


df = rbind(df1, df2)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Link = df$Site 

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Link, color=Link,
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
        "lines")) + scale_color_manual(values=c("Logit" = colors[1], "Probit" = colors[3])) + theme(legend.position="bottom")

p

# Basically identical.
