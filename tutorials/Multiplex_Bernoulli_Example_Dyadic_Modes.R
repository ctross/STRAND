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

# Structure for sample-size sweep
dat_list = NULL
res_list_chol = NULL

Q = 15

# Outcomes stored as a labeled list
outcome = list(
 Give = RICH$Give[1:Q,1:Q], 
 Take = RICH$Take[1:Q,1:Q], 
 Reduce = RICH$Reduce[1:Q,1:Q]
)

# Dyadic data as a labeled list
dyad = list(
 Relatedness = RICH$Relatedness[1:Q,1:Q], 
 Friends = RICH$Friends[1:Q,1:Q],
 Marriage = RICH$Marriage[1:Q,1:Q]
)

# Individual data in data-frame
ind = RICH$Individual[1:Q,]

# Individual blocking measures
groups = data.frame(
 Ethnicity = as.factor(ind$Ethnicity), 
 Sex = as.factor(ind$Sex)
)
rownames(groups) = rownames(ind)

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

############################################### Fast mode
# Model 1 - Full model, all controls
fit_Fast = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = -1,
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.98,
   init=0)
)

res_temp = as.data.frame(fit_Fast$fit$summary(variables ="D_corr"))
res_temp$fit_time = fit_Fast$fit$time()$total
res_temp$node_count = Q
res_temp$method = "cholesky"

############################################### Lasso mode
# Model 2 - Full model, all controls
fit_Lasso = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = 0.02,
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.98,
   init=0)
)

res_Lasso = summarize_strand_results(fit_Lasso)   # 26494 seconds, 24661 seconds

####################################################


# Merged plot
vis_1 = strand_caterpillar_plot(res_Lasso, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)
vis_2 = strand_caterpillar_plot(res_Fast, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)

vis_1$Site = "Base model"
vis_2$Site = "Fast model"

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
df$Variable = gsub("LogWealth", "Log Wealth", df$Variable)

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
        "lines")) + scale_color_manual(values=c("Base model" = colors[1], "Fast model" = colors[3])) + theme(legend.position="bottom")

p

# ggsave("rich_res.pdf",p, width=9, height=4.5)



VPCs_1 = strand_VPCs(fit_Lasso, n_partitions = 4)
VPCs_2 = strand_VPCs(fit_Fast, n_partitions = 4)

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


################################################################### Correlation matrix plots
colors = plvs_vltra("dust_storm", rev=FALSE, elements=c(2,4))
colors = c(colors[1], "grey90", colors[2])

multiplex_plot(fit_Lasso, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE,  height=6, width=7, palette=colors)
multiplex_plot(fit_Fast, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE,  height=6, width=7, palette=colors)

multiplex_plot(fit_Lasso, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE, height=6, width=7, palette=colors)
multiplex_plot(fit_Fast, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE, height=6, width=7, palette=colors)





######################################################################## Reciprocity
VPCs_1 = strand_VPCs(fit_Lasso, n_partitions = 4, include_reciprocity = TRUE)
VPCs_2 = strand_VPCs(fit_Fast, n_partitions = 4, include_reciprocity = TRUE)

df1 = data.frame(VPCs_1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "L2_Norm"
df1$Submodel = rep(c("Generalized","Dyadic"),each=15)

df2 = data.frame(VPCs_2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Pinkney"
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
        "lines")) + scale_color_manual(values=c("L2_Norm" = colors[1], "Pinkney" = colors[3])) + theme(legend.position="bottom")

p
# ggsave("rich_res_recip_probit_v_logit.pdf",p, width=9, height=9)

######################################################################## Reciprocity, adjusted
# Also possible to make amen-style dyadic reciprocity estimates by scaling by ratio of "dyadic variance" to "dyadic variance plus error variance"
VPCs_1 = strand_VPCs(fit_Lasso, n_partitions = 4, include_reciprocity = TRUE, mode="adj")
VPCs_2 = strand_VPCs(fit_Fast, n_partitions = 4, include_reciprocity = TRUE, mode="adj")

df1 = data.frame(VPCs_1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "L2_Norm"
df1$Submodel = rep(c("Generalized","Dyadic"),each=15)

df2 = data.frame(VPCs_2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "Pinkney"
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
        "lines")) + scale_color_manual(values=c("L2_Norm" = colors[1], "Pinkney" = colors[3])) + theme(legend.position="bottom")

p
