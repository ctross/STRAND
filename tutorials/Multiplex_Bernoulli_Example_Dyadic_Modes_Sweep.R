########################################
#
#   Multiplex Bernoulli Analyses  
#
########################################

# Clear working space
#rm(list = ls())

# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors
 colors = plvs_vltra("dust_storm", rev=FALSE, elements=NULL, show=FALSE)
 colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

library(STRAND)
library(stringr)
library(ggplot2)
library(psych)
library(parallel)
library(purrr)

# Load data
 data(RICH_Data)
 RICH = RICH_Data

# Structure for sample-size sweep
dat_list_cholesky = NULL
dat_list_l2norm = NULL
samp_seq = seq(23, 93, by=5)
ticker = 0

for(s in samp_seq){
 ticker = 1 + ticker
 Q = s

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

dat$node_count = s

dat$model_to_deploy = 1
dat_list_cholesky[[ticker]] = dat

dat$model_to_deploy = 2
dat_list_l2norm[[ticker]] = dat
}

dat_list = c(dat_list_cholesky, dat_list_l2norm)

##################################### Now run the simulations

deploy_on_cores = function(dat){
if(dat$model_to_deploy == 1){
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
res_temp$node_count = dat$node_count
res_temp$method = "cholesky"
}

if(dat$model_to_deploy == 2){
############################################### Lasso mode
# Model 2 - Full model, all controls
fit_Lasso = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = 0.01,
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

res_temp = as.data.frame(fit_Lasso$fit$summary(variables ="D_corr"))
res_temp$fit_time = fit_Lasso$fit$time()$total
res_temp$node_count = dat$node_count
res_temp$method = "l2norm"
}

#################################################### Return
 return(res_temp)
}

############################################################### Fit on server
fit = mclapply(dat_list, function(z){
                       deploy_on_cores(z)
                       }, mc.cores = 30)

res = do.call(rbind, fit)
save(res, file = "res_compare.RData")

################################################################### Process
sum_stats = function(x){return(c(mean(x), HPDI(x, prob=0.9)))}

res_thin = res[which(res$mean != 1),]
res_thin$group = paste0(res_thin$method,"_", res_thin$node_count)
res_thin$method_numeric = ifelse(res_thin$method=="l2norm",1,0)

l = list( 
 aggregate(res_thin$method_numeric~res_thin$group, FUN=median),
 aggregate(res_thin$node_count~res_thin$group, FUN=median),
 aggregate(res_thin$fit_time~res_thin$group, FUN=median),
 aggregate(res_thin$ess_bulk~res_thin$group, FUN=sum_stats),
 aggregate(res_thin$ess_tail~res_thin$group, FUN=sum_stats))

res_merged = purrr::reduce(.x = l, merge, by = c('res_thin$group'), all = TRUE)
res_merged = data.frame(as.matrix(data.frame(res_merged)))

colnames(res_merged) = c("group", "method_numeric", "node_count", "fit_time", "ess_bulk_mean", "ess_bulk_05", "ess_bulk_95", "ess_tail_mean", "ess_tail_05", "ess_tail_95")
res_merged$method = ifelse(res_merged$method_numeric==1,"l2norm","cholesky")

for(i in 2:10){
 res_merged[,i] = as.numeric(res_merged[,i])
}

res_long = NULL
legal_set = c("D_corr[1,2]", "D_corr[1,3]", "D_corr[2,3]",  "D_corr[1,4]", "D_corr[2,5]", "D_corr[3,6]",  "D_corr[1,5]", "D_corr[1,6]", "D_corr[2,6]")

for(i in 1:length(legal_set)){
res_thin = res[which(res$variable == legal_set[i]),]
res_thin$group = paste0(res_thin$method,"_", res_thin$node_count)
res_thin$method_numeric = ifelse(res_thin$method=="l2norm",1,0)
res_long[[i]] = res_thin
}

res_thin = do.call(rbind, res_long)
res_thin$grouping = paste0(res_thin$method, res_thin$variable)

################################################################### Plot
# ggplot(data=res_thin, aes(x=node_count, y=ess_bulk, group=grouping, color=method, fill = method)) + geom_path() 
# ggplot(data=res_thin, aes(x=node_count, y=ess_tail, group=grouping, color=method, fill = method)) + geom_path() 

################################################################### Plot
res_merged$Method = res_merged$method

p1 = ggplot(data=res_merged, aes(x=node_count, y=fit_time, group=Method, color=Method, fill = Method)) +
  geom_line(linewidth=1.5) + ylab("Computation time (seconds)") + xlab("Nodes in network") + 
  scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7])) +
  theme(legend.position="bottom")

ggsave("Compare_time.pdf", p1, height=4, width=4)

p2 = ggplot(data=res_merged, aes(x=node_count, y=ess_bulk_mean, group=Method, color=Method, fill = Method)) +
  geom_line(linewidth=1.5) + ylab("Effective samples (bulk)") + xlab("Nodes in network") + ylim(0,1000) +
  scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7])) +
  theme(legend.position="bottom")

ggsave("Compare_Bulk_ESS.pdf", p2, height=4, width=4)

p3 = ggplot(data=res_merged, aes(x=node_count, y=ess_tail_mean, group=Method, color=Method, fill = Method)) +
  geom_line(linewidth=1.5) + ylab("Effective samples (tail)") + xlab("Nodes in network") + ylim(0,1000) +
  scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7])) +
  theme(legend.position="bottom")

ggsave("Compare_Tail_ESS.pdf", p3, height=4, width=4)

############################################################################################################## Coeff plots

deploy_on_cores2 = function(dat){
if(dat$model_to_deploy == 1){
############################################### Fast mode
# Model 1 - Full model, all controls
fit = fit_multiplex_model(
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

}

if(dat$model_to_deploy == 2){
############################################### Lasso mode
# Model 2 - Full model, all controls
fit = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = 0.01,
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

}

#################################################### Return
 return(fit)
}

############################################################### Fit on server
fit_full = mclapply(dat_list[c(15,30)], function(z){
                       deploy_on_cores2(z)
                       }, mc.cores = 2)

res_1 = summarize_strand_results(fit_full[[2]])
res_2 = summarize_strand_results(fit_full[[1]])

vis_1 = strand_caterpillar_plot(res_1, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)
vis_2 = strand_caterpillar_plot(res_2, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)

vis_1$Site = "l2norm"
vis_2$Site = "cholesky"

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

df$Variable = factor(df$Variable)
df$Variable = factor(df$Variable, levels=c("SD", "Age", "Wealth", "Food Insecure", "Depressed", "Relatedness", "Friends", "Marriage"))


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
        "lines")) + scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7])) + theme(legend.position="bottom")

 ggsave("rich_res.pdf",p, width=9, height=4.5)



VPCs_1 = strand_VPCs(fit_full[[2]], n_partitions = 4)
VPCs_2 = strand_VPCs(fit_full[[1]], n_partitions = 4)

df1 = data.frame(do.call(rbind, VPCs_1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "l2norm"
df1$Submodel = rep(c("Give","Take","Reduce"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic","Error"),3)

df2 = data.frame(do.call(rbind, VPCs_2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "cholesky"
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
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(Submodel ~ ., scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7])) + theme(legend.position="bottom")

 ggsave("rich_vpc.pdf",p, width=4.5, height=9)


######################################################################## Reciprocity
VPCs_1 = strand_VPCs(fit_full[[2]], n_partitions = 4, include_reciprocity = TRUE)
VPCs_2 = strand_VPCs(fit_full[[1]], n_partitions = 4, include_reciprocity = TRUE)

df1 = data.frame(VPCs_1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "l2norm"
df1$Submodel = rep(c("Generalized","Dyadic"),each=15)

df2 = data.frame(VPCs_2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "cholesky"
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
        "lines")) + scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7])) + theme(legend.position="bottom")


 ggsave("rich_recip.pdf",p, width=9, height=9)
