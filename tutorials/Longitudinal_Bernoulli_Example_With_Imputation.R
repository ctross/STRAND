######################################################################################
#
#   Longitudinal Bernoulli analyses with imputation
#
######################################################################################
library(STRAND)
library(PlvsVltra)

###################################################### Clean data
data(AMENDS_Data)

d = AMENDS_Data$Individual
f = AMENDS_Data$Friendship
s = AMENDS_Data$Sharing
m = AMENDS_Data$Mask
R = AMENDS_Data$Relatedness

 node_labels = rownames(AMENDS_Data$Sharing[[1]])

 ind_1 = data.frame(Age = standardize_strand(d$Age_Wave_1), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize_strand(d$Wealth_Wave_1))
 block_1 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_1 = list(Friend=f$Wave_1, Relatedness = standardize_strand(R))
 mask_1 = list(Share=m$Wave_1)
 out_1 = list(Share=s$Wave_1)
 rownames(ind_1) = rownames(block_1) = node_labels

 ind_2 = data.frame(Age = standardize_strand(d$Age_Wave_2), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize_strand(d$Wealth_Wave_2))
 block_2 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_2 = list(Friend=f$Wave_2, Relatedness = standardize_strand(R))
 mask_2 = list(Share=m$Wave_2)
 out_2 = list(Share=s$Wave_2)
 rownames(ind_2) = rownames(block_2) = node_labels

 ind_3 = data.frame(Age = standardize_strand(d$Age_Wave_3), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize_strand(d$Wealth_Wave_3))
 block_3 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_3 = list(Friend=f$Wave_3, Relatedness = standardize_strand(R))
 mask_3 = list(Share=m$Wave_3)
 out_3 = list(Share=s$Wave_3)
 rownames(ind_3) = rownames(block_3) = node_labels

 ind_4 = data.frame(Age = standardize_strand(d$Age_Wave_4), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize_strand(d$Wealth_Wave_4))
 block_4 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_4 = list(Friend=f$Wave_4, Relatedness = standardize_strand(R))
 mask_4 = list(Share=m$Wave_4)
 out_4 = list(Share=s$Wave_4)
 rownames(ind_4) = rownames(block_4) = node_labels

 dat_long = NULL

  dat_long[[1]] = make_strand_data(
  outcome = out_1,
  mask = mask_1,
  block_covariates = block_1, 
  individual_covariates = ind_1, 
  dyadic_covariates = dyad_1,
  longitudinal = TRUE,
  outcome_mode="bernoulli",
  link_mode = "logit",
  imputation = TRUE
  )

  dat_long[[2]] = make_strand_data(
  outcome = out_2,
  mask = mask_2,
  block_covariates = block_2, 
  individual_covariates = ind_2, 
  dyadic_covariates = dyad_2,
  longitudinal = TRUE,
  outcome_mode="bernoulli",
  link_mode = "logit",
  imputation = TRUE
  )

  dat_long[[3]] = make_strand_data(
  outcome = out_3,
  mask = mask_3,
  block_covariates = block_3, 
  individual_covariates = ind_3, 
  dyadic_covariates = dyad_3,
  longitudinal = TRUE,
  outcome_mode="bernoulli",
  link_mode = "logit",
  imputation = TRUE
  )

  dat_long[[4]] = make_strand_data(
  outcome = out_4,
  mask = mask_4,
  block_covariates = block_4, 
  individual_covariates = ind_4, 
  dyadic_covariates = dyad_4,
  longitudinal = TRUE,
  outcome_mode="bernoulli",
  link_mode = "logit",
  imputation = TRUE
  )

names(dat_long) = paste("Wave", c(1:4))


###################################################### Make data messy
data(AMENDS_Data)

d = AMENDS_Data$Individual
f = AMENDS_Data$Friendship
s = AMENDS_Data$Sharing
m = AMENDS_Data$Mask
R = AMENDS_Data$Relatedness

 node_labels = rownames(AMENDS_Data$Sharing[[1]])

 d$Age_Wave_1[c(11, 15, 27)] = NA
 d$Age_Wave_2[c(11, 15, 27, 34)] = NA
 d$Age_Wave_3[c(11, 15, 27, 39, 56)] = NA
 d$Age_Wave_4[c(11, 15, 27, 60:67)] = NA
 d$Sex[c(13, 15, 37)] = NA

 f$Wave_1[c(1,5,17,20:26,40:45), c(1,5,17,20:26,40:45)] = NA
 f$Wave_2[c(1,5,17,20:26), c(1,5,17,20:26)] = NA
 f$Wave_3[c(1,5,17,40:45), c(1,5,17,40:45)] = NA
 f$Wave_4[c(20:26,40:45), c(20:26,40:45)] = NA

 s$Wave_1[c(11,15,27,30:36,50:55), c(11,15,27,30:36,50:55)] = NA
 s$Wave_2[c(1,5,17,20:26), c(1,5,17,20:26)] = NA
 s$Wave_3[c(1,5,17,40:45), c(1,5,17,40:45)] = NA
 s$Wave_4[c(11,15,27,30:36,50:55), c(11,15,27,30:36,50:55)] = NA

 ind_1 = data.frame(Age = standardize_strand(d$Age_Wave_1), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize_strand(d$Wealth_Wave_1))
 block_1 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_1 = list(Friend=f$Wave_1, Relatedness = standardize_strand(R))
 mask_1 = list(Share=m$Wave_1)
 out_1 = list(Share=s$Wave_1)
 rownames(ind_1) = rownames(block_1) = node_labels

 ind_2 = data.frame(Age = standardize_strand(d$Age_Wave_2), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize_strand(d$Wealth_Wave_2))
 block_2 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_2 = list(Friend=f$Wave_2, Relatedness = standardize_strand(R))
 mask_2 = list(Share=m$Wave_2)
 out_2 = list(Share=s$Wave_2)
 rownames(ind_2) = rownames(block_2) = node_labels

 ind_3 = data.frame(Age = standardize_strand(d$Age_Wave_3), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize_strand(d$Wealth_Wave_3))
 block_3 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_3 = list(Friend=f$Wave_3, Relatedness = standardize_strand(R))
 mask_3 = list(Share=m$Wave_3)
 out_3 = list(Share=s$Wave_3)
 rownames(ind_3) = rownames(block_3) = node_labels

 ind_4 = data.frame(Age = standardize_strand(d$Age_Wave_4), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize_strand(d$Wealth_Wave_4))
 block_4 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_4 = list(Friend=f$Wave_4, Relatedness = standardize_strand(R))
 mask_4 = list(Share=m$Wave_4)
 out_4 = list(Share=s$Wave_4)
 rownames(ind_4) = rownames(block_4) = node_labels

 dat_long_messy = NULL

  dat_long_messy[[1]] = make_strand_data(
  outcome = out_1,
  mask = mask_1,
  block_covariates = block_1, 
  individual_covariates = ind_1, 
  dyadic_covariates = dyad_1,
  longitudinal = TRUE,
  outcome_mode="bernoulli",
  link_mode = "logit",
  imputation = TRUE
  )

  dat_long_messy[[2]] = make_strand_data(
  outcome = out_2,
  mask = mask_2,
  block_covariates = block_2, 
  individual_covariates = ind_2, 
  dyadic_covariates = dyad_2,
  longitudinal = TRUE,
  outcome_mode="bernoulli",
  link_mode = "logit",
  imputation = TRUE
  )

  dat_long_messy[[3]] = make_strand_data(
  outcome = out_3,
  mask = mask_3,
  block_covariates = block_3, 
  individual_covariates = ind_3, 
  dyadic_covariates = dyad_3,
  longitudinal = TRUE,
  outcome_mode="bernoulli",
  link_mode = "logit",
  imputation = TRUE
  )

  dat_long_messy[[4]] = make_strand_data(
  outcome = out_4,
  mask = mask_4,
  block_covariates = block_4, 
  individual_covariates = ind_4, 
  dyadic_covariates = dyad_4,
  longitudinal = TRUE,
  outcome_mode="bernoulli",
  link_mode = "logit",
  imputation = TRUE
  )

names(dat_long_messy) = paste("Wave", c(1:4))

############################################################################## Fit models
fit_1a = fit_longitudinal_model_missings(long_data=dat_long_messy,
                                  block_regression = ~ Ethnicity + Sex,
                                  focal_regression = ~ Age + Wealth,
                                  target_regression = ~ Age + Wealth,
                                  dyad_regression = ~ Friend + Relatedness,
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  bandage_penalty = -1,
                                  mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500, init=0,
                                                              iter_sampling = 500, max_treedepth = 12, adapt_delta = 0.95),
                                  priors=NULL
                                  )

fit_1b = fit_longitudinal_model(long_data=dat_long,
                                  block_regression = ~ Ethnicity + Sex,
                                  focal_regression = ~ Age + Wealth,
                                  target_regression = ~ Age + Wealth,
                                  dyad_regression = ~ Friend + Relatedness,
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  bandage_penalty = -1,
                                  mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500, init=0,
                                                              iter_sampling = 500, max_treedepth = 12, adapt_delta = 0.95),
                                  priors=NULL
                                  )

fit_1c = fit_longitudinal_model_missings(long_data=dat_long,
                                  block_regression = ~ Ethnicity + Sex,
                                  focal_regression = ~ Age + Wealth,
                                  target_regression = ~ Age + Wealth,
                                  dyad_regression = ~ Friend + Relatedness,
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  bandage_penalty = -1,
                                  mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500, init=0,
                                                              iter_sampling = 500, max_treedepth = 12, adapt_delta = 0.95),
                                  priors=NULL
                                  )


fit_2a = fit_longitudinal_model_missings(long_data=dat_long_messy,
                                  block_regression = ~ Ethnicity + Sex,
                                  focal_regression = ~ Age + Wealth,
                                  target_regression = ~ Age + Wealth,
                                  dyad_regression = ~ Friend + Relatedness,
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  bandage_penalty = 0.1,
                                  mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500, init=0,
                                                              iter_sampling = 500, max_treedepth = 12, adapt_delta = 0.95),
                                  priors=NULL
                                  )

fit_2b = fit_longitudinal_model(long_data=dat_long,
                                  block_regression = ~ Ethnicity + Sex,
                                  focal_regression = ~ Age + Wealth,
                                  target_regression = ~ Age + Wealth,
                                  dyad_regression = ~ Friend + Relatedness,
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  bandage_penalty = 0.1,
                                  mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500, init=0,
                                                              iter_sampling = 500, max_treedepth = 12, adapt_delta = 0.95),
                                  priors=NULL
                                  )

fit_2c = fit_longitudinal_model_missings(long_data=dat_long,
                                  block_regression = ~ Ethnicity + Sex,
                                  focal_regression = ~ Age + Wealth,
                                  target_regression = ~ Age + Wealth,
                                  dyad_regression = ~ Friend + Relatedness,
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  bandage_penalty = 0.1,
                                  mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500, init=0,
                                                              iter_sampling = 500, max_treedepth = 12, adapt_delta = 0.95),
                                  priors=NULL
                                  )


######################################## Estimates
# Summaries
res_1a = summarize_strand_results(fit_1a)
res_1b = summarize_strand_results(fit_1b)
res_1c = summarize_strand_results(fit_1c)

res_2a = summarize_strand_results(fit_2a)
res_2b = summarize_strand_results(fit_2b)
res_2c = summarize_strand_results(fit_2c)

vis_1a = strand_caterpillar_plot(res_1a, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-Imputation-Pinkney")
vis_1b = strand_caterpillar_plot(res_1b, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="Old-No-Imputation-Pinkney")
vis_1c = strand_caterpillar_plot(res_1c, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-No-Imputation-Pinkney")

vis_2a = strand_caterpillar_plot(res_2a, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-Imputation-L2")
vis_2b = strand_caterpillar_plot(res_2b, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="Old-No-Imputation-L2")
vis_2c = strand_caterpillar_plot(res_2c, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-No-Imputation-L2")


df = rbind(vis_1a, vis_1b, vis_1c, vis_2a, vis_2b, vis_2c)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Method = df$Site 

df1 = df[which(df$Submodel!="Other estimates"),]

p1a = ggplot2::ggplot(df1, ggplot2::aes(x = Variable, y = Median, group = Method, color=Method,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(linewidth = 1, position = position_dodge(width = 1)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 1.0)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("New-Imputation-L2" = "darkred", "Old-No-Imputation-L2" = "black", "New-No-Imputation-L2" = "goldenrod", 
                                                "New-Imputation-Pinkney" = "blue", "Old-No-Imputation-Pinkney" = "slateblue", "New-No-Imputation-Pinkney" = "orange")) + 
        theme(legend.position="bottom")

p1a

# Very consistent recovery across models





