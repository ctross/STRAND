######################################### Colombia example
library(STRAND)
library(PlvsVltra)

data(AMENDS_Data)

d = AMENDS_Data$Individual
f = AMENDS_Data$Friendship
s = AMENDS_Data$Sharing
m = AMENDS_Data$Mask
R = AMENDS_Data$Relatedness

 node_labels = rownames(AMENDS_Data$Sharing[[1]])

 ind_1 = data.frame(Age = standardize(d$Age_Wave_1), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize(d$Wealth_Wave_1))
 block_1 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_1 = list(Friend=f$Wave_1, Relatedness = R)
 mask_1 = list(Share=m$Wave_1)
 out_1 = list(Share=s$Wave_1)
 rownames(ind_1) = rownames(block_1) = node_labels

 ind_2 = data.frame(Age = standardize(d$Age_Wave_2), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize(d$Wealth_Wave_2))
 block_2 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_2 = list(Friend=f$Wave_2, Relatedness = R)
 mask_2 = list(Share=m$Wave_2)
 out_2 = list(Share=s$Wave_2)
 rownames(ind_2) = rownames(block_2) = node_labels

 ind_3 = data.frame(Age = standardize(d$Age_Wave_3), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize(d$Wealth_Wave_3))
 block_3 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_3 = list(Friend=f$Wave_3, Relatedness = R)
 mask_3 = list(Share=m$Wave_3)
 out_3 = list(Share=s$Wave_3)
 rownames(ind_3) = rownames(block_3) = node_labels

 ind_4 = data.frame(Age = standardize(d$Age_Wave_4), Sex = d$Sex, Ethnicity=d$Ethnicity, Wealth=standardize(d$Wealth_Wave_4))
 block_4 = data.frame(Sex = factor(d$Sex), Ethnicity=factor(d$Ethnicity))
 dyad_4 = list(Friend=f$Wave_4, Relatedness = R)
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
  outcome_mode="bernoulli"
  )

  dat_long[[2]] = make_strand_data(
  outcome = out_2,
  mask = mask_2,
  block_covariates = block_2, 
  individual_covariates = ind_2, 
  dyadic_covariates = dyad_2,
  longitudinal = TRUE,
  outcome_mode="bernoulli"
  )

  dat_long[[3]] = make_strand_data(
  outcome = out_3,
  mask = mask_3,
  block_covariates = block_3, 
  individual_covariates = ind_3, 
  dyadic_covariates = dyad_3,
  longitudinal = TRUE,
  outcome_mode="bernoulli"
  )

  dat_long[[4]] = make_strand_data(
  outcome = out_4,
  mask = mask_4,
  block_covariates = block_4, 
  individual_covariates = ind_4, 
  dyadic_covariates = dyad_4,
  longitudinal = TRUE,
  outcome_mode="bernoulli"
  )

names(dat_long) = paste("Wave", c(1:4))

##################### Model 0
fit_0 = fit_longitudinal_model(long_data=dat_long,
                                  block_regression = ~ 1,
                                  focal_regression = ~ 1,
                                  target_regression = ~ 1,
                                  dyad_regression = ~ 1,
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500,
                                                              iter_sampling = 500, max_treedepth = 12, adapt_delta = NULL),
                                  priors=NULL
                                  )

res_0 = summarize_longitudinal_bsrm_results(fit_0)

palLH = plvs_vltra("mystic_mausoleum", elements=c(7,5))
palM = plvs_vltra("skinny_dipping", elements=c(4))
pal = c(palLH[1], "grey80", palLH[2])

multiplex_plot(fit_0, type="dyadic", HPDI=0.9, plot = TRUE, height=7, width=9, palette = pal, save_plot="Colombia_dyadic_0.pdf")
multiplex_plot(fit_0, type="generalized", HPDI=0.9, plot = TRUE, height=7, width=9, palette = pal, save_plot="Colombia_generalized_0.pdf")

######################## Model 1
fit_1 = fit_longitudinal_model(long_data=dat_long,
                                  block_regression = ~ Ethnicity,
                                  focal_regression = ~ Age + Sex + Wealth,
                                  target_regression = ~ Age + Sex + Wealth,
                                  dyad_regression = ~ Friend + Relatedness,
                                  coefficient_mode="varying",
                                  random_effects_mode="fixed",
                                  mode="mcmc",
                                  stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500,
                                                              iter_sampling = 500, max_treedepth = 12, adapt_delta = NULL),
                                  priors=NULL
                                  )

res_1 = summarize_longitudinal_bsrm_results(fit_1)

palLH = plvs_vltra("mystic_mausoleum", elements=c(7,5))
palM = plvs_vltra("skinny_dipping", elements=c(4))
pal = c(palLH[1], "grey80", palLH[2])

multiplex_plot(fit_1, type="dyadic", HPDI=0.9, plot = TRUE, height=7, width=9, palette = pal, save_plot="Colombia_dyadic_1.pdf")
multiplex_plot(fit_1, type="generalized", HPDI=0.9, plot = TRUE, height=7, width=9, palette = pal, save_plot="Colombia_generalized_1.pdf")




pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,2))
longitudinal_plot(fit_1, type="dyadic", save_plot="Colombia_dyadic_long_1.pdf", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,2,9,3,4))
longitudinal_plot(fit_1, type="generalized", save_plot="Colombia_generalized_long_1.pdf", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,2,9,7,5,10,8,6))
longitudinal_plot(fit_1,type="coefficient", 
    parameter_set = list(
    focal="Age", target="Age", 
    focal="Wealth", target="Wealth", 
    focal="SexMale", target="SexMale",
    dyadic="Friend", dyadic="Relatedness"),
    palette=pal,
    normalized=TRUE,
    height=4, width=9,
    save_plot="Slopes_Colombia.pdf")

# Model 1 - Plot results
strand_caterpillar_plot(res_1, 
 submodels=c(
 "Focal effects: Out-degree",
 "Target effects: In-degree",
 "Dyadic effects"), 
 export_as_table = FALSE, 
 normalized=FALSE
)


block_pars = rbind(
process_block_parameters(fit_1, "Afrocolombian to Afrocolombian", "Afrocolombian to Embera", HPDI=0.9),
process_block_parameters(fit_1, "Embera to Embera", "Afrocolombian to Embera", HPDI=0.9),
process_block_parameters(fit_1, "Embera to Afrocolombian", "Afrocolombian to Embera", HPDI=0.9)
)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,3,7))
longitudinal_plot(fit_1, type="custom", results=block_pars, plot = TRUE, save_plot = "Block_params.pdf", height=6, width=6, palette=pal)

