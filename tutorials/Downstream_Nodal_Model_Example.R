###########################################################################################################
#
#   Estimating the effect of nodal out- and in-strength on downstream variables
#
###########################################################################################################
set.seed(1)

# Load libraries
library(STRAND)
library(ggplot2)
library(igraph)
library(MASS)

############################################################################## Part 1: Network Analysis
# Make network data
N_id = 180

# Covariates
Mass = rbern(N_id, 0.4)
Age = standardize_strand(runif(N_id, 5, 50))
individual = data.frame(Mass=Mass, Age=Age)

V = 1            # One blocking variable
G = 1            # One category

clique = rep(1, N_id)
B = matrix(-8, nrow=G, ncol=G)
diag(B) = -4.5 # Block matrix

G = simulate_sbm_plus_srm_network(N_id = N_id, B=list(B=B), V=V, 
                         groups=data.frame(clique=factor(clique)),
                         sr_sigma = c(2.4, 1.8), sr_rho = -0.75,
                         dr_sigma = 1.9, dr_rho = 0.8,
                         outcome_mode="bernoulli"
                               )

Net = G$network

################################################### Organize for model fitting
# Add row and colnames
name_vec = paste("Individual", 1:N_id)
rownames(individual) = rownames(Net) = colnames(Net) = name_vec


model_dat = make_strand_data(outcome=list(Gave = Net),  
                             block_covariates=NULL, 
                             individual_covariates=individual, 
                             dyadic_covariates=NULL,  
                             outcome_mode = "bernoulli", 
                             link_mode="logit",
                             check_standardization = FALSE)

########################################## Now fit network STRAND 
fit1 = fit_block_plus_social_relations_model(
    data=model_dat,
      block_regression = ~ 1,
      focal_regression = ~ 1,
      target_regression = ~ 1,
      dyad_regression = ~ 1,
    mode="mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 500 ,
      iter_sampling = 500 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )

res1 = summarize_strand_results(fit1)





############################################################################## Part 2: Downstream Analysis
################################################################### a) Gaussian model
set.seed(1)
########## Simulate an outcome from nodal random effects and other variables
B = c(-2, 0.5, 1.9, 0.6, -1.5)
RS = c()

out_strength = standardize_strand(G$sr[,1])
in_strength = standardize_strand(G$sr[,2])

for(i in 1:length(Age)){
  RS[i] = rnorm(1, (B[1] + B[2]*out_strength[i] + B[3]*in_strength[i] + B[4]*Age[i] + B[5]*Mass[i]), 0.35)
}

model_dat$individual_predictors$RS = RS

######################################## Now downstream regression
fit1d = fit_downstream_nodal_model(
    fit = fit1, 
    data = model_dat,
    outcome = "RS",
    downstream_regression = ~ Age + Mass,
    nodal_effects = "both",
    outcome_mode = "gaussian",
    link_mode = "identity",
    mode = "mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500 ,
      iter_sampling = 1500 ,
      max_treedepth = 13,
      refresh = 1,
      adapt_delta = 0.98)
  )

res1d = summarize_strand_results(fit1d)
B 

####################### Compare to lm
OS = standardize_strand(rowSums(Net))
IS = standardize_strand(colSums(Net))
summary(lm(RS ~ OS + IS + Age + Mass ))
B



################################################################### b) Bernoulli model
set.seed(7)
########## Simulate an outcome from nodal random effects and other variables
B = c(-1.7, 0.7, -1.1, 0.70, -1.25)
RS = c()

out_strength = standardize_strand(G$sr[,1])
in_strength = standardize_strand(G$sr[,2])

for(i in 1:length(Age)){
  RS[i] = rbinom(1, prob=inv_logit(B[1] + B[2]*out_strength[i] + B[3]*in_strength[i] + B[4]*Age[i] + B[5]*Mass[i]), size=1)
}

model_dat$individual_predictors$RS = RS

######################################## Now downstream regression
fit1d = fit_downstream_nodal_model(
    fit = fit1, 
    data = model_dat,
    outcome = "RS",
    downstream_regression = ~ Age + Mass,
    nodal_effects = "both",
    outcome_mode = "bernoulli",
    link_mode = "logit",
    mode = "mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500 ,
      iter_sampling = 1500 ,
      max_treedepth = 13,
      refresh = 1,
      adapt_delta = 0.98)
  )

res1d = summarize_strand_results(fit1d)
B 

####################### Compare to glm
OS = standardize_strand(rowSums(Net))
IS = standardize_strand(colSums(Net))
summary(glm(RS ~ OS + IS + Age + Mass, family="binomial"))
B 




################################################################### c) Binomial model
set.seed(13)
########## Simulate an outcome from nodal random effects and other variables
B = c(-2, 1.7, -1.99, 0.6, -0.75)
RS = c()
Exposure = c()

out_strength = standardize_strand(G$sr[,1])
in_strength = standardize_strand(G$sr[,2])

for(i in 1:length(Age)){
  Exposure[i] = rpois(1, 30)
  RS[i] = rbinom(1, prob=inv_logit(B[1] + B[2]*out_strength[i] + B[3]*in_strength[i] + B[4]*Age[i] + B[5]*Mass[i]), size=Exposure[i])
}

model_dat$individual_predictors$RS = RS
model_dat$individual_predictors$Exposure = Exposure

######################################## Now downstream regression
fit1d = fit_downstream_nodal_model(
    fit = fit1, 
    data = model_dat,
    outcome = "RS",
    exposure = "Exposure",
    downstream_regression = ~ Age + Mass,
    nodal_effects = "both",
    outcome_mode = "binomial",
    link_mode = "logit",
    mode = "mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500 ,
      iter_sampling = 1500 ,
      max_treedepth = 13,
      refresh = 1,
      adapt_delta = 0.98)
  )

res1d = summarize_strand_results(fit1d)
B 

####################### Compare to glm
OS = standardize_strand(rowSums(Net))
IS = standardize_strand(colSums(Net))
summary(glm(formula = cbind(RS, Exposure - RS) ~ OS + IS + Age + Mass, family = "binomial"))
B




################################################################### d) Poisson model
set.seed(42)
########## Simulate an outcome from nodal random effects and other variables
B = c(-1, 0.7, -0.99, 0.6, -0.95)
RS = c()

out_strength = standardize_strand(G$sr[,1])
in_strength = standardize_strand(G$sr[,2])

for(i in 1:length(Age)){
  RS[i] = rpois(1, lambda=exp(B[1] + B[2]*out_strength[i] + B[3]*in_strength[i] + B[4]*Age[i] + B[5]*Mass[i]))
}

model_dat$individual_predictors$RS = RS

######################################## Now downstream regression
fit1d = fit_downstream_nodal_model(
    fit = fit1, 
    data = model_dat,
    outcome = "RS",
    downstream_regression = ~ Age + Mass,
    nodal_effects = "both",
    outcome_mode = "poisson",
    link_mode = "log",
    mode = "mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500 ,
      iter_sampling = 1500 ,
      max_treedepth = 13,
      refresh = 1,
      adapt_delta = 0.98)
  )

res1d = summarize_strand_results(fit1d)
B 

####################### Compare to glm
OS = standardize_strand(rowSums(Net))
IS = standardize_strand(colSums(Net))
summary(glm(RS ~ OS + IS + Age + Mass, family="poisson"))
B



################################################################### e) Negative binomial model
set.seed(666)
########## Simulate an outcome from nodal random effects and other variables
B = c(1.0, 1.1, -0.99, 0.6, -0.75)
RS = c()

rgampois = function (n, mu, scale){
    shape = mu/scale
    prob = 1/(1 + scale)
    rnbinom(n, size = shape, prob = prob)
}

out_strength = standardize_strand(G$sr[,1])
in_strength = standardize_strand(G$sr[,2])

for(i in 1:length(Age)){
  RS[i] = rgampois(1, mu=exp(B[1] + B[2]*out_strength[i] + B[3]*in_strength[i] + B[4]*Age[i] + B[5]*Mass[i]), scale=2)
}

model_dat$individual_predictors$RS = RS

######################################## Now downstream regression
fit1d = fit_downstream_nodal_model(
    fit = fit1, 
    data = model_dat,
    outcome = "RS",
    downstream_regression = ~ Age + Mass,
    nodal_effects = "both",
    outcome_mode = "negative_binomial",
    link_mode = "log",
    mode = "mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500 ,
      iter_sampling = 1500 ,
      max_treedepth = 13,
      refresh = 1,
      adapt_delta = 0.98)
  )

res1d = summarize_strand_results(fit1d)
B 

####################### Compare to glm
library(MASS)
OS = standardize_strand(rowSums(Net))
IS = standardize_strand(colSums(Net))
summary(glm.nb(RS ~ OS + IS + Age + Mass))
B



################################################################### f) Beta model
set.seed(1337)
########## Simulate an outcome from nodal random effects and other variables
B = c(-1, 0.7, -0.8, 0.6, -0.75)
RS = c()

out_strength = standardize_strand(G$sr[,1])
in_strength = standardize_strand(G$sr[,2])
scale = 90

for(i in 1:length(Age)){
  mu = inv_logit(B[1] + B[2]*out_strength[i] + B[3]*in_strength[i] + B[4]*Age[i] + B[5]*Mass[i])
  RS[i] = rbeta(1, mu*scale, (1-mu)*scale)
}

model_dat$individual_predictors$RS = RS


######################################## Now downstream regression
fit1d = fit_downstream_nodal_model(
    fit = fit1, 
    data = model_dat,
    outcome = "RS",
    downstream_regression = ~ Age + Mass,
    nodal_effects = "both",
    outcome_mode = "beta",
    link_mode = "logit",
    mode = "mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500 ,
      iter_sampling = 1500 ,
      max_treedepth = 13,
      refresh = 1,
      adapt_delta = 0.98)
  )

res1d = summarize_strand_results(fit1d)
B 

####################### Compare to glm
library(betareg)
OS = standardize_strand(rowSums(Net))
IS = standardize_strand(colSums(Net))
summary(betareg(RS ~ OS + IS + Age + Mass))
B




################################################################### g) Gamma model
set.seed(8675309)
########## Simulate an outcome from nodal random effects and other variables
B = c(-1, 0.7, -0.4, 0.6, -0.75)
RS = c()

out_strength = standardize_strand(G$sr[,1])
in_strength = standardize_strand(G$sr[,2])
scale = 20

for(i in 1:length(Age)){
  mu = exp(B[1] + B[2]*out_strength[i] + B[3]*in_strength[i] + B[4]*Age[i] + B[5]*Mass[i])
  RS[i] = rgamma(1, mu*scale, scale)
}

model_dat$individual_predictors$RS = RS + 0.01 # Gotta be hacky


######################################## Now downstream regression
fit1d = fit_downstream_nodal_model(
    fit = fit1, 
    data = model_dat,
    outcome = "RS",
    downstream_regression = ~ Age + Mass,
    nodal_effects = "both",
    outcome_mode = "gamma",
    link_mode = "log",
    mode = "mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500 ,
      iter_sampling = 1500 ,
      max_treedepth = 13,
      refresh = 1,
      adapt_delta = 0.98)
  )

res1d = summarize_strand_results(fit1d)
B 

####################### Compare to glm
OS = standardize_strand(rowSums(Net))
IS = standardize_strand(colSums(Net))
summary(glm(RS ~ OS + IS + Age + Mass, family = Gamma(link = "log"))) 
B





