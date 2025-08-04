#####################################################################
#
#   Computing posterior network metrics from STRAND's latent network
#
########################################

# Clear working space
rm(list = ls())
set.seed(1)
# Load libraries
library(rethinking)
library(STRAND)
library(ggplot2)
library(igraph)


# Make data
N_id = 90

# Covariates
Kinship = standardize(rlkjcorr( 1 , N_id , eta=1.5 ))
Dominant = ceiling(rlkjcorr( 1 , N_id , eta=1.5 ) - 0.1)
Mass = rbern(N_id, 0.4)

# Organize into list
dyadic_preds = array(NA,c(N_id,N_id,2))

dyadic_preds[,,1] = Kinship
dyadic_preds[,,2] = Dominant

# Set effect sizes
sr_mu = c(0,0)  
sr_sigma = c(1.2, 1.7) 
sr_rho = 0.55
dr_mu = 0 
dr_sigma = 2.5
dr_rho= 0.85
sr_effects_1 = c(1.9, 1.3)
dr_effects_1 = c(-1.2, 1.7)

# Block structure
group_probs_block_size = c(0.25, c(0.25, 0.25)*(1-0.25))

B_1 = matrix(-10,nrow=1,ncol=1)
B_2 = matrix(rnorm(9,0,3),nrow=3,ncol=3)
B_3 = matrix(rnorm(4,0,3),nrow=2,ncol=2)

diag(B_2) = diag(B_2) + 2
diag(B_3) = diag(B_3) + 3.5

B=list(B_1, B_2, B_3)
 
groups_1 = rep("Any",N_id) 
groups_2 = sample( c("Red","White","Blue") , size=N_id , replace=TRUE , prob=group_probs_block_size )
groups_3 = sample( c("Strange", "Charm") , size=N_id , replace=TRUE , prob=c(0.5,0.5) )

groups = data.frame(Intercept=as.numeric(factor(groups_1)), Merica=as.numeric(factor(groups_2)), Quantum=as.numeric(factor(groups_3)))
groups_f = data.frame(Intercept=factor(groups_1), Merica=factor(groups_2), Quantum=factor(groups_3))
individual = data.frame(Mass=Mass)

#################################################### Simulate SBM + SRM network
G = simulate_sbm_plus_srm_network(N_id = N_id, 
                         B = B, 
                         V=3,
                         groups=groups,                  
                         sr_mu = sr_mu,  
                         sr_sigma = sr_sigma, 
                         sr_rho = sr_rho,
                         dr_mu = dr_mu,  
                         dr_sigma = dr_sigma, 
                         dr_rho = dr_rho,
                         outcome_mode="bernoulli", 
                         link_mode="logit",                 
                         individual_predictors = data.frame(Mass=Mass),
                         dyadic_predictors = dyadic_preds,
                         individual_effects = matrix(sr_effects_1,nrow=2,ncol=1),
                         dyadic_effects = dr_effects_1
                         )   

Net = G$network

################################################### Organize for model fitting
# Add row and colnames
name_vec = paste("Individual", 1:N_id)
rownames(Net) = colnames(Net) = name_vec
rownames(Kinship) = colnames(Kinship) = name_vec
rownames(Dominant) = colnames(Dominant) = name_vec
rownames(groups_f) = name_vec
rownames(individual) = name_vec

model_dat = make_strand_data(outcome=list(Gave = Net),  
                             block_covariates=groups_f, 
                             individual_covariates=individual, 
                             dyadic_covariates=list(Kinship=Kinship, Dominant=Dominant),  
                             outcome_mode = "bernoulli", 
                             link_mode="logit")

########################################## Now fit STRAND 
fit1 = fit_block_plus_social_relations_model(
    data=model_dat,
      block_regression = ~ Merica + Quantum,
      focal_regression = ~ Mass,
      target_regression = ~ Mass,
      dyad_regression = ~ Kinship + Dominant,
    mode="mcmc",
    return_predicted_network = TRUE,   # These samples take up a lot of space, but we need them here
    stan_mcmc_parameters = list(
      chains = 1,
      iter_warmup = 500 ,
      iter_sampling = 500 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )


res1 = summarize_strand_results(fit1)

################################ Extract predicted network of latent tie probabilities
 latent_network = res1$samples$predicted_network_sample
 
######################### Simulate a Bernoulli outcome
 bern_out = function(x, mode){
  y = matrix(0, nrow=nrow(x), ncol=ncol(x)) 
  for(i in 1:nrow(x)){
   for(j in 1:nrow(x)){
    y[i,j] = rbinom(1, size=1, prob=x[i,j]) 
   }
  }

  if(mode=="obs"){
   return(y)
  }

  if(mode=="latent"){
   return(x)
  }
 }

 
#################################### Compute network metrics for each MCMC sample
################### Observed network
 net_lab0 = net_res0 = matrix(NA, nrow=1, ncol=6)

 net = igraph::graph_from_adjacency_matrix(bern_out(Net, mode="obs"), mode = "directed")

 i = 1
 node_select = 49 # Pick a random node to show posterior quantities for that node

# Three graph-level metrics
 net_res0[i,1] = igraph::centr_betw(net)$centralization
 net_lab0[i,1] = "Centralization"

 net_res0[i,2] = igraph::reciprocity(net)
 net_lab0[i,2] = "Reciprocity"

 net_res0[i,3] = igraph::transitivity(net, type = "average")
 net_lab0[i,3] = "Transitivity"

# Three node-level metrics
 net_res0[i,4] = igraph::betweenness(net, normalize = TRUE)[node_select]
 net_lab0[i,4] = "Betweenness"

 net_res0[i,5] = igraph::degree(net, mode = "in")[node_select]
 net_lab0[i,5] = "Degree"

 net_res0[i,6] = igraph::eigen_centrality(net,scale = TRUE)$vector[node_select]
 net_lab0[i,6] = "Centrality"
 
dat_long0 = data.frame(Value = c(net_res0), Metric = factor(c(net_lab0))) 

######################### Posterior
 net_lab = net_res = matrix(NA, nrow=500, ncol=6)

 # Loop over MCMC samples and compute network metrics on each predicted network
 for(i in 1:500){
  net = igraph::graph_from_adjacency_matrix(bern_out(latent_network[i,,],mode="obs"), mode = "directed")

 net_res[i,1] = igraph::centr_betw(net)$centralization
 net_lab[i,1] = "Centralization"

 net_res[i,2] = igraph::reciprocity(net)
 net_lab[i,2] = "Reciprocity"

 net_res[i,3] = igraph::transitivity(net, type = "average")
 net_lab[i,3] = "Transitivity"

 net_res[i,4] = igraph::betweenness(net, normalize = TRUE)[node_select]
 net_lab[i,4] = "Betweenness"

 net_res[i,5] = igraph::degree(net, mode = "in")[node_select]
 net_lab[i,5] = "Degree"

 net_res[i,6] = igraph::eigen_centrality(net,scale = TRUE)$vector[node_select]
 net_lab[i,6] = "Centrality"
 }

dat_long = data.frame(Value = c(net_res), Metric = factor(c(net_lab))) 


############ Plot
dat_long$Metric = factor(dat_long$Metric, levels=c("Centralization","Reciprocity","Transitivity","Betweenness","Degree","Centrality"))
dat_long0$Metric = factor(dat_long0$Metric, levels=c("Centralization","Reciprocity","Transitivity","Betweenness","Degree","Centrality"))

ggplot(dat_long, aes(Value)) +
  geom_density(fill = "grey40", color=NA) + 
  geom_vline(data=dat_long0,aes(xintercept=Value)) +
  facet_wrap(~ Metric, scale="free")
  
