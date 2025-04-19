###################################################
#
#   Multiplex Poisson analyses with data simulation 
#
########################################

# Clear working space
rm(list = ls())
set.seed(50)
# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors
 colors = plvs_vltra("dust_storm", rev=FALSE, elements=NULL, show=FALSE)

library(stringr)
library(ggplot2)
library(psych)
library(rethinking)
library(STRAND)

# Make data
 N_id = 65      # Individuals in network
 N_layers = 3   # Network layers

# Covariates
 Kinship = standardize(rlkjcorr( 1 , N_id , eta=1.5 ))
 Dominance = ceiling(rlkjcorr( 1 , N_id , eta=1.5 ) - 0.1)
 Mass = rbern(N_id, 0.4)
 Age = rnorm(N_id, 0, 1)
 Strength = rnorm(N_id, 0, 1)

# Organize into list
 dyadic_preds = array(NA,c(N_id,N_id,2))

 dyadic_preds[,,1] = Kinship
 dyadic_preds[,,2] = Dominance

# Set effect sizes
sr_mu = rep(0, N_layers*2)                      # Average effect size, should be zero
sr_sigma = c(0.2, 0.7, 0.1, 0.3, 0.7, 0.8)      # Variation in random effects. First 3 are sender effects, one per layer. Last 3 are receiver effects.
sr_Rho = rlkjcorr( 1 , N_layers*2 , eta=1.5 )   # Generalized reciprocity matrix.
dr_mu = rep(0, N_layers)                        # Average effect size, should be zero
dr_sigma = c(0.9, 1.1, 1.2)                     # Variation in dyadic random effects.

# Build dyadic matrix in full
dr_Rho = structure(c(1, -0.226068357721697, 0.362134589509303, 0.247782279255332, 
0.116059835333589, -0.116688891170166, -0.226068357721697, 1, 
0.623976731189517, 0.110006148152744, 0.00712660302928801, 0.60311693469165, 
0.362134589509303, 0.623976731189517, 1, -0.11764761911035, 0.612649591988435, 
0.628827889254145, 0.247782279255332, 0.110006148152744, -0.11764761911035, 
1, -0.226131976783832, 0.362199967761652, 0.116059835333589, 
0.00712660302928801, 0.612649591988435, -0.226131976783832, 1, 
0.624070637620319, -0.116688891170166, 0.60311693469165, 0.628827889254145, 
0.362199967761652, 0.624070637620319, 1), dim = c(6L, 6L))   

chol(dr_Rho)                                    # Check if positive definite

# Covariate effects
sr_1 = matrix(NA, nrow=2, ncol=3)               # Layer 1 
sr_1[1,] = c(-0.5, 1.0, -0.7)                    # Effect of Mass, Age, and Strength on out-degree
sr_1[2,] = c(0.7, -1.1, -1)                      # Effect of Mass, Age, and Strength on in-degree

sr_2 = matrix(NA, nrow=2, ncol=3)               # Layer 2
sr_2[1,] = c(-0.1, -1.0, 0.7)                   # Effect of Mass, Age, and Strength on out-degree
sr_2[2,] = c(-0.7, -0.6, -0.01)                  # Effect of Mass, Age, and Strength on in-degree

sr_3 = matrix(NA, nrow=2, ncol=3)               # Layer 3
sr_3[1,] = c(0.1, 0.3, -0.5)                    # Effect of Mass, Age, and Strength on out-degree
sr_3[2,] = c(0.1, 0.4, 0)                       # Effect of Mass, Age, and Strength on in-degree

sr_effects = list(sr_1, sr_2, sr_3)             # Organize into list

dr_effects = list(c(0.6, 0.3),                  # Layer 1 effect of Kinship and Dominant        
                  c(-0.2, -0.7),                # Layer 2 effect of Kinship and Dominant 
                  c(0.2, -0.7))                 # Layer 3 effect of Kinship and Dominant 

# Block structure
group_probs_block_size = c(0.25, c(0.25, 0.75)*(1-0.25))
groups_1 = rep("Any",N_id) 
groups_2 = sample( c("Mottled","Striped","Spotted") , size=N_id , replace=TRUE , prob=group_probs_block_size )
groups_3 = sample( c("Male", "Female") , size=N_id , replace=TRUE , prob=c(0.5,0.5) )

# Intercept in each layer
B_1a = matrix(-2.7,nrow=1,ncol=1)
B_2a = matrix(0.5,nrow=1,ncol=1)
B_3a = matrix(-2.5,nrow=1,ncol=1)

# Offsets for Pattern
B_1b = matrix(rnorm(9,0,1),nrow=3,ncol=3)
B_2b = matrix(rnorm(9,-2,1),nrow=3,ncol=3)
B_3b = matrix(rnorm(9,0.5,1),nrow=3,ncol=3)

diag(B_1b) = diag(B_1b) - 0
diag(B_2b) = diag(B_2b) + 0.5
diag(B_3b) = diag(B_3b) + 0

# Offset for Sex
B_1c = matrix(rnorm(4,0,1),nrow=2,ncol=2)
B_2c = matrix(rnorm(4,-3,1),nrow=2,ncol=2)
B_3c = matrix(rnorm(4,-2,1),nrow=2,ncol=2)

diag(B_1c) = diag(B_1c) - 0.5
diag(B_2c) = diag(B_2c) + 1
diag(B_3c) = diag(B_3c) + 0.2

# Merge into lists
B1 = list(B_1a, B_1b, B_1c)
B2 = list(B_2a, B_2b, B_2c)
B3 = list(B_3a, B_3b, B_3c)

B = list(B1, B2, B3)
 
groups = data.frame(Intercept=as.numeric(factor(groups_1)), Pattern=as.numeric(factor(groups_2)), Sex=as.numeric(factor(groups_3)))
groups_f = data.frame(Intercept=factor(groups_1), Pattern=factor(groups_2), Sex=factor(groups_3))

#################################################### Simulate network
G = simulate_multiplex_network(
  N_id = N_id,            
  N_layers = N_layers,                   
  B = B,                       
  V = 3,       
  groups = groups,                     
  sr_mu = sr_mu,            
  sr_sigma = sr_sigma,                        
  sr_Rho = sr_Rho,                     
  dr_mu = dr_mu,                            
  dr_sigma = dr_sigma,                         
  dr_Rho = dr_Rho,                          
  outcome_mode="poisson",    
  link_mode = "log",            
  individual_predictors = data.frame(Mass=Mass, 
                                     Age=Age, 
                                     Strength=Strength),    
  dyadic_predictors = dyadic_preds,        
  individual_effects = sr_effects,        
  dyadic_effects = dr_effects           
 )

table(G$network[1,,])
table(G$network[2,,])
table(G$network[3,,])

#################################################### Create the STRAND data object
outcome = list(Feeding = G$network[1,,], Fighting = G$network[2,,], Grooming = G$network[3,,])
exposure = list(Feeding = G$exposure[1,,], Fighting = G$exposure[2,,], Grooming = G$exposure[3,,])

dyad = list(Kinship = Kinship, 
            Dominance = Dominance
            )

groups = data.frame(Pattern=factor(groups_2), 
                    Sex=factor(groups_3)
                    )

indiv =  data.frame(Mass=Mass, 
                    Age=Age, 
                    Strength=Strength
                     )

### col and row names are now a soft requirement
# can turn off with check_data_organization = FALSE, but its reccmended to always run checks on row and col names
labels = paste("Ind", 1:N_id)
colnames(outcome$Feeding) = rownames(outcome$Feeding) = labels
colnames(outcome$Fighting) = rownames(outcome$Fighting) = labels
colnames(outcome$Grooming) = rownames(outcome$Grooming) = labels

colnames(exposure$Feeding) = rownames(exposure$Feeding) = labels
colnames(exposure$Fighting) = rownames(exposure$Fighting) = labels
colnames(exposure$Grooming) = rownames(exposure$Grooming) = labels

colnames(dyad$Kinship) = rownames(dyad$Kinship) = labels
colnames(dyad$Dominance) = rownames(dyad$Dominance) = labels

rownames(indiv) = labels
rownames(groups) = labels


############# Build data
dat0 = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       exposure = exposure,
                       outcome_mode="poisson",
                       link_mode="log",
                       multiplex = TRUE,
                       imputation=TRUE)

############# Add some NAs to data
indiv$Age[c(7, 13, 30:35, 42, 46+2)] = NA
indiv$Mass[c(2,4,6,8,10,12,14,16,18,20)] = NA
dyad$Kinship[c(3, 7, 14, 59), c(15:20, 30:35)] = NA
groups$Sex[c(17, 39, 22:25, 42)] = NA

outcome$Feeding[c(3, 17, 24, 51:60), c(17:22, 23:49)] = NA
outcome$Feeding[c(40:60), c(19:22, 27:39)] = NA

outcome$Fighting[c(51:60), c(17:22, 23:49)] = NA
outcome$Fighting[c(5, 7, 14, 59:60), c(19:22, 27:39)] = NA

outcome$Grooming[c(3, 17, 24, 51:60), c(17:28, 29:49)] = NA
outcome$Grooming[c(5, 7, 14, 59:60), c(19:22, 27:39)] = NA

dat = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       exposure = exposure,
                       outcome_mode="poisson",
                       link_mode="log",
                       multiplex = TRUE,
                       imputation=TRUE)

# Model
fit1a = fit_multiplex_model_missings(data=dat,
                          block_regression = ~ Pattern + Sex,
                          focal_regression = ~ Mass + Age + Strength,
                          target_regression = ~ Mass + Age + Strength,
                          dyad_regression = ~ Kinship + Dominance,
                          bandage_penalty = 0.1,
                          mode="mcmc",
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = NULL, adapt_delta = 0.95)
)

# Model
fit2a = fit_multiplex_model_missings(data=dat0,
                          block_regression = ~ Pattern + Sex,
                          focal_regression = ~ Mass + Age + Strength,
                          target_regression = ~ Mass + Age + Strength,
                          dyad_regression = ~ Kinship + Dominance,
                          bandage_penalty = 0.1,
                          mode="mcmc",
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = NULL, adapt_delta = 0.95)
)

# Model
fit0a = fit_multiplex_model(data=dat0,
                          block_regression = ~ Pattern + Sex,
                          focal_regression = ~ Mass + Age + Strength,
                          target_regression = ~ Mass + Age + Strength,
                          dyad_regression = ~ Kinship + Dominance,
                          bandage_penalty = 0.1,
                          mode="mcmc",
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = NULL, adapt_delta = 0.95)
)


# Model
fit1b = fit_multiplex_model_missings(data=dat,
                          block_regression = ~ Pattern + Sex,
                          focal_regression = ~ Mass + Age + Strength,
                          target_regression = ~ Mass + Age + Strength,
                          dyad_regression = ~ Kinship + Dominance,
                          bandage_penalty = -1,
                          mode="mcmc",
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1, init=0,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = NULL, adapt_delta = 0.95)
)

# Model
fit2b = fit_multiplex_model_missings(data=dat0,
                          block_regression = ~ Pattern + Sex,
                          focal_regression = ~ Mass + Age + Strength,
                          target_regression = ~ Mass + Age + Strength,
                          dyad_regression = ~ Kinship + Dominance,
                          bandage_penalty = -1,
                          mode="mcmc",
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1, init=0,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = NULL, adapt_delta = 0.95)
)

# Model
fit0b = fit_multiplex_model(data=dat0,
                          block_regression = ~ Pattern + Sex,
                          focal_regression = ~ Mass + Age + Strength,
                          target_regression = ~ Mass + Age + Strength,
                          dyad_regression = ~ Kinship + Dominance,
                          bandage_penalty = -1,
                          mode="mcmc",
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1, init=0,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = NULL, adapt_delta = 0.95)
)


######################################## Estimates
# Summaries
res1a = summarize_strand_results(fit1a)
res2a = summarize_strand_results(fit2a)
res0a = summarize_strand_results(fit0a)

res1b = summarize_strand_results(fit1b)
res2b = summarize_strand_results(fit2b)
res0b = summarize_strand_results(fit0b)


# Plots
vis_1a = strand_caterpillar_plot(res1a, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-Imputation-L2")
vis_2a = strand_caterpillar_plot(res2a, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-No-Imputation-L2")
vis_0a = strand_caterpillar_plot(res0a, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="Old-No-Imputation-L2")

vis_1b = strand_caterpillar_plot(res1b, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-Imputation-Pinkney")
vis_2b = strand_caterpillar_plot(res2b, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="New-No-Imputation-Pinkney")
vis_0b = strand_caterpillar_plot(res0b, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, only_slopes=TRUE, export_as_table=TRUE, site="Old-No-Imputation-Pinkney")

df = rbind(vis_1a, vis_2a, vis_0a, vis_1b, vis_2b, vis_0b)
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


