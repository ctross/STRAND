###############################################################################################
#
#   Multiplex analyses with JAX 
#
###############################################################################################

# Load libraries
 library(reticulate)
 library(posterior)
 library(igraph)
 library(STRAND)
 library(ggplot2)
 library(stringr)
 library(psych)

 set.seed(1)

# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors
 colors = plvs_vltra("dust_storm", rev=FALSE, elements=NULL, show=FALSE)

# Make data
 N_id = 50       # Individuals in network
 N_layers = 7    # Network layers

# Covariates
 Kinship = standardize_strand(rlkjcorr( 1 , N_id , eta=1.5 ))
 Dominance = ceiling(rlkjcorr( 1 , N_id , eta=1.5 ) - 0.1)
 Mass = rbern(N_id, 0.4)
 Age = rnorm(N_id, 0, 1)
 Strength = rnorm(N_id, 0, 1)

# Organize into list
 dyadic_preds = array(NA,c(N_id,N_id,2))

 dyadic_preds[,,1] = Kinship
 dyadic_preds[,,2] = Dominance

# Set effect sizes
 sr_mu = rep(0, N_layers*2)                                          # Average effect size, should be zero
 sr_sigma = c(0.2, 0.7, 1.1, 0.3, 0.7, 1.8, 1.5,                     # Variation in random effects. First 7 are sender effects, one per layer. 
              0.5, 1.7, 0.7, 0.9, 0.7, 0.8, 1.5)                     # Last 7 are receiver effects.
 sr_Rho = rlkjcorr( 1 , N_layers*2 , eta=1.5 )                       # Generalized reciprocity matrix.
 dr_mu = rep(0, N_layers)                                            # Average effect size, should be zero
 dr_sigma = c(0.9, 1.1, 1.2, 0.5, 1.6, 2.1, 2.9)                     # Variation in dyadic random effects.
 error_sigma = c(0.35, 2.1, 1.2, 0.2, 0.9, 1.5, 2.5)                 # Variation in noise effects.

 data = NULL
 data$N_responses = N_layers
 data$bandage_penalty = 0.01

 dr_gen = generate_structured_correlation_matrix(
  data,
  eta = 1,
  setting = "multiplex_dyadic_reciprocity",
  mode = "l2norm",
  mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1, seed=67, init=0,
                                                        iter_warmup = 2000, iter_sampling = 20,
                                                        max_treedepth = 12, adapt_delta = 0.95))

 stanfit = posterior::as_draws_rvars(dr_gen$fit$draws())
 D_corr = posterior::draws_of(stanfit$"D_corr")

# Build dyadic matrix in full
 dr_Rho = D_corr[1,,]   

 chol(dr_Rho)                                    # Check if positive definite

# Covariate effects
 sr_1 = matrix(NA, nrow=2, ncol=3)                # Layer 1 
 sr_1[1,] = c(-0.5, 1.0, -0.7)                    # Effect of Mass, Age, and Strength on out-degree
 sr_1[2,] = c(0.7, -1.1, -1)                      # Effect of Mass, Age, and Strength on in-degree

 sr_2 = matrix(NA, nrow=2, ncol=3)                # Layer 2
 sr_2[1,] = c(-0.1, -1.0, 0.7)                    # Effect of Mass, Age, and Strength on out-degree
 sr_2[2,] = c(-0.7, -0.6, -0.01)                  # Effect of Mass, Age, and Strength on in-degree

 sr_3 = matrix(NA, nrow=2, ncol=3)                # Layer 3
 sr_3[1,] = c(0.1, 0.3, -0.5)                     # Effect of Mass, Age, and Strength on out-degree
 sr_3[2,] = c(0.1, 0.4, 0)                        # Effect of Mass, Age, and Strength on in-degree

 sr_4 = matrix(NA, nrow=2, ncol=3)                # Layer 4
 sr_4[1,] = c(0.7, 1.3, -1.5)                     # Effect of Mass, Age, and Strength on out-degree
 sr_4[2,] = c(0.2, 0.4, 1.0)                      # Effect of Mass, Age, and Strength on in-degree

 sr_5 = matrix(NA, nrow=2, ncol=3)                # Layer 5
 sr_5[1,] = c(0.1, -0.3, 0.5)                     # Effect of Mass, Age, and Strength on out-degree
 sr_5[2,] = c(-0.6, 2.4, -1.3)                    # Effect of Mass, Age, and Strength on in-degree

 sr_6 = matrix(NA, nrow=2, ncol=3)                # Layer 6
 sr_6[1,] = c(1.1, 0.8, -2.5)                     # Effect of Mass, Age, and Strength on out-degree
 sr_6[2,] = c(2.1, 0.8, 0)                        # Effect of Mass, Age, and Strength on in-degree

 sr_7 = matrix(NA, nrow=2, ncol=3)                # Layer 7
 sr_7[1,] = c(1.1, 0, -0.9)                       # Effect of Mass, Age, and Strength on out-degree
 sr_7[2,] = c(1.1, 0, 0.9)                        # Effect of Mass, Age, and Strength on in-degree

 sr_effects = list(sr_1, sr_2, sr_3, sr_4, sr_5, sr_6, sr_7)  # Organize into list

 dr_effects = list(c(0.6, 0.3),                 # Layer 1 effect of Kinship and Dominant        
                  c(-0.2, -0.7),                # Layer 2 effect of Kinship and Dominant
                  c(-1.1, 1.7),                 # Layer 3 effect of Kinship and Dominant 
                  c(1.2, -0.3),                 # Layer 4 effect of Kinship and Dominant 
                  c(0.2, 0.7),                  # Layer 5 effect of Kinship and Dominant 
                  c(-1.2, -1.7),                # Layer 6 effect of Kinship and Dominant  
                  c(0.2, -0.7))                 # Layer 7 effect of Kinship and Dominant 

# Block structure
 group_probs_block_size = c(0.25, c(0.25, 0.75)*(1-0.25))
 groups_1 = rep("Any",N_id) 
 groups_2 = sample( c("Mottled","Striped","Spotted") , size=N_id , replace=TRUE , prob=group_probs_block_size )
 groups_3 = sample( c("Male", "Female") , size=N_id , replace=TRUE , prob=c(0.5,0.5) )

# Intercept in each layer
 B_1a = matrix(-0.7,nrow=1,ncol=1)
 B_2a = matrix(0.5,nrow=1,ncol=1)
 B_3a = matrix(-3.5,nrow=1,ncol=1)
 B_4a = matrix(-2.5,nrow=1,ncol=1)
 B_5a = matrix(-1.5,nrow=1,ncol=1)
 B_6a = matrix(-3.5,nrow=1,ncol=1)
 B_7a = matrix(-2.5,nrow=1,ncol=1)

# Offsets for Pattern
 B_1b = matrix(rnorm(9,0,1),nrow=3,ncol=3)
 B_2b = matrix(rnorm(9,-2,1),nrow=3,ncol=3)
 B_3b = matrix(rnorm(9,0.5,1),nrow=3,ncol=3)
 B_4b = matrix(rnorm(9,0.5,1),nrow=3,ncol=3)
 B_5b = matrix(rnorm(9,0.5,1),nrow=3,ncol=3)
 B_6b = matrix(rnorm(9,0.5,1),nrow=3,ncol=3)
 B_7b = matrix(rnorm(9,0.5,1),nrow=3,ncol=3)

 diag(B_1b) = diag(B_1b) - 0
 diag(B_2b) = diag(B_2b) + 1.5
 diag(B_3b) = diag(B_3b) + 0
 diag(B_4b) = diag(B_4b) + 0
 diag(B_5b) = diag(B_5b) + 0
 diag(B_6b) = diag(B_6b) + 0
 diag(B_7b) = diag(B_7b) + 0

# Offset for Sex
 B_1c = matrix(rnorm(4,0,1),nrow=2,ncol=2)
 B_2c = matrix(rnorm(4,-3,1),nrow=2,ncol=2)
 B_3c = matrix(rnorm(4,-2,1),nrow=2,ncol=2)
 B_4c = matrix(rnorm(4,-2,1),nrow=2,ncol=2)
 B_5c = matrix(rnorm(4,-2,1),nrow=2,ncol=2)
 B_6c = matrix(rnorm(4,-2,1),nrow=2,ncol=2)
 B_7c = matrix(rnorm(4,-2,1),nrow=2,ncol=2)

 diag(B_1c) = diag(B_1c) - 0.5
 diag(B_2c) = diag(B_2c) + 1
 diag(B_3c) = diag(B_3c) + 0.2
 diag(B_4c) = diag(B_4c) - 0.2
 diag(B_5c) = diag(B_5c) - 0.2
 diag(B_6c) = diag(B_6c) + 0.2
 diag(B_7c) = diag(B_7c) - 0.2

# Merge into lists
 B1 = list(B_1a, B_1b, B_1c)
 B2 = list(B_2a, B_2b, B_2c)
 B3 = list(B_3a, B_3b, B_3c)
 B4 = list(B_4a, B_4b, B_4c)
 B5 = list(B_5a, B_5b, B_5c)
 B6 = list(B_6a, B_6b, B_6c)
 B7 = list(B_7a, B_7b, B_7c)

 B = list(B1, B2, B3, B4, B5, B6, B7)
 
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
  outcome_mode="gaussian",    
  link_mode = "identity",            
  individual_predictors = data.frame(Mass=Mass, 
                                     Age=Age, 
                                     Strength=Strength),    
  dyadic_predictors = dyadic_preds,        
  individual_effects = sr_effects,        
  dyadic_effects = dr_effects,
  error_sigma = error_sigma           
 )

 par(mfrow=c(1,7))
 hist(G$network[1,,])
 hist(G$network[2,,])
 hist(G$network[3,,])
 hist(G$network[4,,])
 hist(G$network[5,,])
 hist(G$network[6,,])
 hist(G$network[7,,])

#################################################### Create the STRAND data object
 outcome = list(Feeding = G$network[1,,], Fighting = G$network[2,,], Grooming = G$network[3,,],
                Hunting = G$network[4,,], Smilling = G$network[5,,], Barking = G$network[6,,],
                Killing = G$network[7,,])

 dyad = list(Kinship = Kinship, 
             Dominance = Dominance)

 groups = data.frame(Pattern=factor(groups_2), 
                     Sex=factor(groups_3))

 indiv = data.frame(Mass=Mass, 
                    Age=Age, 
                    Strength=Strength)

### col and row names are now a soft requirement
# can turn off with check_data_organization = FALSE, but its reccmended to always run checks on row and col names
 labels = paste("Ind", 1:N_id)
 colnames(outcome$Feeding) = rownames(outcome$Feeding) = labels
 colnames(outcome$Fighting) = rownames(outcome$Fighting) = labels
 colnames(outcome$Grooming) = rownames(outcome$Grooming) = labels
 colnames(outcome$Hunting) = rownames(outcome$Hunting) = labels
 colnames(outcome$Smilling) = rownames(outcome$Smilling) = labels
 colnames(outcome$Barking) = rownames(outcome$Barking) = labels
 colnames(outcome$Killing) = rownames(outcome$Killing) = labels

 colnames(dyad$Kinship) = rownames(dyad$Kinship) = labels
 colnames(dyad$Dominance) = rownames(dyad$Dominance) = labels

 rownames(indiv) = labels
 rownames(groups) = labels

############# Build data object
 dat = make_strand_data(outcome = outcome,
                        block_covariates = groups, 
                        individual_covariates = indiv, 
                        dyadic_covariates = dyad,
                        outcome_mode="gaussian",
                        link_mode="identity",
                        multiplex = TRUE)

####################################################### Fit with JAX
 fit_numpyro = fit_multiplex_model(data=dat,
                           block_regression = ~ Pattern + Sex,
                           focal_regression = ~ Mass + Age + Strength,
                           target_regression = ~ Mass + Age + Strength,
                           dyad_regression = ~ Kinship + Dominance,
                           mode="numpyro",
                           mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1, seed=42,
                                                 iter_warmup = 1000, iter_sampling = 1000, init = 0.25,
                                                 max_treedepth = 12, adapt_delta = 0.95,
                                                 cores=4, chain_method = "vectorized"))

 res_numpyro = summarize_strand_results(fit_numpyro)

 ####################################################### Fit with Stan
 fit_stan = fit_multiplex_model(data=dat,
                           block_regression = ~ Pattern + Sex,
                           focal_regression = ~ Mass + Age + Strength,
                           target_regression = ~ Mass + Age + Strength,
                           dyad_regression = ~ Kinship + Dominance,
                           mode="mcmc",
                           mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1, seed=67,
                                                 iter_warmup = 1000, iter_sampling = 1000,
                                                 max_treedepth = 12, adapt_delta = 0.95))

 res_stan = summarize_strand_results(fit_stan)


######################################################### Visualize results

recip_to_long = function(X){
  len_X = nrow(X)
  res = c()
  ticker = 0
     for(m in 1:(len_X-1)){
     for(n in (m+1):len_X){
        ticker = ticker + 1 
        res[ticker] = X[m,n]
        }}
 return(res)
      }

####################### First NumPyro      
df_plt = res_numpyro$summary

df_plt$TrueValues = c(
 sr_sigma[1:7],
 sr_1[1,],
 sr_2[1,],
 sr_3[1,],
 sr_4[1,],
 sr_5[1,],
 sr_6[1,],
 sr_7[1,],
 sr_sigma[8:14],
 sr_1[2,],
 sr_2[2,],
 sr_3[2,],
 sr_4[2,],
 sr_5[2,],
 sr_6[2,],
 sr_7[2,],
 dr_sigma,
 dr_effects[[1]],
 dr_effects[[2]],
 dr_effects[[3]],
 dr_effects[[4]],
 dr_effects[[5]],
 dr_effects[[6]],
 dr_effects[[7]],
 error_sigma,
 recip_to_long(sr_Rho),
 recip_to_long(dr_Rho),
 c(t(B_1a)) - mean(c(B_1a)),
 c(t(B_1b)) - mean(c(B_1b)), 
 c(t(B_1c)) - mean(c(B_1c)),
 c(t(B_2a)) - mean(c(B_2a)),
 c(t(B_2b)) - mean(c(B_2b)), 
 c(t(B_2c)) - mean(c(B_2c)),
 c(t(B_3a)) - mean(c(B_3a)),
 c(t(B_3b)) - mean(c(B_3b)), 
 c(t(B_3c)) - mean(c(B_3c)),
 c(t(B_4a)) - mean(c(B_4a)),
 c(t(B_4b)) - mean(c(B_4b)), 
 c(t(B_4c)) - mean(c(B_4c)),
 c(t(B_5a)) - mean(c(B_5a)),
 c(t(B_5b)) - mean(c(B_5b)), 
 c(t(B_5c)) - mean(c(B_5c)),
 c(t(B_6a)) - mean(c(B_6a)),
 c(t(B_6b)) - mean(c(B_6b)), 
 c(t(B_6c)) - mean(c(B_6c)),
 c(t(B_7a)) - mean(c(B_7a)),
 c(t(B_7b)) - mean(c(B_7b)), 
 c(t(B_7c)) - mean(c(B_7c))
  )

df_plt$Outcome2 = c(
 c("Feeding","Fighting","Grooming", "Hunting", "Smilling", "Barking", "Killing"),
 c("Feeding","Feeding","Feeding"),
 c("Fighting","Fighting","Fighting"),
 c("Grooming","Grooming","Grooming"),
 c("Hunting","Hunting","Hunting"),
 c("Smilling","Smilling","Smilling"),
 c("Barking","Barking","Barking"),
 c("Killing","Killing","Killing"),
 c("Feeding","Fighting","Grooming","Hunting", "Smilling", "Barking", "Killing"),
 c("Feeding","Feeding","Feeding"),
 c("Fighting","Fighting","Fighting"),
 c("Grooming","Grooming","Grooming"),
 c("Hunting","Hunting","Hunting"),
 c("Smilling","Smilling","Smilling"),
 c("Barking","Barking","Barking"),
 c("Killing","Killing","Killing"),
 c("Feeding","Fighting","Grooming","Hunting", "Smilling", "Barking", "Killing"),
 c("Feeding","Fighting","Grooming","Hunting", "Smilling", "Barking", "Killing"),
 c("Feeding","Feeding"),
 c("Fighting","Fighting"),
 c("Grooming","Grooming"),
 c("Hunting","Hunting"),
 c("Smilling","Smilling"),
 c("Barking","Barking"),
 c("Killing","Killing"),
 rep("Recip", 182),
 rep("Other", 98)
  )

df_plt$Variable2 = df_plt$Variable

colnames(df_plt) = c("Variable", "Median", "LI", "HI", "Mean", "SD", "P", "TrueValues","Outcome2", "Variable2")

df_plt$Median = as.numeric(df_plt$Median)
df_plt$LI = as.numeric(df_plt$LI)
df_plt$HI = as.numeric(df_plt$HI)


# Indentify block effect offsets
X = mean(df_plt$Median[267])
df_plt$Median[267] = df_plt$Median[267] - X
df_plt$LI[267] = df_plt$LI[267] - X
df_plt$HI[267] = df_plt$HI[267] - X

X = mean(df_plt$Median[268:276])
df_plt$Median[268:276] = df_plt$Median[268:276] - X
df_plt$LI[268:276] = df_plt$LI[268:276] - X
df_plt$HI[268:276] = df_plt$HI[268:276] - X

X = mean(df_plt$Median[277:280])
df_plt$Median[277:280] = df_plt$Median[277:280] - X
df_plt$LI[277:280] = df_plt$LI[277:280] - X
df_plt$HI[277:280] = df_plt$HI[277:280] - X


X = mean(df_plt$Median[281])
df_plt$Median[281] = df_plt$Median[281] - X
df_plt$LI[281] = df_plt$LI[281] - X
df_plt$HI[281] = df_plt$HI[281] - X

X = mean(df_plt$Median[282:290])
df_plt$Median[282:290] = df_plt$Median[282:290] - X
df_plt$LI[282:290] = df_plt$LI[282:290] - X
df_plt$HI[282:290] = df_plt$HI[282:290] - X

X = mean(df_plt$Median[291:294])
df_plt$Median[291:294] = df_plt$Median[291:294] - X
df_plt$LI[291:294] = df_plt$LI[291:294] - X
df_plt$HI[291:294] = df_plt$HI[291:294] - X


X = mean(df_plt$Median[295])
df_plt$Median[295] = df_plt$Median[295] - X
df_plt$LI[295] = df_plt$LI[295] - X
df_plt$HI[295] = df_plt$HI[295] - X

X = mean(df_plt$Median[296:304])
df_plt$Median[296:304] = df_plt$Median[296:304] - X
df_plt$LI[296:304] = df_plt$LI[296:304] - X
df_plt$HI[296:304] = df_plt$HI[296:304] - X

X = mean(df_plt$Median[305:308])
df_plt$Median[305:308] = df_plt$Median[305:308] - X
df_plt$LI[305:308] = df_plt$LI[305:308] - X
df_plt$HI[305:308] = df_plt$HI[305:308] - X


X = mean(df_plt$Median[309])
df_plt$Median[309] = df_plt$Median[309] - X
df_plt$LI[309] = df_plt$LI[309] - X
df_plt$HI[309] = df_plt$HI[309] - X

X = mean(df_plt$Median[310:318])
df_plt$Median[310:318] = df_plt$Median[310:318] - X
df_plt$LI[310:318] = df_plt$LI[310:318] - X
df_plt$HI[310:318] = df_plt$HI[310:318] - X

X = mean(df_plt$Median[319:322])
df_plt$Median[319:322] = df_plt$Median[319:322] - X
df_plt$LI[319:322] = df_plt$LI[319:322] - X
df_plt$HI[319:322] = df_plt$HI[319:322] - X


X = mean(df_plt$Median[323])
df_plt$Median[323] = df_plt$Median[323] - X
df_plt$LI[323] = df_plt$LI[323] - X
df_plt$HI[323] = df_plt$HI[323] - X

X = mean(df_plt$Median[324:332])
df_plt$Median[324:332] = df_plt$Median[324:332] - X
df_plt$LI[324:332] = df_plt$LI[324:332] - X
df_plt$HI[324:332] = df_plt$HI[324:332] - X

X = mean(df_plt$Median[333:336])
df_plt$Median[333:336] = df_plt$Median[333:336] - X
df_plt$LI[333:336] = df_plt$LI[333:336] - X
df_plt$HI[333:336] = df_plt$HI[333:336] - X


X = mean(df_plt$Median[337])
df_plt$Median[337] = df_plt$Median[337] - X
df_plt$LI[337] = df_plt$LI[337] - X
df_plt$HI[337] = df_plt$HI[337] - X

X = mean(df_plt$Median[338:346])
df_plt$Median[338:346] = df_plt$Median[338:346] - X
df_plt$LI[338:346] = df_plt$LI[338:346] - X
df_plt$HI[338:346] = df_plt$HI[338:346] - X

X = mean(df_plt$Median[347:350])
df_plt$Median[347:350] = df_plt$Median[347:350] - X
df_plt$LI[347:350] = df_plt$LI[347:350] - X
df_plt$HI[347:350] = df_plt$HI[347:350] - X


X = mean(df_plt$Median[351])
df_plt$Median[351] = df_plt$Median[351] - X
df_plt$LI[351] = df_plt$LI[351] - X
df_plt$HI[351] = df_plt$HI[351] - X

X = mean(df_plt$Median[352:360])
df_plt$Median[352:360] = df_plt$Median[352:360] - X
df_plt$LI[352:360] = df_plt$LI[352:360] - X
df_plt$HI[352:360] = df_plt$HI[352:360] - X

X = mean(df_plt$Median[361:364])
df_plt$Median[361:364] = df_plt$Median[361:364] - X
df_plt$LI[361:364] = df_plt$LI[361:364] - X
df_plt$HI[361:364] = df_plt$HI[361:364] - X


############### Merge
df_plt$Type = ifelse(str_detect(df_plt$Variable, "focal"), "Focal",
             ifelse(str_detect(df_plt$Variable, "target"), "Target",
             ifelse(str_detect(df_plt$Variable, "dyadic"), "Dyadic",
             ifelse(str_detect(df_plt$Variable, "error"), "Error",
                    NA))))

df_plt$Type = factor(df_plt$Type)
df_plt$Type = factor(df_plt$Type, levels=c("Focal", "Target", "Dyadic", "Error"))

df_plt$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", df_plt$Variable)
df_plt$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", df_plt$Variable)
df_plt$Variable = gsub("dyadic effects coeffs, ", "", df_plt$Variable)
df_plt$Variable = gsub("dyadic effects ", "", df_plt$Variable)
df_plt$Variable = gsub("focal effects ", "", df_plt$Variable)
df_plt$Variable = gsub("target effects ", "", df_plt$Variable)


df_plt$Outcome = ifelse(str_detect(df_plt$Variable, "Feeding"), "Feeding",
             ifelse(str_detect(df_plt$Variable, "Fighting"), "Fighting",
             ifelse(str_detect(df_plt$Variable, "Grooming"), "Grooming",
                ifelse(str_detect(df_plt$Variable, "Hunting"), "Hunting",
                    ifelse(str_detect(df_plt$Variable, "Smilling"), "Smilling",
                        ifelse(str_detect(df_plt$Variable, "Barking"), "Barking",
                            ifelse(str_detect(df_plt$Variable, "Killing"), "Killing",
                    NA)))))))

df_plt$Outcome = factor(df_plt$Outcome)
df_plt$Outcome = factor(df_plt$Outcome, levels=c("Fighting", "Feeding", "Grooming","Hunting","Smilling","Barking","Killing"))

df_plt$Variable = gsub("Fighting - ", "", df_plt$Variable)
df_plt$Variable = gsub("Feeding - ", "", df_plt$Variable)
df_plt$Variable = gsub("Grooming - ", "", df_plt$Variable)
df_plt$Variable = gsub("Hunting - ", "", df_plt$Variable)
df_plt$Variable = gsub("Smilling - ", "", df_plt$Variable)
df_plt$Variable = gsub("Barking - ", "", df_plt$Variable)
df_plt$Variable = gsub("Killing - ", "", df_plt$Variable)

df_plt$Variable = gsub(" - Fighting", "", df_plt$Variable)
df_plt$Variable = gsub(" - Feeding", "", df_plt$Variable)
df_plt$Variable = gsub(" - Grooming", "", df_plt$Variable)
df_plt$Variable = gsub(" - Hunting", "", df_plt$Variable)
df_plt$Variable = gsub(" - Smilling", "", df_plt$Variable)
df_plt$Variable = gsub(" - Barking", "", df_plt$Variable)
df_plt$Variable = gsub(" - Killing", "", df_plt$Variable)

df_plt$Variable = gsub("offset, ", "", df_plt$Variable)

df_plt$Variable = gsub("sd", "SD", df_plt$Variable)


df_plt$Block = ifelse(str_detect(df_plt$Variable, "Any"), "Intercept",
             ifelse(str_detect(df_plt$Variable, "Mottled"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Striped"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Spotted"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Male"), "Sex",
             ifelse(str_detect(df_plt$Variable, "Female"), "Sex",
                    NA))))))

df_plt_numpyro = df_plt
df_plt_numpyro$Method = "NumPyro"

####################### Now Stan  
df_plt = res_stan$summary

df_plt$TrueValues = c(
 sr_sigma[1:7],
 sr_1[1,],
 sr_2[1,],
 sr_3[1,],
 sr_4[1,],
 sr_5[1,],
 sr_6[1,],
 sr_7[1,],
 sr_sigma[8:14],
 sr_1[2,],
 sr_2[2,],
 sr_3[2,],
 sr_4[2,],
 sr_5[2,],
 sr_6[2,],
 sr_7[2,],
 dr_sigma,
 dr_effects[[1]],
 dr_effects[[2]],
 dr_effects[[3]],
 dr_effects[[4]],
 dr_effects[[5]],
 dr_effects[[6]],
 dr_effects[[7]],
 error_sigma,
 recip_to_long(sr_Rho),
 recip_to_long(dr_Rho),
 c(t(B_1a)) - mean(c(B_1a)),
 c(t(B_1b)) - mean(c(B_1b)), 
 c(t(B_1c)) - mean(c(B_1c)),
 c(t(B_2a)) - mean(c(B_2a)),
 c(t(B_2b)) - mean(c(B_2b)), 
 c(t(B_2c)) - mean(c(B_2c)),
 c(t(B_3a)) - mean(c(B_3a)),
 c(t(B_3b)) - mean(c(B_3b)), 
 c(t(B_3c)) - mean(c(B_3c)),
 c(t(B_4a)) - mean(c(B_4a)),
 c(t(B_4b)) - mean(c(B_4b)), 
 c(t(B_4c)) - mean(c(B_4c)),
 c(t(B_5a)) - mean(c(B_5a)),
 c(t(B_5b)) - mean(c(B_5b)), 
 c(t(B_5c)) - mean(c(B_5c)),
 c(t(B_6a)) - mean(c(B_6a)),
 c(t(B_6b)) - mean(c(B_6b)), 
 c(t(B_6c)) - mean(c(B_6c)),
 c(t(B_7a)) - mean(c(B_7a)),
 c(t(B_7b)) - mean(c(B_7b)), 
 c(t(B_7c)) - mean(c(B_7c))
  )

df_plt$Outcome2 = c(
 c("Feeding","Fighting","Grooming", "Hunting", "Smilling", "Barking", "Killing"),
 c("Feeding","Feeding","Feeding"),
 c("Fighting","Fighting","Fighting"),
 c("Grooming","Grooming","Grooming"),
 c("Hunting","Hunting","Hunting"),
 c("Smilling","Smilling","Smilling"),
 c("Barking","Barking","Barking"),
 c("Killing","Killing","Killing"),
 c("Feeding","Fighting","Grooming","Hunting", "Smilling", "Barking", "Killing"),
 c("Feeding","Feeding","Feeding"),
 c("Fighting","Fighting","Fighting"),
 c("Grooming","Grooming","Grooming"),
 c("Hunting","Hunting","Hunting"),
 c("Smilling","Smilling","Smilling"),
 c("Barking","Barking","Barking"),
 c("Killing","Killing","Killing"),
 c("Feeding","Fighting","Grooming","Hunting", "Smilling", "Barking", "Killing"),
 c("Feeding","Fighting","Grooming","Hunting", "Smilling", "Barking", "Killing"),
 c("Feeding","Feeding"),
 c("Fighting","Fighting"),
 c("Grooming","Grooming"),
 c("Hunting","Hunting"),
 c("Smilling","Smilling"),
 c("Barking","Barking"),
 c("Killing","Killing"),
 rep("Recip", 182),
 rep("Other", 98)
  )

df_plt$Variable2 = df_plt$Variable

colnames(df_plt) = c("Variable", "Median", "LI", "HI", "Mean", "SD", "P", "TrueValues","Outcome2", "Variable2")

df_plt$Median = as.numeric(df_plt$Median)
df_plt$LI = as.numeric(df_plt$LI)
df_plt$HI = as.numeric(df_plt$HI)


# Indentify block effect offsets
X = mean(df_plt$Median[267])
df_plt$Median[267] = df_plt$Median[267] - X
df_plt$LI[267] = df_plt$LI[267] - X
df_plt$HI[267] = df_plt$HI[267] - X

X = mean(df_plt$Median[268:276])
df_plt$Median[268:276] = df_plt$Median[268:276] - X
df_plt$LI[268:276] = df_plt$LI[268:276] - X
df_plt$HI[268:276] = df_plt$HI[268:276] - X

X = mean(df_plt$Median[277:280])
df_plt$Median[277:280] = df_plt$Median[277:280] - X
df_plt$LI[277:280] = df_plt$LI[277:280] - X
df_plt$HI[277:280] = df_plt$HI[277:280] - X


X = mean(df_plt$Median[281])
df_plt$Median[281] = df_plt$Median[281] - X
df_plt$LI[281] = df_plt$LI[281] - X
df_plt$HI[281] = df_plt$HI[281] - X

X = mean(df_plt$Median[282:290])
df_plt$Median[282:290] = df_plt$Median[282:290] - X
df_plt$LI[282:290] = df_plt$LI[282:290] - X
df_plt$HI[282:290] = df_plt$HI[282:290] - X

X = mean(df_plt$Median[291:294])
df_plt$Median[291:294] = df_plt$Median[291:294] - X
df_plt$LI[291:294] = df_plt$LI[291:294] - X
df_plt$HI[291:294] = df_plt$HI[291:294] - X


X = mean(df_plt$Median[295])
df_plt$Median[295] = df_plt$Median[295] - X
df_plt$LI[295] = df_plt$LI[295] - X
df_plt$HI[295] = df_plt$HI[295] - X

X = mean(df_plt$Median[296:304])
df_plt$Median[296:304] = df_plt$Median[296:304] - X
df_plt$LI[296:304] = df_plt$LI[296:304] - X
df_plt$HI[296:304] = df_plt$HI[296:304] - X

X = mean(df_plt$Median[305:308])
df_plt$Median[305:308] = df_plt$Median[305:308] - X
df_plt$LI[305:308] = df_plt$LI[305:308] - X
df_plt$HI[305:308] = df_plt$HI[305:308] - X


X = mean(df_plt$Median[309])
df_plt$Median[309] = df_plt$Median[309] - X
df_plt$LI[309] = df_plt$LI[309] - X
df_plt$HI[309] = df_plt$HI[309] - X

X = mean(df_plt$Median[310:318])
df_plt$Median[310:318] = df_plt$Median[310:318] - X
df_plt$LI[310:318] = df_plt$LI[310:318] - X
df_plt$HI[310:318] = df_plt$HI[310:318] - X

X = mean(df_plt$Median[319:322])
df_plt$Median[319:322] = df_plt$Median[319:322] - X
df_plt$LI[319:322] = df_plt$LI[319:322] - X
df_plt$HI[319:322] = df_plt$HI[319:322] - X


X = mean(df_plt$Median[323])
df_plt$Median[323] = df_plt$Median[323] - X
df_plt$LI[323] = df_plt$LI[323] - X
df_plt$HI[323] = df_plt$HI[323] - X

X = mean(df_plt$Median[324:332])
df_plt$Median[324:332] = df_plt$Median[324:332] - X
df_plt$LI[324:332] = df_plt$LI[324:332] - X
df_plt$HI[324:332] = df_plt$HI[324:332] - X

X = mean(df_plt$Median[333:336])
df_plt$Median[333:336] = df_plt$Median[333:336] - X
df_plt$LI[333:336] = df_plt$LI[333:336] - X
df_plt$HI[333:336] = df_plt$HI[333:336] - X


X = mean(df_plt$Median[337])
df_plt$Median[337] = df_plt$Median[337] - X
df_plt$LI[337] = df_plt$LI[337] - X
df_plt$HI[337] = df_plt$HI[337] - X

X = mean(df_plt$Median[338:346])
df_plt$Median[338:346] = df_plt$Median[338:346] - X
df_plt$LI[338:346] = df_plt$LI[338:346] - X
df_plt$HI[338:346] = df_plt$HI[338:346] - X

X = mean(df_plt$Median[347:350])
df_plt$Median[347:350] = df_plt$Median[347:350] - X
df_plt$LI[347:350] = df_plt$LI[347:350] - X
df_plt$HI[347:350] = df_plt$HI[347:350] - X


X = mean(df_plt$Median[351])
df_plt$Median[351] = df_plt$Median[351] - X
df_plt$LI[351] = df_plt$LI[351] - X
df_plt$HI[351] = df_plt$HI[351] - X

X = mean(df_plt$Median[352:360])
df_plt$Median[352:360] = df_plt$Median[352:360] - X
df_plt$LI[352:360] = df_plt$LI[352:360] - X
df_plt$HI[352:360] = df_plt$HI[352:360] - X

X = mean(df_plt$Median[361:364])
df_plt$Median[361:364] = df_plt$Median[361:364] - X
df_plt$LI[361:364] = df_plt$LI[361:364] - X
df_plt$HI[361:364] = df_plt$HI[361:364] - X


############### Merge
df_plt$Type = ifelse(str_detect(df_plt$Variable, "focal"), "Focal",
             ifelse(str_detect(df_plt$Variable, "target"), "Target",
             ifelse(str_detect(df_plt$Variable, "dyadic"), "Dyadic",
             ifelse(str_detect(df_plt$Variable, "error"), "Error",
                    NA))))

df_plt$Type = factor(df_plt$Type)
df_plt$Type = factor(df_plt$Type, levels=c("Focal", "Target", "Dyadic", "Error"))

df_plt$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", df_plt$Variable)
df_plt$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", df_plt$Variable)
df_plt$Variable = gsub("dyadic effects coeffs, ", "", df_plt$Variable)
df_plt$Variable = gsub("dyadic effects ", "", df_plt$Variable)
df_plt$Variable = gsub("focal effects ", "", df_plt$Variable)
df_plt$Variable = gsub("target effects ", "", df_plt$Variable)


df_plt$Outcome = ifelse(str_detect(df_plt$Variable, "Feeding"), "Feeding",
             ifelse(str_detect(df_plt$Variable, "Fighting"), "Fighting",
             ifelse(str_detect(df_plt$Variable, "Grooming"), "Grooming",
                ifelse(str_detect(df_plt$Variable, "Hunting"), "Hunting",
                    ifelse(str_detect(df_plt$Variable, "Smilling"), "Smilling",
                        ifelse(str_detect(df_plt$Variable, "Barking"), "Barking",
                            ifelse(str_detect(df_plt$Variable, "Killing"), "Killing",
                    NA)))))))

df_plt$Outcome = factor(df_plt$Outcome)
df_plt$Outcome = factor(df_plt$Outcome, levels=c("Fighting", "Feeding", "Grooming","Hunting","Smilling","Barking","Killing"))

df_plt$Variable = gsub("Fighting - ", "", df_plt$Variable)
df_plt$Variable = gsub("Feeding - ", "", df_plt$Variable)
df_plt$Variable = gsub("Grooming - ", "", df_plt$Variable)
df_plt$Variable = gsub("Hunting - ", "", df_plt$Variable)
df_plt$Variable = gsub("Smilling - ", "", df_plt$Variable)
df_plt$Variable = gsub("Barking - ", "", df_plt$Variable)
df_plt$Variable = gsub("Killing - ", "", df_plt$Variable)

df_plt$Variable = gsub(" - Fighting", "", df_plt$Variable)
df_plt$Variable = gsub(" - Feeding", "", df_plt$Variable)
df_plt$Variable = gsub(" - Grooming", "", df_plt$Variable)
df_plt$Variable = gsub(" - Hunting", "", df_plt$Variable)
df_plt$Variable = gsub(" - Smilling", "", df_plt$Variable)
df_plt$Variable = gsub(" - Barking", "", df_plt$Variable)
df_plt$Variable = gsub(" - Killing", "", df_plt$Variable)

df_plt$Variable = gsub("offset, ", "", df_plt$Variable)

df_plt$Variable = gsub("sd", "SD", df_plt$Variable)


df_plt$Block = ifelse(str_detect(df_plt$Variable, "Any"), "Intercept",
             ifelse(str_detect(df_plt$Variable, "Mottled"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Striped"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Spotted"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Male"), "Sex",
             ifelse(str_detect(df_plt$Variable, "Female"), "Sex",
                    NA))))))

df_plt_stan = df_plt
df_plt_stan$Method = "Stan"

########################## Plot 1
df_plt = rbind(df_plt_numpyro, df_plt_stan)
main_df = df_plt[which(df_plt$Outcome %in% c("Feeding","Fighting","Grooming","Hunting","Smilling","Barking","Killing")),]
main_df = df_plt[which(df_plt$Type %in% c("Focal","Target","Dyadic","Error")),]

main_df$Variable = factor(main_df$Variable)
main_df$Variable = factor(main_df$Variable, levels=rev(c("error SD", "SD", "Age", "Mass", "Strength", "Dominance", "Kinship")))

p = ggplot(main_df, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Method, color=Method)) + 
           geom_linerange(size = 1, position = position_dodge(width = 0.3)) + 
           geom_point(size = 2, position = position_dodge(width = 0.3)) +
           geom_point(size = 2, aes(x = Variable, y = TrueValues, group=Outcome), color="darkred", shape=18, position = position_dodge(width = 0.3)) +
           facet_grid(Type ~ Outcome, scales = "free_y", space = "free_y") + 
           geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") + 
           labs(y = "Regression parameters", x = "") + 
           theme(strip.text.x = element_text(size = 12, face = "bold"), 
                 strip.text.y = element_text(size = 12, face = "bold"), 
                 axis.text = element_text(size = 12), 
                 axis.title.y = element_text(size = 14, face = "bold"), 
                 axis.title.x = element_blank()) + 
           theme(strip.text.y = element_text(angle = 360)) + 
           coord_flip() + 
           theme(panel.spacing = unit(1,"lines")) + 
           theme(legend.position="bottom")
p

# We recover main effects


########################## Plot 2
block_df = df_plt[which(df_plt$Outcome2 == "Other" & df_plt$Block != "Intercept"),]

p = ggplot(block_df, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Method, color=Method)) + 
           geom_linerange(size = 1, position = position_dodge(width = 0.3)) + 
           geom_point(size = 2, position = position_dodge(width = 0.3)) +
           geom_point(size = 2, aes(x = Variable, y = TrueValues, group=Outcome), color="darkred", shape=18, position = position_dodge(width = 0.3)) +
           facet_grid(Block ~ Outcome, scales = "free", space = "free") + 
           geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") + 
           labs(y = "Regression parameters", x = "") + 
           theme(strip.text.x = element_text(size = 12, face = "bold"), 
                 strip.text.y = element_text(size = 12, face = "bold"), 
                 axis.text = element_text(size = 12), 
                 axis.title.y = element_text(size = 14, face = "bold"), 
                 axis.title.x = element_blank()) + 
           theme(strip.text.y = element_text(angle = 360)) + 
           coord_flip() + 
           theme(panel.spacing = unit(1,"lines")) + 
           theme(legend.position="bottom")
p


########################## Plot 3
recip_df = df_plt[which(df_plt$Outcome2 == "Recip"),]

recip_df$Type = ifelse(str_detect(recip_df$Variable2, "Generalized"), "Generalized",
                ifelse(str_detect(recip_df$Variable2, "Dyadic"), "Dyadic",
                    NA))

recip_df$Variable2 = gsub("Dyadic reciprocity - ", "", recip_df$Variable2)
recip_df$Variable2 = gsub("Generalized reciprocity - ", "", recip_df$Variable2)

p1 = ggplot(recip_df[which(recip_df$Type == "Generalized"),], aes(x = Variable2, y = Median, ymin = LI, ymax = HI, group=Method, color=Method)) + 
           geom_linerange(size = 1, position = position_dodge(width = 0.3)) + 
           geom_point(size = 2, position = position_dodge(width = 0.3)) +
           geom_point(size = 2, aes(x = Variable2, y = TrueValues, group=Outcome), color="darkred", shape=18, position = position_dodge(width = 0.3)) +
           facet_grid(. ~ Type, scales = "free", space = "free") + 
           geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") + 
           labs(y = "Regression parameters", x = "") + 
           theme(strip.text.x = element_text(size = 12, face = "bold"), 
                 strip.text.y = element_text(size = 12, face = "bold"), 
                 axis.text = element_blank(), 
                 axis.title.y = element_text(size = 14, face = "bold"), 
                 axis.title.x = element_blank()) + 
           theme(strip.text.y = element_text(angle = 90)) + 
           theme(panel.spacing = unit(1,"lines")) + 
           theme(legend.position="bottom")
p1



p2 = ggplot(recip_df[which(recip_df$Type == "Dyadic"),], aes(x = Variable2, y = Median, ymin = LI, ymax = HI, group=Method, color=Method)) + 
           geom_linerange(size = 1, position = position_dodge(width = 0.3)) + 
           geom_point(size = 2, position = position_dodge(width = 0.3)) +
           geom_point(size = 2, aes(x = Variable2, y = TrueValues, group=Outcome), color="darkred", shape=18, position = position_dodge(width = 0.3)) +
           facet_grid(. ~ Type, scales = "free", space = "free") + 
           geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") + 
           labs(y = "Regression parameters", x = "") + 
           theme(strip.text.x = element_text(size = 12, face = "bold"), 
                 strip.text.y = element_text(size = 12, face = "bold"), 
                 axis.text = element_blank(), 
                 axis.title.y = element_text(size = 14, face = "bold"), 
                 axis.title.x = element_blank()) + 
           theme(strip.text.y = element_text(angle = 90)) + 
           theme(panel.spacing = unit(1,"lines")) + 
           theme(legend.position="bottom")
p2

# We recover dyadic reciprocity, within and between layers, too 


###############################################
# Dont forget to check rhat and ess

# Stan
 fit_stan$fit$summary("dyad_effects")
 fit_stan$fit$summary("focal_effects")
 fit_stan$fit$summary("target_effects")
 fit_stan$fit$summary("block_effects")

# JAX
 samples = fit_numpyro$fit$get_samples()
 jax_summary(samples$dyad_effects)
 jax_summary(samples$focal_effects)
 jax_summary(samples$target_effects)
 jax_summary(samples$block_effects)



