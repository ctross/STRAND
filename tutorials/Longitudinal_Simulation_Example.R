##################################################### Simulate longitudinal networks
###################################### Load packages
# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors
 color_set = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

 library(stringr)
 library(ggplot2)
 library(psych)
 library(rethinking)
 library(Matrix)
 library(igraph)
 library(STRAND)

###################################### Make simulated longitudinal data
 set.seed(666)

# Base settings
 N_id = 125       # Individuals in network
 N_timesteps = 5  # Network layers
 labels = paste("ind", 1:N_id)

# Individual covariates - Fixed
 Mass = rbern(N_id, 0.4)
 Age = rnorm(N_id, 0, 1)
 Strength = rnorm(N_id, 0, 1)

# Dyadic covariates - Fixed
 Distance = standardize(as.matrix(rlkjcorr(1, N_id, eta=1.5)))
 Dominance = as.matrix(ceiling(rlkjcorr(1, N_id, eta=1.5) - 0.1))
 colnames(Distance) = rownames(Distance) = labels
 colnames(Dominance) = rownames(Dominance) = labels

# Individual covariates - Time varying
 Food = Size = NULL

 SizeBase = rnorm(N_id, 0, 1)
 FoodBase = rnorm(N_id, 0, 1)
 
 for(t in 1:N_timesteps){
 Size[[t]] =  SizeBase + 0.1*rnorm(N_id, 0, 1)
 Food[[t]] =  FoodBase + 0.2*rnorm(N_id, 0, 1)
 }

# Dyadic covariates - Time varying
 Agression = Friendship = NULL

 Friendship_prob = simulate_srm_network(N_id=N_id, B=-3)$tie_strength
 Agression_prob = simulate_srm_network(N_id=N_id, B=-1)$tie_strength
 scrap1 = scrap2 = Friendship_prob

 for(t in 1:N_timesteps){
    for(i in 1:N_id){
        for(j in 1:N_id){
             scrap1[i,j] = rbinom(1, size=1, prob=Friendship_prob[i,j])
             scrap2[i,j] = rbinom(1, size=1, prob=Agression_prob[i,j])
        }
    }
   Friendship[[t]] =  scrap1
   Agression[[t]] =  scrap2
   colnames(Friendship[[t]]) = rownames(Friendship[[t]]) = labels
   colnames(Agression[[t]]) = rownames(Agression[[t]]) = labels
 }

# Block structure - Fixed
 groups_1 = rep("Any",N_id) 
 groups_2 = sample( c("Male", "Female"), size=N_id, replace=TRUE, prob=c(0.5, 0.5))

# Block structure - Varying
 groups_3 = NULL

 groups_3[[1]] = sample( c("Mottled", "Striped", "Spotted"), size=N_id, replace=TRUE, prob=c(0.2, 0.25, 0.55))

 for(t in 2:N_timesteps){
  rando = rbinom(N_id, size = 1, prob = 0.25)  
  temp = rep(NA, N_id)
  types = groups_3[[t-1]]

  for(i in 1:length(temp)){
    if(rando[i] == 1){
      temp[i] = sample( c("Mottled", "Striped", "Spotted"), size=1, replace=TRUE, prob=c(0.2, 0.25, 0.55))   
    } else{
      temp[i] = types[i]   
    }
   }
  groups_3[[t]] = temp
 }

####################################### Organize data lists
 block_predictors_f = block_predictors = dyadic_predictors = individual_predictors = NULL
 dyadic_preds = array(NA, c(N_id, N_id, 4))
 labels = paste("Ind", 1:N_id)

 for(t in 1:N_timesteps){
   individual_predictors[[t]] = data.frame(Age = Age, Food = Food[[t]], Mass = Mass, Size = Size[[t]], Strength = Strength)
   rownames(individual_predictors[[t]]) = labels

   dyadic_preds[,,1] = Agression[[t]]
   dyadic_preds[,,2] = Distance
   dyadic_preds[,,3] = Dominance
   dyadic_preds[,,4] = Friendship[[t]]

   dyadic_predictors[[t]] = dyadic_preds

   block_predictors[[t]] = data.frame(Intercept = as.numeric(factor(groups_1)), Coloration = as.numeric(factor(groups_3[[t]])), Sex = as.numeric(factor(groups_2)))
   block_predictors_f[[t]] = data.frame(Intercept = factor(groups_1), Coloration = factor(groups_3[[t]]), Sex = factor(groups_2))
   rownames(block_predictors_f[[t]]) = labels
 }

###################################### Make simulation parameters
# Set simple random effects sizes
 sr_mu = rep(0, N_timesteps*2)                                          # Average effect size, should be zero
 sr_sigma = c(2.2, 1.7, 1.3, 0.5, 1.7,                                  # Variation in random effects. First 5 are sender effects, one per layer. Last 5 are receiver effects.
              1.1, 0.3, 1.6, 0.9, 1.5)                                  #
 dr_mu = rep(0, N_timesteps)                                            # Average effect size, should be zero
 dr_sigma = c(0.9, 3.1, 2.5, 0.4, 1.3)                                  # Variation in dyadic random effects.

# Dyadic reciprocity matrix
  D = matrix(NA, nrow=(N_timesteps*2), ncol=(N_timesteps*2))
  D[1,] = c(c(1, 0.62187, 0.541227, 0.42,  0.31), 0.7*c(0.62187, 0.561227, 0.41,  0.31, 0.155443))  
    
  for(k in 1:(N_timesteps-1)){
    K =  N_timesteps - k
    for(m in 1:K){
     D[m, m+k] = D[1, k+1]
     D[m, N_timesteps + m + k] = D[1, N_timesteps + k + 1]
    }
   }
     
  for(m in 1:N_timesteps){
     D[m, m+N_timesteps] = D[1, N_timesteps + 1]
    }

  for(m in 1:(N_timesteps-1)){
   for(n in (m+1):N_timesteps){
     D[m+N_timesteps, n+N_timesteps] = D[m, n] 
    }}

  for(m in 1:(N_timesteps-1)){
   for(n in (m+1):N_timesteps){
     D[n, m+N_timesteps] = D[m, n+N_timesteps]
    }}

  for(i in 1:((2*N_timesteps)-1)){
   for(j in (i+1):(N_timesteps*2)){
     D[j, i] = D[i, j]
    }}

  for(i in 1:((2*N_timesteps))){
     D[i, i] = 1
    }

    D2 = as.matrix(round(D,3))

    D3 = nearPD(D2, corr=TRUE)

    dr_Rho = as.matrix(D3$mat)  


# Generalized reciprocity matrix
  G = matrix(NA, nrow=(N_timesteps*2), ncol=(N_timesteps*2))
  
    Vals = c(0.4, 0.35, 0.28, 0.15)
    Vals2 = c(0.19, 0.13, 0.07, 0.02)
    Vals3 = c(0.61, 0.61, 0.61, 0.61, 0.61)
    Vals4 = c(0.5, 0.4, 0.3, 0.25)
    Vals5 = c(0.2, 0.15, 0.07, 0.01)

    target = 0
    for(m in 1:(N_timesteps-1)){
    for(n in (m+1):N_timesteps){
     target = target + 1
     G[m, n] = Vals[target]
     G[m+N_timesteps, n+N_timesteps] = Vals4[target]

     G[n, m+N_timesteps] = Vals2[target]
     G[m, n+N_timesteps] = Vals5[target]
    }}

  for(m in 1:N_timesteps){
     G[m, m+N_timesteps] = Vals3[m]
    }

  for(k in 1:(N_timesteps-1)){
    for(m in 1:(N_timesteps - k)){
     G[m, m+k] = G[1, k+1]
     G[m, N_timesteps + m + k] = G[1, N_timesteps + k + 1]
     G[N_timesteps + m, N_timesteps + m+k] = G[N_timesteps + 1, N_timesteps + k+1]
     G[m + k, N_timesteps + m] = G[k + 1, N_timesteps + 1]
    }
   }

  for(m in 1:N_timesteps){
     G[m, m+N_timesteps] = G[1, N_timesteps + 1]
    }

  for(i in 1:((2*N_timesteps)-1)){
   for(j in (i+1):(N_timesteps*2)){
     G[j, i] = G[i, j]
    }}

  for(i in 1:((2*N_timesteps))){
     G[i, i] = 1
    }

    G2 = as.matrix(round(G,3))

    G3 = nearPD(G2, corr=TRUE)

    sr_Rho = as.matrix(G3$mat)  


# Covariate effects                                  
 sr_1 = matrix(NA, nrow=2, ncol=5)                        # Time 1 
 sr_1[1,] = c(1.1, 1.3, -0.1, 0.4, 0.7)                   # Effect of Age, Food, Mass, Size, and Strength on out-degree
 sr_1[2,] = c(0.7, 1.4, -0.1, 0.7, -0.31)                 # Effect of Age, Food, Mass, Size, and Strength on in-degree

 sr_2 = sr_1 + matrix(rnorm(10, 0, 0.22), nrow=2, ncol=5) # Time 2 - After that, let them Brownian walk
 sr_3 = sr_2 + matrix(rnorm(10, 0, 0.22), nrow=2, ncol=5) # Time 3
 sr_4 = sr_3 + matrix(rnorm(10, 0, 0.22), nrow=2, ncol=5) # Time 4
 sr_5 = sr_4 + matrix(rnorm(10, 0, 0.22), nrow=2, ncol=5) # Time 5
                
 sr_effects = list(sr_1, sr_2, sr_3, sr_4, sr_5)          # Organize into list

 dr_1 = c(-1.4, -0.7, 2.5, -1.1)                          # Layer 1 Agression, Distance, Dominance, and Friendship     
 dr_2 = dr_1 + rnorm(4, 0, 0.22)                          # Time 2 - After that, let them Brownian walk
 dr_3 = dr_2 + rnorm(4, 0, 0.22)                          # Time 3
 dr_4 = dr_3 + rnorm(4, 0, 0.22)                          # Time 4
 dr_5 = dr_4 + rnorm(4, 0, 0.22)                          # Time 5

 dr_effects = list(dr_1, dr_2, dr_3, dr_4, dr_5)         # Organize into list


# Intercept in time-step
 B_1a = matrix(-12.2, nrow=1, ncol=1)
 B_2a = B_1a + rnorm(1, 0, 0.2)
 B_3a = B_2a + rnorm(1, 0, 0.2)
 B_4a = B_3a + rnorm(1, 0, 0.2)
 B_5a = B_4a + rnorm(1, 0, 0.2)

 # Offsets for Coloration
 B_1b = matrix(rnorm(9, 0, 2), nrow=3, ncol=3)
 diag(B_1b) = diag(B_1b) + 4.5

 B_2b = B_1b + matrix(rnorm(9, 0, 1.2), nrow=3, ncol=3)
 B_3b = B_2b + matrix(rnorm(9, 0, 1.2), nrow=3, ncol=3)
 B_4b = B_3b + matrix(rnorm(9, 0, 1.2), nrow=3, ncol=3)
 B_5b = B_4b + matrix(rnorm(9, 0, 1.2), nrow=3, ncol=3)

# Offset for Sex
 B_1c = matrix(rnorm(4, 0, 1), nrow=2, ncol=2)
 diag(B_1c) = diag(B_1c) + 2.5

 B_2c = B_1c + matrix(rnorm(4, 0, 1.2), nrow=2, ncol=2)
 B_3c = B_2c + matrix(rnorm(4, 0, 1.2), nrow=2, ncol=2)
 B_4c = B_3c + matrix(rnorm(4, 0, 1.2), nrow=2, ncol=2)
 B_5c = B_4c + matrix(rnorm(4, 0, 1.2), nrow=2, ncol=2)

# Merge into lists
 B1 = list(B_1a, B_1b, B_1c)
 B2 = list(B_2a, B_2b, B_2c)
 B3 = list(B_3a, B_3b, B_3c)
 B4 = list(B_4a, B_4b, B_4c)
 B5 = list(B_5a, B_5b, B_5c)

 B = list(B1, B2, B3, B4, B5)


#################################################### Simulate network
A = simulate_longitudinal_network(
  N_id = N_id,            
  N_timesteps = N_timesteps,                   
  B = B,                       
  V = 3,       
  groups = block_predictors,                    
  sr_mu = sr_mu,          
  sr_sigma = sr_sigma,                        
  sr_Rho = sr_Rho,                    
  dr_mu = dr_mu,                            
  dr_sigma = dr_sigma,                        
  dr_Rho = dr_Rho,                         
  outcome_mode="bernoulli",   
  link_mode="logit",              
  individual_predictors = individual_predictors,  
  dyadic_predictors = dyadic_predictors,    
  individual_effects = sr_effects,  
  dyadic_effects = dr_effects
 )




#################################################### Create the STRAND data object
 long_dat = NULL
 for(t in 1:5){
  dyadic_preds = dyadic_predictors[[t]]

  dyadic_predictors_temp = list(Agression = dyadic_preds[,,1],
                                Distance = dyadic_preds[,,2], 
                                Dominance = dyadic_preds[,,3],
                                Friendship = dyadic_preds[,,4]
                            )

  temp_outcome = A$network[t,,]
  rownames(temp_outcome) = colnames(temp_outcome) = labels

  rownames(dyadic_predictors_temp$Agression) = colnames(dyadic_predictors_temp$Agression) = labels
  rownames(dyadic_predictors_temp$Distance) = colnames(dyadic_predictors_temp$Distance) = labels
  rownames(dyadic_predictors_temp$Dominance) = colnames(dyadic_predictors_temp$Dominance) = labels
  rownames(dyadic_predictors_temp$Friendship) = colnames(dyadic_predictors_temp$Friendship) = labels

  long_dat[[t]] = make_strand_data(outcome = list(Tolerance=temp_outcome),
                       block_covariates = block_predictors_f[[t]], 
                       individual_covariates = individual_predictors[[t]], 
                       dyadic_covariates = dyadic_predictors_temp,
                       longitudinal = TRUE,
                       outcome_mode = "bernoulli",
                       link_mode = "logit")
 }

 names(long_dat) = paste("Time", c(1:5))

# Model
fit_4 = fit_longitudinal_model(long_data=long_dat,
                               block_regression = ~ Coloration + Sex,
                               focal_regression = ~ Age + Food + Mass + Size + Strength,
                               target_regression = ~ Age + Food + Mass + Size + Strength,
                               dyad_regression = ~ Agression + Distance + Dominance + Friendship,
                               coefficient_mode="varying",
                               random_effects_mode="fixed",
                               mode="mcmc",
                               stan_mcmc_parameters = list(seed = 1, chains = 1, parallel_chains = 1, refresh = 1, iter_warmup = 500,
                                                           iter_sampling = 500, max_treedepth = 12, adapt_delta = NULL),
                               priors=NULL
                               )

res_4 = summarize_longitudinal_bsrm_results(fit_4)


########################################################### Coeff recovery
out = longitudinal_plot(fit_4, type="coefficient", 
    parameter_set = list(
    focal="Mass", focal="Strength", focal="Age", focal="Food", focal="Size",
    target="Mass", target="Strength", target="Age", target="Food", target="Size",
    dyadic="Agression", dyadic="Distance", dyadic="Dominance", dyadic="Friendship"),
    normalized=FALSE,
    plot = FALSE,
    export_as_table=TRUE
    )

 df_plt = res_4$summary
 df_plt$Block = ifelse(str_detect(df_plt$Variable, "effects sd"), "SD", "Other")

 df_plt$type_set = ifelse(str_detect(df_plt$Variable, "dyadic"), "Dyadic", 
                  ifelse(str_detect(df_plt$Variable, "focal"),"Focal","Target"))
                   

 df_plt$Wave = ifelse(str_detect(df_plt$Variable, "Time 1"), 1,
             ifelse(str_detect(df_plt$Variable, "Time 2"), 2,
             ifelse(str_detect(df_plt$Variable, "Time 3"), 3,
             ifelse(str_detect(df_plt$Variable, "Time 4"), 4,
             ifelse(str_detect(df_plt$Variable, "Time 5"), 5,
                    NA)))))

 df_plt$Wave2 = ifelse(str_detect(df_plt$Variable, "Time 1"), "Time 1",
             ifelse(str_detect(df_plt$Variable, "Time 2"), "Time 2",
             ifelse(str_detect(df_plt$Variable, "Time 3"), "Time 3",
             ifelse(str_detect(df_plt$Variable, "Time 4"), "Time 4",
             ifelse(str_detect(df_plt$Variable, "Time 5"), "Time 5",
                    NA)))))

 df_plt_2 = df_plt[which(df_plt$Block=="SD"),]
 df_plt_2$True = c(sr_sigma, dr_sigma)
 df_plt_2$extra_short_names = paste0(df_plt_2$type_set, " - SD")
 

 s_temp = r_temp = NULL
 for(i in 1:5){
  s_temp[[i]] = sr_effects[[i]][1,]
  r_temp[[i]] = sr_effects[[i]][2,]
 }

 true_set = c(unlist(dr_effects), unlist(s_temp), unlist(r_temp))

 out$True = true_set


 df_plt_3 = df_plt_2[,c("extra_short_names", "Median", "HPDI:0.05", "HPDI:0.95", "Block", "Wave2", "True")]
 out2 = out[,c("extra_short_names", "Median", "HPDI:0.05", "HPDI:0.95", "type_set", "time_point", "True")]


 colnames(out2) = colnames(df_plt_3) = c("extra_short_names", "Median", "l", "h","type","time_point","True")

 out3 = rbind(out2, df_plt_3)

 out3$type = factor(out3$type)
 out3$type = factor(out3$type, levels=c("Focal", "Target", "Dyadic", "SD"))

 pal = plvs_vltra("mystic_mausoleum",elements=c(10,11,12,13,14,1,3,5,7,9,2,4,6,8))
 pal = c(pal[1:4], "grey10", pal[5:7], "grey25", pal[8:12], "grey40", pal[13:14])

 p0 = ggplot2::ggplot(out3, ggplot2::aes(x=extra_short_names, y=as.numeric(Median), ymin=as.numeric(l), ymax=as.numeric(h), group=extra_short_names, color=extra_short_names))+ 
     ggplot2::geom_linerange(size=1, position = ggplot2::position_dodge(width = 0.3)) + ggplot2::facet_grid(type~time_point, scales="free", space="free") +
     ggplot2::geom_point(size=2, position = ggplot2::position_dodge(width = 0.3))+ ggplot2::geom_point(aes(x=extra_short_names, y=as.numeric(True)),size=2, position = ggplot2::position_dodge(width = 0.3),color="deeppink3",shape=18)+
     ggplot2::geom_hline(ggplot2::aes(yintercept=0),color="black",linetype="dashed")+
     ggplot2::labs(y="Effect size", x="Variable") + 
     ggplot2::theme(strip.text.x = ggplot2::element_text(size=12,face="bold"), 
      strip.text.y = ggplot2::element_text(size=12,face="bold"),
      axis.text = ggplot2::element_text(size=12),
      axis.title = ggplot2::element_text(size=14, face="bold"))+
     ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
     ggplot2::coord_flip() + 
     ggplot2::theme(panel.spacing = grid::unit(1, "lines")) + ggplot2::scale_color_manual(values = pal) + 
     ggplot2::theme(legend.position="none") + ggplot2::theme(legend.title = ggplot2::element_blank())
p0
     ggsave("Covariate_Recovery.pdf", p0, height=5, width=12)

########################################################### Dyadic recovery
 out = longitudinal_plot(fit_4, type="dyadic", 
    palette=pal,
    normalized=FALSE,
    plot = FALSE,
    export_as_table=TRUE
    )

  out$True = dr_Rho[1,2:10]

 pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,2))

 p1 = ggplot2::ggplot(out, ggplot2::aes(x=Time2, y=as.numeric(rs_m), ymin=as.numeric(l), ymax=as.numeric(h), group=rs_type, color=rs_type))+ 
     ggplot2::geom_linerange(size=1, position = ggplot2::position_dodge(width = 0.3)) +
     ggplot2::geom_point(size=2, position = ggplot2::position_dodge(width = 0.3))+ ggplot2::geom_point(aes(x=Time2, y=as.numeric(True)),size=3, position = ggplot2::position_dodge(width = 0.3),color="deeppink3",shape=18)+
     ggplot2::geom_hline(ggplot2::aes(yintercept=0),color="black",linetype="dashed")+
     ggplot2::labs(y="Effect size", x="Time step") + 
     ggplot2::theme(strip.text.x = ggplot2::element_text(size=12,face="bold"), 
      strip.text.y = ggplot2::element_text(size=12,face="bold"),
      axis.text = ggplot2::element_text(size=12),
      axis.title = ggplot2::element_text(size=14, face="bold"))+
     ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
    # ggplot2::coord_flip() + 
     ggplot2::theme(panel.spacing = grid::unit(1, "lines")) + ggplot2::scale_color_manual(values = pal) + 
     ggplot2::theme(legend.position="bottom") + ggplot2::theme(legend.title = ggplot2::element_blank())
p1
 ggsave("Dyadic_Recovery.pdf", p1, height=6.5, width=6.5)

########################################################### Generalized recovery
out = longitudinal_plot(fit_4, type="generalized", 
    palette=pal,
    normalized=FALSE,
    plot = FALSE,
    export_as_table=TRUE
    )

 out$True = c(sr_Rho[1,2:5], sr_Rho[6,7:10], sr_Rho[1,7:10], sr_Rho[5,9:6], sr_Rho[1,6])

 pal = plvs_vltra("mystic_mausoleum", elements=c(1,2,9,3,4))

 p2 = ggplot2::ggplot(out, ggplot2::aes(x=Time, y=as.numeric(rs_m), ymin=as.numeric(l), ymax=as.numeric(h), group=Type, color=Type))+ 
     ggplot2::geom_linerange(size=1, position = ggplot2::position_dodge(width = 0.3)) +
     ggplot2::geom_point(size=2, position = ggplot2::position_dodge(width = 0.3))+ ggplot2::geom_point(aes(x=Time, y=as.numeric(True)),size=3, position = ggplot2::position_dodge(width = 0.3),color="deeppink3",shape=18)+
     ggplot2::geom_hline(ggplot2::aes(yintercept=0),color="black",linetype="dashed")+
     ggplot2::labs(y="Effect size", x="Time step") + 
     ggplot2::theme(strip.text.x = ggplot2::element_text(size=12,face="bold"), 
      strip.text.y = ggplot2::element_text(size=12,face="bold"),
      axis.text = ggplot2::element_text(size=12),
      axis.title = ggplot2::element_text(size=14, face="bold"))+
     ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
    # ggplot2::coord_flip() + 
     ggplot2::theme(panel.spacing = grid::unit(1, "lines")) + ggplot2::scale_color_manual(values = pal) + 
     ggplot2::theme(legend.position="bottom") + ggplot2::theme(legend.title = ggplot2::element_blank())  
p2
 ggsave("Generalized_Recovery.pdf", p2, height=6.5, width=6.5)

######################################################### Visualize results
df_plt = res_4$summary
df_plt$Block = ifelse(str_detect(df_plt$Variable, "Any"), "Intercept",
             ifelse(str_detect(df_plt$Variable, "Mottled"), "Coloration",
             ifelse(str_detect(df_plt$Variable, "Striped"), "Coloration",
             ifelse(str_detect(df_plt$Variable, "Spotted"), "Coloration",
             ifelse(str_detect(df_plt$Variable, "Male"), "Sex",
             ifelse(str_detect(df_plt$Variable, "Female"), "Sex",
                    NA))))))

df_plt$Wave = ifelse(str_detect(df_plt$Variable, "Time 1"), 1,
             ifelse(str_detect(df_plt$Variable, "Time 2"), 2,
             ifelse(str_detect(df_plt$Variable, "Time 3"), 3,
             ifelse(str_detect(df_plt$Variable, "Time 4"), 4,
             ifelse(str_detect(df_plt$Variable, "Time 5"), 5,
                    NA)))))

df_plt$Wave2 = ifelse(str_detect(df_plt$Variable, "Time 1"), "Time 1",
             ifelse(str_detect(df_plt$Variable, "Time 2"), "Time 2",
             ifelse(str_detect(df_plt$Variable, "Time 3"), "Time 3",
             ifelse(str_detect(df_plt$Variable, "Time 4"), "Time 4",
             ifelse(str_detect(df_plt$Variable, "Time 5"), "Time 5",
                    NA)))))

################################################################# Coloration
make_blocks = function(df,B,t,Par,par){
  pat = df[which(df$Block==Par & df$Wave == t),]
  colnames(pat) = c("Variable", "Median", "L", "H", "Mean", "SD", "Block", "Wave", "Wave2")
  pat$Variable = gsub(paste0("offset, Time ", t, " - "), "", pat$Variable)

  X = mean(as.numeric(pat$Mean))
  Xt = mean(B[[t]][[par]])

  pat$True = as.numeric(c(t(B[[t]][[par]]))) - Xt
  pat$M = as.numeric(pat$Median) - X
  pat$L = as.numeric(pat$L) - X
  pat$H = as.numeric(pat$H) - X

  return(pat)               
}

res_blocks_1 = NULL
for(t in 1:5)
res_blocks_1[[t]] = make_blocks(df=df_plt, B=B, t=t, Par="Coloration", par=2)

res_blocks_1 = do.call(rbind, res_blocks_1)
res_blocks_1$Var = "Coloration"

res_blocks_2 = NULL
for(t in 1:5)
res_blocks_2[[t]] = make_blocks(df=df_plt, B=B, t=t, Par="Sex", par=3)

res_blocks_2 = do.call(rbind, res_blocks_2)
res_blocks_2$Var = "Sex"


res_blocks = rbind(res_blocks_1, res_blocks_2)

pal = plvs_vltra("mystic_mausoleum",elements=c(1:12,14))

p3 = ggplot(res_blocks, aes(x = Variable, y = M, ymin = L, ymax = H, color=Variable)) + 
           geom_linerange(size = 1) + 
           geom_point(size = 2) +
           geom_point(size = 3, aes(x = Variable, y = True), color="deeppink3", shape=18) +
           facet_grid(Var ~ Wave2, scales = "free", space = "free") + 
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
           theme(legend.position="bottom")  + scale_color_manual(values = pal) + theme(legend.position="none")
p3

ggsave("Block_Recovery.pdf", p3, height=4, width=12)



