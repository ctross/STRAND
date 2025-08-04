###############################################################
#
#   Undirected layer in multiplex analyses with data simulation 
#
########################################

# Clear working space
rm(list = ls())
set.seed(1111)

# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors
 colors = plvs_vltra("dust_storm", rev=FALSE, elements=NULL, show=FALSE)

library(rethinking)
library(STRAND)
library(stringr)
library(ggplot2)
library(psych)


# Make data
 N_id = 100      # Individuals in network
 N_layers = 3    # Network layers

# Covariates
 Kinship = STRAND::standardize(rlkjcorr( 1 , N_id , eta=1.5 ))      # Dyadic covariate, should be undirected since it predicts undirected ties
 Dominance = ceiling(rlkjcorr( 1 , N_id , eta=1.5 ) - 0.1)          # Dyadic covariate, should be undirected since it predicts undirected ties
 Mass = rbern(N_id, 0.4)
 Age = rnorm(N_id, 0, 1)
 Strength = rnorm(N_id, 0, 1)

# Organize into list
 dyadic_preds = array(NA,c(N_id,N_id,2))

 dyadic_preds[,,1] = Kinship
 dyadic_preds[,,2] = Dominance

# Set effect sizes
sr_mu = rep(0, N_layers*2)                                                              # Average random effect size, should be zero
sr_sigma = c(2.2, 0.7, 2.8, 1.3, 1.7, 2.8)                                              # Variation in random effects. First 3 are sender effects, one per layer. Last 3 are receiver effects.
sr_Rho = structure( c(1, -0.0950240087167285, -0.0877205103265762, -0.0124498813127008, # Generalized reciprocity matrix. # Force layer 3 to have near symmetric sender and receiver random effects with cor=0.98
    0.213083065621364, -0.000940306650023561, -0.0950240087167285, 
    1, -0.266151640364901, -0.00700343702170813, -0.266687239331943, 
    -0.210308163966482, -0.0877205103265762, -0.266151640364901, 
    1, -0.0873769988222759, -0.266064125271499, 0.98, -0.0124498813127008, 
    -0.00700343702170813, -0.0873769988222759, 1, -0.322752782831108, 
    -0.00126130028222611, 0.213083065621364, -0.266687239331943, 
    -0.266064125271499, -0.322752782831108, 1, -0.210389942229842, 
    -0.000940306650023561, -0.210308163966482, 0.98, -0.00126130028222611, 
    -0.210389942229842, 1), dim = c(6L, 6L))          
chol(sr_Rho) 

dr_mu = rep(0, N_layers)                        # Average random effect size, should be zero
dr_sigma = c(1.9, 1.1, 2.9)                     # Variation in dyadic random effects.

# Build dyadic matrix. Force layer 3 to be near symmetric/undirected by setting corr=0.98
dr_Rho = structure(c(1, -0.104394026927657, 0.251667701975812, 0.267417782629568, 
-0.034981842119401, 0.251667701975812, -0.104394026927657, 1, 
-0.0855352102622197, -0.034981842119401, -0.0675187589362978, 
-0.0855352102622197, 0.251667701975812, -0.0855352102622197, 
1, 0.251667701975812, -0.0855352102622197, 0.98, 0.267417782629568, 
-0.034981842119401, 0.251667701975812, 1, -0.104394026927657, 
0.251667701975812, -0.034981842119401, -0.0675187589362978, -0.0855352102622197, 
-0.104394026927657, 1, -0.0855352102622197, 0.251667701975812, 
-0.0855352102622197, 0.98, 0.251667701975812, -0.0855352102622197, 
1), dim = c(6L, 6L))
chol(dr_Rho)                                    # Check if positive definite

# Covariate effects
sr_1 = matrix(NA, nrow=2, ncol=3)               # Layer 1 - Directed, so free to differ across in- and out-
sr_1[1,] = c(-0.5, 1.0, -0.7)                    # Effect of Mass, Age, and Strength on out-degree
sr_1[2,] = c(0.7, -1.1, -1)                      # Effect of Mass, Age, and Strength on in-degree

sr_2 = matrix(NA, nrow=2, ncol=3)               # Layer 2 - Directed, so free to differ across in- and out-
sr_2[1,] = c(-0.1, -1.0, 0.7)                    # Effect of Mass, Age, and Strength on out-degree
sr_2[2,] = c(-0.7, -0.6, -0.01)                  # Effect of Mass, Age, and Strength on in-degree

sr_3 = matrix(NA, nrow=2, ncol=3)               # Layer 3 - Undirected, so not free to differ across in- and out-
sr_3[1,] = c(1.1, 0.3, -0.95)                    # Effect of Mass, Age, and Strength on out-degree
sr_3[2,] = c(1.1, 0.3, -0.95)                    # Effect of Mass, Age, and Strength on in-degree

sr_effects = list(sr_1, sr_2, sr_3)             # Organize into list

dr_effects = list(c(0.6, 0.3),                  # Layer 1 effect of Kinship and Dominant        
                  c(-0.2, -0.7),                # Layer 2 effect of Kinship and Dominant 
                  c(1.62, 0.7))                 # Layer 3 effect of Kinship and Dominant 

# Block structure
group_probs_block_size = c(0.25, c(0.25, 0.75)*(1-0.25))
groups_1 = rep("Any",N_id) 

# Intercept in each layer
B_1a = matrix(-0.7,nrow=1,ncol=1)
B_2a = matrix(0.5,nrow=1,ncol=1)
B_3a = matrix(-4.9,nrow=1,ncol=1)

# Merge into lists
B1 = list(B_1a)
B2 = list(B_2a)
B3 = list(B_3a)

B = list(B1, B2, B3)
 
groups = data.frame(Intercept=as.numeric(factor(groups_1)))
groups_f = data.frame(Intercept=factor(groups_1))

#################################################### Simulate network
G = simulate_multiplex_network(
  N_id = N_id,            
  N_layers = N_layers,                   
  B = B,                       
  V = 1,       
  groups = groups,                     
  sr_mu = sr_mu,            
  sr_sigma = sr_sigma,                        
  sr_Rho = sr_Rho,                     
  dr_mu = dr_mu,                            
  dr_sigma = dr_sigma,                         
  dr_Rho = dr_Rho,                          
  outcome_mode="binomial",    
  link_mode = "logit",            
  individual_predictors = data.frame(Mass=Mass, 
                                     Age=Age, 
                                     Strength=Strength),    
  dyadic_predictors = dyadic_preds,        
  individual_effects = sr_effects,        
  dyadic_effects = dr_effects           
 )

 # Tie strength in layer 3 is now essentially undirected up to a small amount of rng error from 0.98 instead of 1.0 in generative model   
par(mfrow=c(1,3)) 
image(G$tie_strength[1,,]-t(G$tie_strength[1,,]))
image(G$tie_strength[2,,]-t(G$tie_strength[2,,]))
image(G$tie_strength[3,,]-t(G$tie_strength[3,,]))

# As a hack to make layer 3 fully undirected, while keeping the current multiplex strucuture, we can add each layer to itself, but transpose layer 3:
G$network2 = G$network
G$samps2 = G$exposure

G$network2[1,,] = G$network[1,,] + (G$network[1,,])  # Dont symmeterize
G$network2[2,,] = G$network[2,,] + (G$network[2,,])  # Dont symmeterize
G$network2[3,,] = G$network[3,,] + t(G$network[3,,]) # Fully symmeterize

G$samps2[1,,] = G$exposure[1,,] + (G$exposure[1,,])        # Dont symmeterize
G$samps2[2,,] = G$exposure[2,,] + (G$exposure[2,,])        # Dont symmeterize
G$samps2[3,,] = G$exposure[3,,] + t(G$exposure[3,,])       # Fully symmeterize

######################## Pre-symmeterized
par(mfrow=c(1,3))
image(G$network[1,,])
image(G$network[2,,])
image(G$network[3,,])

######################## Post-symmeterized
par(mfrow=c(1,3))
image(G$network2[1,,])
image(G$network2[2,,])
image(G$network2[3,,])

#################################################### Create the STRAND data object
outcome = list(Feeding = G$network2[1,,], Fighting = G$network2[2,,], Coresidence = G$network2[3,,])
exposure = list(Feeding = G$samps2[1,,], Fighting = G$samps2[2,,], Coresidence = G$samps2[3,,])

dyad = list(Kinship = Kinship, 
            Dominance = Dominance
            )

indiv =  data.frame(Mass=Mass, 
                    Age=Age, 
                    Strength=Strength
                     )

### Col and row names are now a soft requirement
# Ccan turn off with check_data_organization = FALSE, but its reccomended to always run checks on row and col names 
# when applying STRAND models on real emprical data, just to be sure that everything is correctly organized.
labels = paste("Ind", 1:N_id)
colnames(outcome$Feeding) = rownames(outcome$Feeding) = labels
colnames(outcome$Fighting) = rownames(outcome$Fighting) = labels
colnames(outcome$Coresidence) = rownames(outcome$Coresidence) = labels

colnames(exposure$Feeding) = rownames(exposure$Feeding) = labels
colnames(exposure$Fighting) = rownames(exposure$Fighting) = labels
colnames(exposure$Coresidence) = rownames(exposure$Coresidence) = labels

colnames(dyad$Kinship) = rownames(dyad$Kinship) = labels
colnames(dyad$Dominance) = rownames(dyad$Dominance) = labels

rownames(indiv) = labels


############# Build data
# STRAND will give a warning that a layer is undirected. This doesnt mean that anything is wrong, it is just
# to alert users to check a few features of their model. 
dat = make_strand_data(outcome = outcome,
                       exposure = exposure,
                       block_covariates = NULL, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       outcome_mode="binomial",
                       link_mode="logit",
                       multiplex = TRUE)


# Model
fit = fit_multiplex_model(data=dat,
                          block_regression = ~ 1,
                          focal_regression = ~ Mass + Age + Strength,   # These two regressions (focal and target) should probably be the same for the undirected model to make sense
                          target_regression = ~ Mass + Age + Strength,  #
                          dyad_regression = ~ Kinship + Dominance,      # Dyadic predictors should probably be undirected if used to predict other undirected layers.
                          mode="mcmc",
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 700, iter_sampling = 700,
                                                        max_treedepth = NULL, adapt_delta = 0.95)
)

res = summarize_strand_results(fit)


######################################################### Visualize results
df_plt = res$summary
df_plt = df_plt[which(!df_plt$Variable %in% c("error sd - Feeding", "error sd - Fighting", "error sd - Coresidence")),]

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

df_plt$TrueValues = c(
 sr_sigma[1:3],
 sr_1[1,],
 sr_2[1,],
 sr_3[1,],
 sr_sigma[4:6],
 sr_1[2,],
 sr_2[2,],
 sr_3[2,],
 dr_sigma,
 dr_effects[[1]],
 dr_effects[[2]],
 dr_effects[[3]],
 recip_to_long(sr_Rho),
 recip_to_long(dr_Rho),
 c(t(B_1a)) - mean(c(B_1a)),
 c(t(B_2a)) - mean(c(B_2a)),
 c(t(B_3a)) - mean(c(B_3a))
  )

df_plt$Outcome2 = c(
 c("Feeding","Fighting","Coresidence"),
 c("Feeding","Feeding","Feeding"),
 c("Fighting","Fighting","Fighting"),
 c("Coresidence","Coresidence","Coresidence"),
 c("Feeding","Fighting","Coresidence"),
 c("Feeding","Feeding","Feeding"),
 c("Fighting","Fighting","Fighting"),
 c("Coresidence","Coresidence","Coresidence"),
 c("Feeding","Fighting","Coresidence"),
 c("Feeding","Feeding"),
 c("Fighting","Fighting"),
 c("Coresidence","Coresidence"),
 rep("Recip", 30),
 rep("Other", 3)
  )

df_plt$Variable2 = df_plt$Variable

colnames(df_plt) = c("Variable", "Median", "LI", "HI", "Mean", "SD", "P", "TrueValues","Outcome2", "Variable2")

df_plt$Median = as.numeric(df_plt$Median)
df_plt$LI = as.numeric(df_plt$LI)
df_plt$HI = as.numeric(df_plt$HI)

X = mean(df_plt$Median[64])
df_plt$Median[64] = df_plt$Median[64] - X
df_plt$LI[64] = df_plt$LI[64] - X
df_plt$HI[64] = df_plt$HI[64] - X


X = mean(df_plt$Median[65])
df_plt$Median[65] = df_plt$Median[65] - X
df_plt$LI[65] = df_plt$LI[65] - X
df_plt$HI[65] = df_plt$HI[65] - X


X = mean(df_plt$Median[66])
df_plt$Median[66] = df_plt$Median[66] - X
df_plt$LI[66] = df_plt$LI[66] - X
df_plt$HI[66] = df_plt$HI[66] - X



df_plt$Type = ifelse(str_detect(df_plt$Variable, "focal"), "Focal",
             ifelse(str_detect(df_plt$Variable, "target"), "Target",
             ifelse(str_detect(df_plt$Variable, "dyadic"), "Dyadic",
                    NA)))

df_plt$Type = factor(df_plt$Type)
df_plt$Type = factor(df_plt$Type, levels=c("Focal", "Target", "Dyadic"))

df_plt$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", df_plt$Variable)
df_plt$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", df_plt$Variable)
df_plt$Variable = gsub("dyadic effects coeffs, ", "", df_plt$Variable)
df_plt$Variable = gsub("dyadic effects ", "", df_plt$Variable)
df_plt$Variable = gsub("focal effects ", "", df_plt$Variable)
df_plt$Variable = gsub("target effects ", "", df_plt$Variable)


df_plt$Outcome = ifelse(str_detect(df_plt$Variable, "Feeding"), "Feeding",
             ifelse(str_detect(df_plt$Variable, "Fighting"), "Fighting",
             ifelse(str_detect(df_plt$Variable, "Coresidence"), "Coresidence",
                    NA)))

df_plt$Outcome = factor(df_plt$Outcome)
df_plt$Outcome = factor(df_plt$Outcome, levels=c("Fighting", "Feeding", "Coresidence"))

df_plt$Variable = gsub("Fighting - ", "", df_plt$Variable)
df_plt$Variable = gsub("Feeding - ", "", df_plt$Variable)
df_plt$Variable = gsub("Coresidence - ", "", df_plt$Variable)

df_plt$Variable = gsub(" - Fighting", "", df_plt$Variable)
df_plt$Variable = gsub(" - Feeding", "", df_plt$Variable)
df_plt$Variable = gsub(" - Coresidence", "", df_plt$Variable)

df_plt$Variable = gsub("offset, ", "", df_plt$Variable)

df_plt$Variable = gsub("sd", "SD", df_plt$Variable)


df_plt$Block = ifelse(str_detect(df_plt$Variable, "Any"), "Intercept",
             ifelse(str_detect(df_plt$Variable, "Mottled"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Striped"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Spotted"), "Pattern",
             ifelse(str_detect(df_plt$Variable, "Male"), "Sex",
             ifelse(str_detect(df_plt$Variable, "Female"), "Sex",
                    NA))))))


########################## Plot 1
main_df = df_plt[which(df_plt$Outcome %in% c("Feeding","Fighting","Coresidence")),]
main_df = df_plt[which(df_plt$Type %in% c("Focal","Target","Dyadic")),]

main_df$Variable = factor(main_df$Variable)
main_df$Variable = factor(main_df$Variable, levels=rev(c("SD", "Age", "Mass", "Strength", "Dominance", "Kinship")))

p = ggplot(main_df, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Outcome)) + 
           geom_linerange(size = 1, color=colors[4]) + 
           geom_point(size = 2, color=colors[4]) +
           geom_point(size = 2, aes(x = Variable, y = TrueValues, group=Outcome), color=colors[2], shape=18) +
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

########################## Plot 3
recip_df = df_plt[which(df_plt$Outcome2 == "Recip"),]

recip_df$Type = ifelse(str_detect(recip_df$Variable2, "Generalized"), "Generalized",
                ifelse(str_detect(recip_df$Variable2, "Dyadic"), "Dyadic",
                    NA))

recip_df$Variable2 = gsub("Dyadic reciprocity - ", "", recip_df$Variable2)
recip_df$Variable2 = gsub("Generalized reciprocity - ", "", recip_df$Variable2)

p1 = ggplot(recip_df[which(recip_df$Type == "Generalized"),], aes(x = Variable2, y = Median, ymin = LI, ymax = HI, group=Outcome)) + 
           geom_linerange(size = 1, color=colors[4]) + 
           geom_point(size = 2, color=colors[4]) +
           geom_point(size = 2, aes(x = Variable2, y = TrueValues, group=Outcome), color=colors[2], shape=18) +
           facet_grid(. ~ Type, scales = "free", space = "free") + 
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
p1
# We recover generalized reciprocity, within and between layers, too 


p2 = ggplot(recip_df[which(recip_df$Type == "Dyadic"),], aes(x = Variable2, y = Median, ymin = LI, ymax = HI, group=Outcome)) + 
           geom_linerange(size = 1, color=colors[4]) + 
           geom_point(size = 2, color=colors[4]) +
           geom_point(size = 2, aes(x = Variable2, y = TrueValues, group=Outcome), color=colors[2], shape=18) +
           facet_grid(. ~ Type, scales = "free", space = "free") + 
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
p2
# We recover dyadic reciprocity, within and between layers, too 






