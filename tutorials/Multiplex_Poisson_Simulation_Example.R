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

library(STRAND)
library(stringr)
library(ggplot2)
library(psych)
library(rethinking)

# Make data
 N_id = 150     # Individuals in network
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
B_1a = matrix(-3.7,nrow=1,ncol=1)
B_2a = matrix(0.5,nrow=1,ncol=1)
B_3a = matrix(-3.5,nrow=1,ncol=1)

# Offsets for Pattern
B_1b = matrix(rnorm(9,0,1),nrow=3,ncol=3)
B_2b = matrix(rnorm(9,-2,1),nrow=3,ncol=3)
B_3b = matrix(rnorm(9,0.5,1),nrow=3,ncol=3)

diag(B_1b) = diag(B_1b) - 0
diag(B_2b) = diag(B_2b) + 1.5
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

dat = make_strand_data(outcome = outcome,
                       block_covariates = groups, 
                       individual_covariates = indiv, 
                       dyadic_covariates = dyad,
                       exposure = exposure,
                       outcome_mode="poisson",
                       link_mode="log",
                       multiplex = TRUE)


# Model
fit = fit_multiplex_model(data=dat,
                          block_regression = ~ Pattern + Sex,
                          focal_regression = ~ Mass + Age + Strength,
                          target_regression = ~ Mass + Age + Strength,
                          dyad_regression = ~ Kinship + Dominance,
                          mode="mcmc",
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = NULL, adapt_delta = 0.95)
)

res = summarize_multiplex_bsrm_results(fit)


######################################################### Visualize results
df_plt = res$summary
df_plt = df_plt[which(!df_plt$Variable %in% c("error sd - Feeding", "error sd - Fighting", "error sd - Grooming")),]

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
 c(t(B_1b)) - mean(c(B_1b)), 
 c(t(B_1c)) - mean(c(B_1c)),
 c(t(B_2a)) - mean(c(B_2a)),
 c(t(B_2b)) - mean(c(B_2b)), 
 c(t(B_2c)) - mean(c(B_2c)),
 c(t(B_3a)) - mean(c(B_3a)),
 c(t(B_3b)) - mean(c(B_3b)), 
 c(t(B_3c)) - mean(c(B_3c))
  )

df_plt$Outcome2 = c(
 c("Feeding","Fighting","Grooming"),
 c("Feeding","Feeding","Feeding"),
 c("Fighting","Fighting","Fighting"),
 c("Grooming","Grooming","Grooming"),
 c("Feeding","Fighting","Grooming"),
 c("Feeding","Feeding","Feeding"),
 c("Fighting","Fighting","Fighting"),
 c("Grooming","Grooming","Grooming"),
 c("Feeding","Fighting","Grooming"),
 c("Feeding","Feeding"),
 c("Fighting","Fighting"),
 c("Grooming","Grooming"),
 rep("Recip", 30),
 rep("Other", 42)
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

X = mean(df_plt$Median[65:73])
df_plt$Median[65:73] = df_plt$Median[65:73] - X
df_plt$LI[65:73] = df_plt$LI[65:73] - X
df_plt$HI[65:73] = df_plt$HI[65:73] - X

X = mean(df_plt$Median[74:77])
df_plt$Median[74:77] = df_plt$Median[74:77] - X
df_plt$LI[74:77] = df_plt$LI[74:77] - X
df_plt$HI[74:77] = df_plt$HI[74:77] - X


X = mean(df_plt$Median[78])
df_plt$Median[78] = df_plt$Median[78] - X
df_plt$LI[78] = df_plt$LI[78] - X
df_plt$HI[78] = df_plt$HI[78] - X

X = mean(df_plt$Median[79:87])
df_plt$Median[79:87] = df_plt$Median[79:87] - X
df_plt$LI[79:87] = df_plt$LI[79:87] - X
df_plt$HI[79:87] = df_plt$HI[79:87] - X

X = mean(df_plt$Median[88:91])
df_plt$Median[88:91] = df_plt$Median[88:91] - X
df_plt$LI[88:91] = df_plt$LI[88:91] - X
df_plt$HI[88:91] = df_plt$HI[88:91] - X


X = mean(df_plt$Median[92])
df_plt$Median[92] = df_plt$Median[92] - X
df_plt$LI[92] = df_plt$LI[92] - X
df_plt$HI[92] = df_plt$HI[92] - X

X = mean(df_plt$Median[93:101])
df_plt$Median[93:101] = df_plt$Median[93:101] - X
df_plt$LI[93:101] = df_plt$LI[93:101] - X
df_plt$HI[93:101] = df_plt$HI[93:101] - X

X = mean(df_plt$Median[102:105])
df_plt$Median[102:105] = df_plt$Median[102:105] - X
df_plt$LI[102:105] = df_plt$LI[102:105] - X
df_plt$HI[102:105] = df_plt$HI[102:105] - X


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
             ifelse(str_detect(df_plt$Variable, "Grooming"), "Grooming",
                    NA)))

df_plt$Outcome = factor(df_plt$Outcome)
df_plt$Outcome = factor(df_plt$Outcome, levels=c("Fighting", "Feeding", "Grooming"))

df_plt$Variable = gsub("Fighting - ", "", df_plt$Variable)
df_plt$Variable = gsub("Feeding - ", "", df_plt$Variable)
df_plt$Variable = gsub("Grooming - ", "", df_plt$Variable)

df_plt$Variable = gsub(" - Fighting", "", df_plt$Variable)
df_plt$Variable = gsub(" - Feeding", "", df_plt$Variable)
df_plt$Variable = gsub(" - Grooming", "", df_plt$Variable)

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
main_df = df_plt[which(df_plt$Outcome %in% c("Feeding","Fighting","Grooming")),]
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

# ggsave("sim_res.pdf",p, width=9, height=4.5)

########################## Plot 2
block_df = df_plt[which(df_plt$Outcome2 == "Other" & df_plt$Block != "Intercept"),]

p = ggplot(block_df, aes(x = Variable, y = Median, ymin = LI, ymax = HI, group=Outcome)) + 
           geom_linerange(size = 1, color=colors[4]) + 
           geom_point(size = 2, color=colors[4]) +
           geom_point(size = 2, aes(x = Variable, y = TrueValues, group=Outcome), color=colors[2], shape=18) +
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

# ggsave("sim_res_block.pdf",p, width=9, height=4.5)


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

# ggsave("sim_res_gen.pdf",p1, width=6, height=6)

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

# ggsave("sim_res_dyad.pdf",p2, width=6, height=6)






