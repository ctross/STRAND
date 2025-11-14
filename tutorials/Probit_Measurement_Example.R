################################################################################################################
#
# How does a latent network probit SRM compare to a standard probit SRM ?
#
################################################################################################################

# STRAND is designed to estimate the parameters of the true network, even from single-layer data with measurement error. 
# Treating observations/reports as representing the "true" latent social network
# leads to biased estimates of dyadic reciprocity. The bias is: dr_rho^2 / (dr_rho^2 + error_sigma^2).
# In the special case of error_sigma=0, there is no bias. Otherwise, it is better to estimate the error and account for it.

# Lets look at this in detail and compare STRAND to AMEN in the case of a probit model for binary outcomes.
# We will see that STRAND is equivalent to the AMEN biprobit model when we put very strong priors on the relative level of measurement error.
# STRAND is robust to measurement error. AMEN and other biprobit models that don't incorporate latent network estimation are not.

# Helper function to shape posteriors
 Q = function(x){
  return(c(median(x),HPDI(x, 0.99)))
 }

# Load libraries
 library(amen)
 library(igraph)
 library(STRAND)
 library(ggplot2)

set.seed(1)
 N_id = 70
 B = -2.5
 sr_sigma = c(0.6, 0.9)    # Sender and reciever effect SDs in the true network 
 sr_rho = 0.6              # Sender and reciever correlation in the true network
 dr_sigma = 1              # Dyadic effect SD in the true network 
 dr_rho = 0.8              # Dyadic correlation in the true network 

 y_binary_e = y_binary_0 = y_star = dr = matrix(0, nrow=N_id, ncol=N_id)

 # Create correlation matrices 
  Rho_sr = Rho_dr = diag(c(1,1))
  Rho_sr[1,2] = Rho_sr[2,1] = sr_rho
  Rho_dr[1,2] = Rho_dr[2,1] = dr_rho

 # Varying effects on individuals
  sr = matrix(NA, nrow=N_id, ncol=2)

  for( i in 1:N_id){
   sr[i,] = rmvnorm2(1, Mu=c(0,0), sigma=sr_sigma, Rho=Rho_sr)
   } 

 # Dyadic effects
  for ( i in 1:(N_id-1) ){
    for ( j in (i+1):N_id){
     dr_scrap = rmvnorm2(1, Mu=rep(0, 2), sigma=rep(dr_sigma,2), Rho=Rho_dr)

     dr[i,j] = dr_scrap[1] + B
     dr[j,i] = dr_scrap[2] + B

 # Simulate outcomes via thresholding the latent variable
  # True latent tie strength
  y_star[i,j] = sr[i,1] + sr[j,2] + dr[i,j]
  y_star[j,i] = sr[j,1] + sr[i,2] + dr[j,i]
  
  # Binary probit based on latent tie strength with no measurement noise
  y_binary_0[i,j] = ifelse(y_star[i,j] < 0, 0, 1)
  y_binary_0[j,i] = ifelse(y_star[j,i] < 0, 0, 1)

  # Binary probit based on latent tie strength with correlation-0 measurement noise at level sd_noise
  sd_noise = 0.95 # About the same scale as sender, receiver, and dyadic effects

  y_binary_e[i,j] = ifelse(y_star[i,j] + rnorm(1, 0, sd_noise) < 0, 0, 1)
  y_binary_e[j,i] = ifelse(y_star[j,i] + rnorm(1, 0, sd_noise) < 0, 0, 1)
          }
        }

Net = graph_from_adjacency_matrix(y_binary_0, mode = c("directed"))
par(mfrow=c(1,2))
plot(Net, edge.arrow.size = 0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
image(y_binary_0)

######################################################## True variance partitions in the no noise model
# Note: we use the Koster and Leckie definition of VPCs here (but the dyadic VPC should really be: dr_rho^2 * dr_sigma^2,
# since only a fraction dr_rho^2 of the variance is "explained" by dyadic effects. The other (1-dr_rho^2) can be 
# incorpoarted into the noise term. But we will leave this for another tutorial.)

vcs = c(sr_sigma, dr_sigma, 0)^2
vpcs = vcs/(vcs[1] + vcs[2] + vcs[3] + vcs[4])
d_e = sum(vcs[3:4]/(vcs[1] + vcs[2] + vcs[3] + vcs[4]))

True_Measures_Alt_0 = c(sr_rho,                                            # Sender and reciever correlation in the true network
                        dr_rho,                                            # Dyadic correlation in the true network 
                        vpcs,                                              # Actual VPCs
                        d_e,                                               # Sum of dyadic and error VCPs
                        dr_rho*(dr_sigma^2 /(dr_sigma^2 + 0^2))            # Dyadic correlation in the measured/observed responses 
                        )

######################################################## True variance partitions in the noise model
vcs = c(sr_sigma, dr_sigma, sd_noise)^2
vpcs = vcs/(vcs[1] + vcs[2] + vcs[3] + vcs[4])
d_e = sum(vcs[3:4]/(vcs[1] + vcs[2] + vcs[3] + vcs[4]))

True_Measures_Alt_e = c(sr_rho,                                            # Sender and reciever correlation in the true network
                        dr_rho,                                            # Dyadic correlation in the true network 
                        vpcs,                                              # Actual VPCs
                        d_e,                                               # Sum of dyadic and error VCPs
                        dr_rho*(dr_sigma^2 /(dr_sigma^2 + sd_noise^2))     # Dyadic correlation in the measured/observed responses 
                        )

################################################################################################ Fit AMEN 0
fit_SRM = ame(y_binary_0, family = "bin")
summary(fit_SRM)

## The variance parameters for nodes are reported as variances and covariances
ame_sender_var = fit_SRM$VC[,1]
ame_covariance = fit_SRM$VC[,2]
ame_target_var = fit_SRM$VC[,3]
ame_dyadic_cor = fit_SRM$VC[,4] 

ame_generalized_reciprocity = ame_covariance / (sqrt (ame_sender_var * ame_target_var))

## One can also calculate the Variance Partition Coefficients
ame_dyadic_var = 1 
ame_total_var = ame_sender_var + ame_target_var + ame_dyadic_var 

ame_sender_VPC = ame_sender_var / ame_total_var
ame_target_VPC = ame_target_var / ame_total_var
ame_dyadic_VPC = ame_dyadic_var / ame_total_var

ame_dyadic_cor_meas = ame_dyadic_cor*(ame_dyadic_var /(ame_dyadic_var + 0^2))

AME_Measures_Alt_0 = rbind(Q(ame_generalized_reciprocity), rep(NA,3), Q(ame_sender_VPC), Q(ame_target_VPC), rep(NA,3), rep(NA,3), Q(ame_dyadic_VPC), Q(ame_dyadic_cor_meas))

colnames(AME_Measures_Alt_0) = c("M","L","H")
rownames(AME_Measures_Alt_0) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")


################################################################################################ Fit AMEN e
fit_SRM = ame(y_binary_e, family = "bin")
summary(fit_SRM)

## The variance parameters for nodes are reported as variances and covariances
ame_sender_var = fit_SRM$VC[,1]
ame_covariance = fit_SRM$VC[,2]
ame_target_var = fit_SRM$VC[,3]
ame_dyadic_cor = fit_SRM$VC[,4] 

ame_generalized_reciprocity = ame_covariance / (sqrt (ame_sender_var * ame_target_var))

## One can also calculate the Variance Partition Coefficients
ame_dyadic_var = 1
ame_total_var = ame_sender_var + ame_target_var + ame_dyadic_var 

ame_sender_VPC = ame_sender_var / ame_total_var
ame_target_VPC = ame_target_var / ame_total_var
ame_dyadic_VPC = ame_dyadic_var / ame_total_var

ame_dyadic_cor_meas = ame_dyadic_cor*(ame_dyadic_var /(ame_dyadic_var + 0^2)) # zero is by fiat here

AME_Measures_Alt_e = rbind(Q(ame_generalized_reciprocity), rep(NA,3), Q(ame_sender_VPC), Q(ame_target_VPC), rep(NA,3), rep(NA,3), Q(ame_dyadic_VPC), Q(ame_dyadic_cor_meas))

colnames(AME_Measures_Alt_e) = c("M","L","H")
rownames(AME_Measures_Alt_e) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")

################################################################################################ Fit STRAND 0
########################################################################### Organize data
# Create the STRAND data object
rownames(y_binary_0) = colnames(y_binary_0) = paste0("ID", 1:N_id)
outcome = list(trade = y_binary_0)

dat1 = make_strand_data(
  self_report = outcome,
  individual_covariates = NULL,
  dyadic_covariates = NULL,
  outcome_mode = "bernoulli",
  link_mode = "probit"
)

########################################## Fit STRAND with amen-like priors on measurement error
fit1 = fit_block_plus_social_relations_model(
    data=dat1 ,
    block_regression = ~ 1,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ 1,
    priors=make_priors(             # We can replicate the AMEN model by setting priors that sender, receiver, and dyadic sds are enormous relative to the noise
       priors_to_change=list(       # Error has a value of 1, so priors on vpcs are now hella strong: c(10^2, 10^2, 10^2, 1)/sum(c(10^2, 10^2, 10^2, 1))
         "B_ingroup"=c(0.01, 10),   # Let intercept slide
         "sr_sigma"=c(10, 3),       # Location of 10
         "dr_sigma"=c(10, 3)        # Location of 10
       )),
    mode="mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 500 ,
      iter_sampling = 500 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )

res = summarize_strand_results(fit1)

### VPCs
STRAND_dyadic_reciprocity = res$samples$srm_model_samples$dyadic_L[,2,1]
STRAND_generalized_reciprocity = res$samples$srm_model_samples$focal_target_L[,2,1]

STRAND_sender_var = (res$samples$srm_model_samples$focal_target_sd[,1])^2
STRAND_target_var = (res$samples$srm_model_samples$focal_target_sd[,2])^2
STRAND_dyadic_var = (res$samples$srm_model_samples$dyadic_sd)^2
STRAND_error_var = 1
STRAND_total_var = STRAND_sender_var + STRAND_target_var + STRAND_dyadic_var + STRAND_error_var 

STRAND_sender_VPC = STRAND_sender_var / STRAND_total_var
STRAND_target_VPC = STRAND_target_var / STRAND_total_var
STRAND_dyadic_VPC = STRAND_dyadic_var / STRAND_total_var
STRAND_error_VPC = STRAND_error_var / STRAND_total_var

STRAND_dyadicerror_VPC = (STRAND_dyadic_var+STRAND_error_var) / STRAND_total_var

STRAND_dyadic_reciprocity_meas = STRAND_dyadic_reciprocity*(STRAND_dyadic_var /(STRAND_dyadic_var + STRAND_error_var))

STRAND_Measures_StrongPrior_0 = rbind(Q(STRAND_generalized_reciprocity), Q(STRAND_dyadic_reciprocity), Q(STRAND_sender_VPC), Q(STRAND_target_VPC), Q(STRAND_dyadic_VPC), Q(STRAND_error_VPC), Q(STRAND_dyadicerror_VPC), Q(STRAND_dyadic_reciprocity_meas))

colnames(STRAND_Measures_StrongPrior_0) = c("M","L","H")
rownames(STRAND_Measures_StrongPrior_0) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")

strand_VPCs(fit1, n_partitions = 3)
strand_VPCs(fit1, n_partitions = 3, include_reciprocity=TRUE)

########################################## Now fit STRAND with default priors which treat measurement error as having the same scale as the other effects a priori
fit2 = fit_block_plus_social_relations_model(
    data=dat1 ,
    block_regression = ~ 1,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ 1,
        priors=make_priors(         # More reasonable priors that treat measurement error as having similar scale to other effects a priori
       priors_to_change=list(       # 
         "B_ingroup"=c(0.001, 10),  # 
         "sr_sigma"=c(0, 1),        # 
         "dr_sigma"=c(0, 1)         # 
       )),
    mode="mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 500 ,
      iter_sampling = 500 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )

res = summarize_strand_results(fit2)

### VPCs
STRAND_dyadic_reciprocity = res$samples$srm_model_samples$dyadic_L[,2,1]
STRAND_generalized_reciprocity = res$samples$srm_model_samples$focal_target_L[,2,1]

STRAND_sender_var = (res$samples$srm_model_samples$focal_target_sd[,1])^2
STRAND_target_var = (res$samples$srm_model_samples$focal_target_sd[,2])^2
STRAND_dyadic_var = (res$samples$srm_model_samples$dyadic_sd)^2
STRAND_error_var = 1
STRAND_total_var = STRAND_sender_var + STRAND_target_var + STRAND_dyadic_var + STRAND_error_var 

STRAND_sender_VPC = STRAND_sender_var / STRAND_total_var
STRAND_target_VPC = STRAND_target_var / STRAND_total_var
STRAND_dyadic_VPC = STRAND_dyadic_var / STRAND_total_var
STRAND_error_VPC = STRAND_error_var / STRAND_total_var

STRAND_dyadicerror_VPC = (STRAND_dyadic_var+STRAND_error_var) / STRAND_total_var

STRAND_dyadic_reciprocity_meas = STRAND_dyadic_reciprocity*(STRAND_dyadic_var /(STRAND_dyadic_var + STRAND_error_var))

STRAND_Measures_WeakPrior_0 = rbind(Q(STRAND_generalized_reciprocity), Q(STRAND_dyadic_reciprocity), Q(STRAND_sender_VPC), Q(STRAND_target_VPC), Q(STRAND_dyadic_VPC), Q(STRAND_error_VPC), Q(STRAND_dyadicerror_VPC), Q(STRAND_dyadic_reciprocity_meas))

colnames(STRAND_Measures_WeakPrior_0) = c("M","L","H")
rownames(STRAND_Measures_WeakPrior_0) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")

strand_VPCs(fit2, n_partitions = 3)
strand_VPCs(fit2, n_partitions = 3, include_reciprocity=TRUE)

################################################################################################ Fit STRAND e
########################################################################### Organize data
# Create the STRAND data object
rownames(y_binary_e) = colnames(y_binary_e) = paste0("ID", 1:N_id)
outcome = list(trade = y_binary_e)

dat2 = make_strand_data(
  self_report = outcome,
  individual_covariates = NULL,
  dyadic_covariates = NULL,
  outcome_mode = "bernoulli",
  link_mode = "probit"
)

########################################## Fit STRAND with amen-like priors on measurement error
fit3 = fit_block_plus_social_relations_model(
    data=dat2 ,
    block_regression = ~ 1,,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ 1,
    priors=make_priors(            # We can replicate the AMEN model by setting priors that sender, receiver, and dyadic sds are enormous relative to the noise
       priors_to_change=list(      # Error has a value of 1, so priors on vpcs are now hella strong: c(10^2, 10^2, 10^2, 1)/sum(c(10^2, 10^2, 10^2, 1))
         "B_ingroup"=c(0.01, 10),  # Let intercept slide
         "sr_sigma"=c(10, 3),      # Location of 10
         "dr_sigma"=c(10, 3)       # Location of 10
       )),
    mode="mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 500 ,
      iter_sampling = 500 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )

res = summarize_strand_results(fit3)

### VPCs
STRAND_dyadic_reciprocity = res$samples$srm_model_samples$dyadic_L[,2,1]
STRAND_generalized_reciprocity = res$samples$srm_model_samples$focal_target_L[,2,1]

STRAND_sender_var = (res$samples$srm_model_samples$focal_target_sd[,1])^2
STRAND_target_var = (res$samples$srm_model_samples$focal_target_sd[,2])^2
STRAND_dyadic_var = (res$samples$srm_model_samples$dyadic_sd)^2
STRAND_error_var = 1
STRAND_total_var = STRAND_sender_var + STRAND_target_var + STRAND_dyadic_var + STRAND_error_var 

STRAND_sender_VPC = STRAND_sender_var / STRAND_total_var
STRAND_target_VPC = STRAND_target_var / STRAND_total_var
STRAND_dyadic_VPC = STRAND_dyadic_var / STRAND_total_var
STRAND_error_VPC = STRAND_error_var / STRAND_total_var

STRAND_dyadicerror_VPC = (STRAND_dyadic_var+STRAND_error_var) / STRAND_total_var

STRAND_dyadic_reciprocity_meas = STRAND_dyadic_reciprocity*(STRAND_dyadic_var /(STRAND_dyadic_var + STRAND_error_var))

STRAND_Measures_StrongPrior_e = rbind(Q(STRAND_generalized_reciprocity), Q(STRAND_dyadic_reciprocity), Q(STRAND_sender_VPC), Q(STRAND_target_VPC), Q(STRAND_dyadic_VPC), Q(STRAND_error_VPC), Q(STRAND_dyadicerror_VPC), Q(STRAND_dyadic_reciprocity_meas))

colnames(STRAND_Measures_StrongPrior_e) = c("M","L","H")
rownames(STRAND_Measures_StrongPrior_e) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")

strand_VPCs(fit3, n_partitions = 3)
strand_VPCs(fit3, n_partitions = 3, include_reciprocity=TRUE)

########################################## Now fit STRAND with default priors on measurement error which assume it to be of the same scale as other effects
fit4 = fit_block_plus_social_relations_model(
    data=dat2 ,
    block_regression = ~ 1,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ 1,
        priors=make_priors(         # More reasonable priors that treat measurement error as having similar scale to other effects a priori
       priors_to_change=list(       # 
         "B_ingroup"=c(0.01, 10),   # 
         "sr_sigma"=c(0, 1),        # 
         "dr_sigma"=c(0, 1)         # 
       )),
    mode="mcmc",
    mcmc_parameters = list(
      chains = 1,
      iter_warmup = 1500 ,
      iter_sampling = 1500 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )

res = summarize_strand_results(fit4)

### VPCs
STRAND_dyadic_reciprocity = res$samples$srm_model_samples$dyadic_L[,2,1]
STRAND_generalized_reciprocity = res$samples$srm_model_samples$focal_target_L[,2,1]

STRAND_sender_var = (res$samples$srm_model_samples$focal_target_sd[,1])^2
STRAND_target_var = (res$samples$srm_model_samples$focal_target_sd[,2])^2
STRAND_dyadic_var = (res$samples$srm_model_samples$dyadic_sd)^2
STRAND_error_var = 1
STRAND_total_var = STRAND_sender_var + STRAND_target_var + STRAND_dyadic_var + STRAND_error_var 

STRAND_sender_VPC = STRAND_sender_var / STRAND_total_var
STRAND_target_VPC = STRAND_target_var / STRAND_total_var
STRAND_dyadic_VPC = STRAND_dyadic_var / STRAND_total_var
STRAND_error_VPC = STRAND_error_var / STRAND_total_var

STRAND_dyadicerror_VPC = (STRAND_dyadic_var+STRAND_error_var) / STRAND_total_var

STRAND_dyadic_reciprocity_meas = STRAND_dyadic_reciprocity*(STRAND_dyadic_var /(STRAND_dyadic_var + STRAND_error_var))

STRAND_Measures_WeakPrior_e = rbind(Q(STRAND_generalized_reciprocity), Q(STRAND_dyadic_reciprocity), Q(STRAND_sender_VPC), Q(STRAND_target_VPC), Q(STRAND_dyadic_VPC), Q(STRAND_error_VPC), Q(STRAND_dyadicerror_VPC), Q(STRAND_dyadic_reciprocity_meas))

colnames(STRAND_Measures_WeakPrior_e) = c("M","L","H")
rownames(STRAND_Measures_WeakPrior_e) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")

strand_VPCs(fit4, n_partitions = 3)
strand_VPCs(fit4, n_partitions = 3, include_reciprocity=TRUE)
  
###################################################################################################### Plot all
True_Measures_Alt_0_merge = data.frame(M=True_Measures_Alt_0,L=NA,H=NA)
True_Measures_Alt_e_merge = data.frame(M=True_Measures_Alt_e,L=NA,H=NA)

AME_Measures_Alt_0_merge = as.data.frame(AME_Measures_Alt_0)
AME_Measures_Alt_e_merge = as.data.frame(AME_Measures_Alt_e)

STRAND_Measures_StrongPrior_0_merge = as.data.frame(STRAND_Measures_StrongPrior_0)
STRAND_Measures_StrongPrior_e_merge = as.data.frame(STRAND_Measures_StrongPrior_e)

STRAND_Measures_WeakPrior_0_merge = as.data.frame(STRAND_Measures_WeakPrior_0)
STRAND_Measures_WeakPrior_e_merge = as.data.frame(STRAND_Measures_WeakPrior_e)

AME_Measures_Alt_0_merge$Package="AMEN"
STRAND_Measures_StrongPrior_0_merge$Package="STRAND_StrongPrior"
STRAND_Measures_WeakPrior_0_merge$Package="STRAND_WeakPrior"
True_Measures_Alt_0_merge$Package="True Value"

AME_Measures_Alt_e_merge$Package="AMEN"
STRAND_Measures_StrongPrior_e_merge$Package="STRAND_StrongPrior"
STRAND_Measures_WeakPrior_e_merge$Package="STRAND_WeakPrior"
True_Measures_Alt_e_merge$Package="True Value"

AME_Measures_Alt_0_merge$NoiseLevel="0"
STRAND_Measures_StrongPrior_0_merge$NoiseLevel="0"
STRAND_Measures_WeakPrior_0_merge$NoiseLevel="0"
True_Measures_Alt_0_merge$NoiseLevel="0"

AME_Measures_Alt_e_merge$NoiseLevel="1.5"
STRAND_Measures_StrongPrior_e_merge$NoiseLevel="1.5"
STRAND_Measures_WeakPrior_e_merge$NoiseLevel="1.5"
True_Measures_Alt_e_merge$NoiseLevel="1.5"


Measures_Alt = rbind(AME_Measures_Alt_0_merge, 
                 STRAND_Measures_StrongPrior_0_merge,
                 STRAND_Measures_WeakPrior_0_merge, 
                 True_Measures_Alt_0_merge,
                 AME_Measures_Alt_e_merge,
                 STRAND_Measures_StrongPrior_e_merge,
                 STRAND_Measures_WeakPrior_e_merge,
                 True_Measures_Alt_e_merge
                 )

Measures_Alt$Variable = rep(rownames(AME_Measures_Alt_0_merge),8)

Measures_Alt$Type = ifelse(Measures_Alt$Variable %in% c("Gen. Recip.","Dyad. Recip.","Meas. Dyad. Recip."), "Reciprocity", "VPC")

pal = c("#0D205C", "#507BAF", "#C5570E", "#DF913E", "#13667B", "#35ABC0", 
"#E39913", "#F3BB57", "#721712", "#AF5249", "#0B231A", "#4E7A65", 
"#431008", "#A56137")

Measures_Alt$Variable = factor(Measures_Alt$Variable)
Measures_Alt$Variable = factor(Measures_Alt$Variable, levels=levels(Measures_Alt$Variable)[rev(c(6,1,5,7,2,4,3,8))])

p1 = ggplot(Measures_Alt,aes(x=Variable,y=M,ymin=L,ymax=H,color=Package))+ 
     geom_linerange(size=1,aes(color=Package), position = position_dodge2(width=0.5,preserve = "total",padding = 1))+
     geom_point(size=2,aes(color=Package), position = position_dodge2(width=0.5,preserve = "total",padding = 1))+ facet_grid(Type~NoiseLevel,scales="free",space = "free_y")+
     labs(y="Estimated value", x="") + theme(strip.text.x = element_text(size=17,face="bold"), 
     strip.text.y = element_text(size=17,face="bold"),axis.text=element_text(size=17),axis.title=element_text(size=17,
     face="bold"))+theme(strip.text.y = element_text(angle = 360),legend.title=element_text(size=17), 
    legend.text=element_text(size=17))  + coord_flip() + theme(panel.spacing = unit(1, "lines")) +
     scale_color_manual(values=c("True Value"="black", "AMEN"="darkred", "STRAND_StrongPrior"=pal[7], "STRAND_WeakPrior"=pal[5]
                 )) +
          scale_linetype_manual(values=c("None"="solid", "ME"="dashed")) +
    theme(legend.text=element_text(size=17)) 

 p1

###############################################################################################################################

# Results:

# True dyadic reciprocity was set at 0.8 in the generative model, but amen, and the STRAND model with priors that force it to behave like amen, estimate dyadic reciprocity to be about 0.4. Why?
# Measurement error introduces correlation-zero noise into the dyadic observations, masking the true level of reciprocity in the real latent network. STRAND integrates over all causal paths of producing
# observations with a dyadic correlation of 0.4. This includes the extreme endpoint of a true dyadic reciprocity of 0.4, and no measurement error, as well as generative processes with higher levels of dyadic 
# reciprocity coupled with higher levels of measurement error. As such, STRANDs posteriors on dyadic_rho and dyadic_sigma are much wider, better reflecting the true level of uncertainty we have in the 
# parameters, when noise in measurement is possible.

# Because measurement error can mask dyadic reciprocity, we can only learn a lower bound on the magnitude of dyadic reciprocity from the data. STRAND represents this, while amen doesnt.

# Note that STRAND with either strong or weak priors for the level of measurement error recovers the correct measurement-level reciprocity (whether we set the sd of noise in the generative model to 0 or 1.5). 
# With weak priors, STRAND actually gets the latent-network reciprocity estimate correct too (amen doesnt). Using strong priors that measurement error is zero is normally a bad idea.

# Other notes: 

# 1) One can always recover amen-style 3-way variance partitions from a STRAND model just by summing the dyadic and error VPCs: Or by using strand_VPCs(fit2, n_partitions = 3)

# 2) One can always recover the amen-style, measurement-biased, dyadic reciprocity estimate from a STRAND model just by multiplying dr_rho by: dr_sigma^2 / (dr_sigma^2 + error_sigma^2) 
#  or by running strand_VPCs(fit2, n_partitions = 3, include_reciprocity=TRUE, mode="adj")

# 3) If users have special reason to believe that measurement error is zero a priori (e.g., data are from twitter ties), then its fine to use strong, amen-like, priors: "sr_sigma"=c(10, 3) and "dr_sigma"=c(10, 3).
# In most cases, using weak priors is a better idea.
