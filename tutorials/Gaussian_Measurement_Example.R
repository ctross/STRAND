################################################################################################################
#
# How does a latent network SRM compare to a standard SRM ?
#
################################################################################################################

# Some have argued that STRAND SRM models are not identifiable because we have unit-level random effects and unit-level residual error. 
# We have three reasons for rejecting this claim:
 
# 1) In general, Bayes' theorem states that if Pr(A) is proper, then (Pr(B|A)Pr(A))/Pr(B) = P(A|B) is proper with probability 1. 
# In STRAND, we use proper priors, and so our models have proper posteriors. Thus, our models are always identified. ■

# 2) Setting aside the guarantees of Bayes' theorem and proper priors, univariate outcome models of the form: y_[i] = a + b_[i] + e_[i], with y_[i] being a 1-vector/scalar 
# outcome, and with b_[i] ~ normal(0, sd_b) being random effects and e_[i] ~ normal(0, sd_e) being residuals, do suffer a clear invariance if priors are improper, since:
#  (b_[i] + e_[i]) ~ normal(0, sqrt(sd_b^2 + sd_e^2)), and sd_b^2 + sd_e^2 are fully exchangeable (or indeed continuously mixable) without affecting the shape of the outcome distribution.
sigma_e = 13
sigma_b = 3.14

b = rnorm(10000, 0, sigma_b)
e = rnorm(10000, 0, sigma_e)

y1 = b + e

# swap the sigmas
sigma_e = 3.14
sigma_b = 13

b = rnorm(10000, 0, sigma_b)
e = rnorm(10000, 0, sigma_e)

y2 = b + e

plot(density(y1))
lines(density(y2), col="blue") # Same densities

# This fact, however, does not immediately generalize to the 2-dimensional outcome case. Bivariate outcome models of the form: y_[i] = a + b_[i] + e_[i], with y_[i] being a 2-vector
# outcome, and with b_[i] ~ MultivariateNormal(c(0,0), sd_b %*% Rho %*% sd_b) being correlated random effects, and e_[i] ~ MultivariateNormal(c(0,0), sd_e %*% I %*% sd_e) being uncorrelated residuals, 
# do not suffer the same invariance, since (b_[i] + e_[i]) ~ MultivariateNormal(c(0,0), sd_b %*% Rho %*% sd_b  +  sd_e %*% I %*% sd_e). Here, sd_b affects shape of the outcome distribution 
# by hitting the off-diagonal of Rho, while sd_e doesn't affect the off-diagonal of I (which is the identity matrix, and thus has an off-diagonal of 0).


# ########################## Can you exchange sigma_e and sigma_b now?
sigma_e = c(13, 13)
sigma_b = c(3.14, 3.14)
rho = 0.8
Sigma_b = diag(sigma_b^2)
Sigma_b[1,2] = Sigma_b[2,1] = rho*sigma_b[1]^2
Sigma_e = diag(sigma_e^2)
y1 = amen::rmvnorm(10000, mu = c(0,0), Sigma = Sigma_b + Sigma_e)

# ########################## swap the sigmas
sigma_e = c(3.14, 3.14)
sigma_b = c(13, 13)
rho = 0.9
Sigma_b = diag(sigma_b^2)
Sigma_b[1,2] = Sigma_b[2,1] = rho*sigma_b[1]^2
Sigma_e = diag(c(sigma_e^2))
y2 = amen::rmvnorm(10000, mu = c(0,0), Sigma = Sigma_b + Sigma_e)

plot(y1[,1], y1[,2])
points(y2[,1], y2[,2], col="blue") # Not the same joint densities

# There is still an invariance lurking though, since the observation-level correlation implied by (b_[i] + e_[i]) ~ MultivariateNormal(c(0,0),  sd_b %*% Rho %*% sd_b  +  sd_e %*% I %*% sd_e)
# is rho*(sigma_b^2)/(sigma_b^2 + sigma_e^2). This invariance isn't quite so severe, since (sigma_b^2)/(sigma_b^2 + sigma_e^2) is bounded to (0,1). There is some leverage for inference here.
# For example, if the outcomes are correlated, e.g., rho*(sigma_b^2)/(sigma_b^2 + sigma_e^2) = 0.8. Then, we learn something about both rho, and sigma_e relative to sigma_b from the data. Namely,
# it is exceedingly unlikely that rho (the latent network dyadic reciprocity parameter) is less than 0.8, since even if sigma_e goes to zero, rho can't be smaller than 0.8 and yield the observed data. However, 
# if sigma_e is non-zero then rho might be larger than 0.8, with measurement error masking its true magnitude. Even so, we can also bound the measurement error, since rho is capped at 1 by theory. Rearranging 
# (sigma_b^2)/(sigma_b^2 + sigma_e^2) = 0.8,  to sigma_b^2 / 4 = sigma_e^2, shows that sigma_e^2 is at most about 25% as large as sigma_b^2 if the observations have correlation of 0.8.

# However, beyond these simple bounds, we cannot hope learn much from a single layer of network data. All combinations of rho in (0.8, 1.0) and sigma_e^2 in (0, sigma_b^2/4) are equally plausible absent prior information.  
# This invariance, however, is specifically what STRAND was designed to represent. From the observations alone, we cannot precisely conclude what rho is. We want a model that represents this posterior uncertainty.


# 3) We can recover parameters using STRAND, which would not be possible if our model were unidentifiable.
# Lets look at this in detail and compare STRAND to AMEN.

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
 N_id = 60
 B = -2.5
 sr_sigma = c(1.25, 1.25)  # Sender and receiver effect SDs in the true network 
 sr_rho = 0.6              # Sender and receiver correlation in the true network
 dr_sigma = 2.5            # Dyadic effect SD in the true network 
 dr_rho = 0.8              # Dyadic correlation in the true network 
 error_sigma = 2.5         # SD of dyadically uncorrelated measurement error/noise. (Noise in measurements is here of same scale as true dyadic effects).

A = simulate_srm_network(N_id = N_id, 
                         B=B, 
                         sr_sigma = sr_sigma, 
                         sr_rho = sr_rho,
                         dr_sigma = dr_sigma, 
                         dr_rho = dr_rho,
                         error_sigma = error_sigma,
                         outcome_mode="gaussian",
                         link_mode="identity"
                         )

Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
par(mfrow=c(1,2))
plot(Net, edge.arrow.size = 0.1, edge.curved = 0.3, vertex.label=NA, vertex.size = 5)
image(A$network)

######################################################## True variance partitions
# We use the Koster and Leckie definition of VPCs here
vcs = c(sr_sigma, dr_sigma, error_sigma)^2
vpcs = vcs/(vcs[1] + vcs[2] + vcs[3] + vcs[4])
d_e = sum(vcs[3:4]/(vcs[1] + vcs[2] + vcs[3] + vcs[4]))

True_Measures_Alt = c(sr_rho,                                            # Sender and receiver correlation in the true network
                      dr_rho,                                            # Dyadic correlation in the true network 
                      vpcs,                                              # Actual VPCs
                      d_e,                                               # Sum of dyadic and error VCPs
                      dr_rho*(dr_sigma^2 /(dr_sigma^2 + error_sigma^2))  # Dyadic correlation in the measured/observed responses 
                      )

################################################################################################ Fit AMEN
fit_SRM = ame(A$network, family = "nrm")
summary(fit_SRM)

## The variance parameters for nodes are reported as variances and covariances
ame_sender_var = fit_SRM$VC[,1]
ame_covariance = fit_SRM$VC[,2]
ame_target_var = fit_SRM$VC[,3]
ame_dyadic_cor = fit_SRM$VC[,4] 


ame_generalized_reciprocity = ame_covariance / (sqrt (ame_sender_var * ame_target_var))

## One can also calculate the Variance Partition Coefficients
ame_dyadic_var = fit_SRM$VC[,5] 
ame_total_var = ame_sender_var + ame_target_var + ame_dyadic_var 

ame_sender_VPC = ame_sender_var / ame_total_var
ame_target_VPC = ame_target_var / ame_total_var
ame_dyadic_VPC = ame_dyadic_var / ame_total_var

ame_dyadic_cor_meas = ame_dyadic_cor*(ame_dyadic_var /(ame_dyadic_var + 0^2))

AME_Measures_Alt = rbind(Q(ame_generalized_reciprocity), rep(NA,3), Q(ame_sender_VPC), Q(ame_target_VPC), rep(NA,3), rep(NA,3), Q(ame_dyadic_VPC), Q(ame_dyadic_cor_meas))

colnames(AME_Measures_Alt) = c("M","L","H")
rownames(AME_Measures_Alt) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")

################################################################################################ Fit STRAND
########################################################################### Organize data
# Create the STRAND data object
rownames(A$network) = colnames(A$network) = paste0("ID", 1:N_id)
outcome = list(trade = A$network)

dat2 = make_strand_data(
  self_report = outcome,
  individual_covariates = NULL,
  dyadic_covariates = NULL,
  outcome_mode = "gaussian",
  link_mode = "identity"
)

########################################## Fit STRAND with priors on measurement error which assume it can be of the same scale as other effects
fit1 =
  fit_block_plus_social_relations_model(
    data=dat2,
    block_regression = ~ 1,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ 1, 
    mode="mcmc",
    priors=make_priors(
     priors_to_change=list(
      "B_ingroup"=c(0.01, 5),          # Weak prior on intercept
      "sr_sigma"=c(1, 5), 
      "dr_sigma"=c(1, 5),
      "gaussian_error_priors"=c(1, 5)  # Weak prior that permits measurement error to be comparable in scale to sender/receiver/dyad effects
      )
    ),
    mcmc_parameters = list(
      chains = 1,
      seed=1,
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
STRAND_error_var = (res$samples$srm_model_samples$error_sd)^2
STRAND_total_var = STRAND_sender_var + STRAND_target_var + STRAND_dyadic_var + STRAND_error_var 

STRAND_sender_VPC = STRAND_sender_var / STRAND_total_var
STRAND_target_VPC = STRAND_target_var / STRAND_total_var
STRAND_dyadic_VPC = STRAND_dyadic_var / STRAND_total_var
STRAND_error_VPC = STRAND_error_var / STRAND_total_var

STRAND_dyadicerror_VPC = (STRAND_dyadic_var+STRAND_error_var) / STRAND_total_var

STRAND_dyadic_reciprocity_meas = STRAND_dyadic_reciprocity*(STRAND_dyadic_var /(STRAND_dyadic_var + STRAND_error_var))

STRAND_Measures_Weak_Alt = rbind(Q(STRAND_generalized_reciprocity), Q(STRAND_dyadic_reciprocity), Q(STRAND_sender_VPC), Q(STRAND_target_VPC), Q(STRAND_dyadic_VPC), Q(STRAND_error_VPC), Q(STRAND_dyadicerror_VPC), Q(STRAND_dyadic_reciprocity_meas))

colnames(STRAND_Measures_Weak_Alt) = c("M","L","H")
rownames(STRAND_Measures_Weak_Alt) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")

strand_VPCs(fit1, n_partitions = 3)
strand_VPCs(fit1, n_partitions = 3, include_reciprocity=TRUE)

########################################## Now fit STRAND with priors on measurement error which assume it is near zero
fit2 =
  fit_block_plus_social_relations_model(
    data=dat2,
    block_regression = ~ 1,
    focal_regression = ~ 1,
    target_regression = ~ 1,
    dyad_regression = ~ 1,
    mode="mcmc",
    priors=make_priors(
     priors_to_change=list(
      "B_ingroup"=c(0.01, 5),             # Weak prior on intercept
      "sr_sigma"=c(1, 5), 
      "dr_sigma"=c(1, 5),
      "gaussian_error_priors"=c(0, 0.1)   # Strong prior: normal(0, 0.1), that measurement error is near zero a priori
      )
    ),
    mcmc_parameters = list(
      chains = 1,
      seed=1,
      iter_warmup = 500 ,
      iter_sampling = 500 ,
      max_treedepth = 12 ,
      refresh = 1,
      adapt_delta = 0.96)
  )

res2 = summarize_strand_results(fit2)

### VPCs
STRAND_dyadic_reciprocity = res2$samples$srm_model_samples$dyadic_L[,2,1]
STRAND_generalized_reciprocity = res2$samples$srm_model_samples$focal_target_L[,2,1]

STRAND_sender_var = (res2$samples$srm_model_samples$focal_target_sd[,1])^2
STRAND_target_var = (res2$samples$srm_model_samples$focal_target_sd[,2])^2
STRAND_dyadic_var = (res2$samples$srm_model_samples$dyadic_sd)^2
STRAND_error_var = (res2$samples$srm_model_samples$error_sd)^2
STRAND_total_var = STRAND_sender_var + STRAND_target_var + STRAND_dyadic_var + STRAND_error_var 

STRAND_sender_VPC = STRAND_sender_var / STRAND_total_var
STRAND_target_VPC = STRAND_target_var / STRAND_total_var
STRAND_dyadic_VPC = STRAND_dyadic_var / STRAND_total_var
STRAND_error_VPC = STRAND_error_var / STRAND_total_var

STRAND_dyadicerror_VPC = (STRAND_dyadic_var+STRAND_error_var) / STRAND_total_var

STRAND_dyadic_reciprocity_meas = STRAND_dyadic_reciprocity*(STRAND_dyadic_var /(STRAND_dyadic_var + STRAND_error_var))

STRAND_Measures_Moderate_Alt = rbind(Q(STRAND_generalized_reciprocity), Q(STRAND_dyadic_reciprocity), Q(STRAND_sender_VPC), Q(STRAND_target_VPC), Q(STRAND_dyadic_VPC), Q(STRAND_error_VPC), Q(STRAND_dyadicerror_VPC), Q(STRAND_dyadic_reciprocity_meas))

colnames(STRAND_Measures_Moderate_Alt) = c("M","L","H")
rownames(STRAND_Measures_Moderate_Alt) = c("Gen. Recip.","Dyad. Recip.","Focal VPC","Target VPC","Dyadic VPC", "Error VPC", "DyadicError VPC", "Meas. Dyad. Recip.")

strand_VPCs(fit2, n_partitions = 3)
strand_VPCs(fit2, n_partitions = 3, include_reciprocity=TRUE)
  
###################################################################################################### Plot all
AME_Measures_Alt = as.data.frame(AME_Measures_Alt)
STRAND_Measures_Weak_Alt = as.data.frame(STRAND_Measures_Weak_Alt)
STRAND_Measures_Moderate_Alt = as.data.frame(STRAND_Measures_Moderate_Alt)
TRUE_Measures_Alt = data.frame(M=True_Measures_Alt,L=NA,H=NA)

AME_Measures_Alt$Package="AMEN"
STRAND_Measures_Weak_Alt$Package="STRAND (Weak prior)"
STRAND_Measures_Moderate_Alt$Package="STRAND (Strong prior)"
TRUE_Measures_Alt$Package="TRUE"

Measures_Alt = rbind(TRUE_Measures_Alt, 
                 AME_Measures_Alt,
                 STRAND_Measures_Weak_Alt, 
                 STRAND_Measures_Moderate_Alt
                 )

Measures_Alt$Variable = rep(rownames(AME_Measures_Alt),4)

Measures_Alt$Type = ifelse(Measures_Alt$Variable %in% c("Gen. Recip.","Dyad. Recip.","Meas. Dyad. Recip."), "Reciprocity", "VPC")

Measures_Alt$Method = rep(c("True Value","AMEN", "STRAND (Weak prior)","STRAND (Strong prior)"),each=8)

pal = c("#0D205C", "#507BAF", "#C5570E", "#DF913E", "#13667B", "#35ABC0", 
"#E39913", "#F3BB57", "#721712", "#AF5249", "#0B231A", "#4E7A65", 
"#431008", "#A56137")

Measures_Alt$Variable = factor(Measures_Alt$Variable)
Measures_Alt$Variable = factor(Measures_Alt$Variable, levels=levels(Measures_Alt$Variable)[rev(c(6,1,5,7,2,4,3,8))])

p2 = ggplot(Measures_Alt,aes(x=Variable,y=M,ymin=L,ymax=H,color=Method))+ 
     geom_linerange(size=1,aes(color=Method), position = position_dodge2(width=0.5,preserve = "total",padding = 1))+
     geom_point(size=2,aes(color=Method), position = position_dodge2(width=0.5,preserve = "total",padding = 1))+ facet_grid(Type~.,scales="free",space = "free_y")+
     labs(y="Estimated value", x="") + theme(strip.text.x = element_text(size=17,face="bold"), 
     strip.text.y = element_text(size=17,face="bold"),axis.text=element_text(size=17),axis.title=element_text(size=17,
     face="bold"))+theme(strip.text.y = element_text(angle = 360),legend.title=element_text(size=17), 
    legend.text=element_text(size=17))  + coord_flip() + theme(panel.spacing = unit(1, "lines")) +
     scale_color_manual(values=c("True Value"="black", "AMEN"="darkred", "STRAND (Weak prior)"=pal[7], "STRAND (Strong prior)"=pal[5]
                 )) +
          scale_linetype_manual(values=c("None"="solid", "ME"="dashed")) +
    theme(legend.text=element_text(size=17)) 

 p2

###############################################################################################################################

# Results:

# Note that STRAND estimates both dyadic reciprocity in the true network (Dyadic Recip. in plot p2 above) and dyadic reciprocity in the reports/measurements (Meas. Dyadic Recip.). In AMEN, and similar
# models, users are taught that the strength of dyadic reciprocity in measurements directly reflects the level of dyadic reciprocity social ties. This, however, is only true if measurement error is zero,
# otherwise measurement error leads to a directional bias in estimation of dyadic reciprocity. STRAND corrects for the bias, which is why we return two different versions of dyadic reciprocity (one uncorrected, 
# and one bias-corrrected). 

# True dyadic reciprocity was set at 0.8 in the generative model, but amen (and the STRAND model with priors that force it to behave like amen) estimate measurement-level dyadic reciprocity to be about 0.4. Why?
# Measurement error introduces correlation-zero noise into the dyadic observations, masking the true level of reciprocity in the real latent network. STRAND integrates over all causal paths of producing
# observations with a dyadic correlation of 0.4. This includes the extreme endpoint of a true dyadic reciprocity of 0.4, and no measurement error, as well as generative processes with higher levels of dyadic 
# reciprocity coupled with higher levels of measurement error. As such, STRAND's posteriors on dyadic_rho and dyadic_sigma are much wider, better reflecting the true level of uncertainty we have in the 
# parameters, when noise in measurement is possible.

# Because measurement error can mask dyadic reciprocity, we can only learn a lower bound on the magnitude of dyadic reciprocity from the data. STRAND represents this, while amen doesn't.

# Notes: 

# 1) One can always recover amen-style 3-way variance partitions from a STRAND model just by summing the dyadic and error VPCs: Or by using strand_VPCs(fit2, n_partitions = 3)

# 2) One can always recover the amen-style, measurement-biased, dyadic reciprocity estimate from a STRAND model just by multiplying dr_rho by: dr_sigma^2 / (dr_sigma^2 + error_sigma^2) 
#  or by running strand_VPCs(fit2, n_partitions = 3, include_reciprocity=TRUE, mode="adj")

# 3) If users have special reason to believe that measurement error is zero a priori, then its fine to use strong, amen-like, priors: e.g., gaussian_error_priors = c(0, 0.1)
# otherwise, consider using weaker priors like: gaussian_error_priors = c(1, 2.5), in order to get estimates of the posterior distribution on dr_rho that are less biased by measurement error, and better 
# reflect the substantial uncertainty we have in the value of dr_rho when only single layer data are supplied and measurements are potentially noisy.



