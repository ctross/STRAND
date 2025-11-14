####################################################################################################################
#
#   Bernoulli analysis - with between-level interaction  
#
####################################################################################################################

# STRAND models are usually written in a way to make predictors of node-level and dyad-level features stand out.
# That is, we have separate equations for focal_regression, target_regression, and dyad_regression. However,
# researchers sometime want to let an individual feature (like a focal's age) interact with a dyadic feature,
# and this is also possible in STRAND. This is because individual data can always be coded equivalently as dyadic.
# Just let some individual-level feature like Age, be written as: Focal_Age[,i] = Age, for a focal effect,
# or Target_Age[i,] = Age, for target effects. Then, either of these can be included in the dyad regression,
# interacting with other variables if desired.

############# Load libraries
library(STRAND)
library(ggplot2)

# Load package data
data(Colombia_Data)

N = nrow(Colombia_Data$Relatedness)
Focal_Age = matrix(0,nrow=N, ncol=N)
rownames(Focal_Age) = colnames(Focal_Age) = rownames(Colombia_Data$Individual)

# Create Focal_Age as a dyadic variable by cloning across columns
for(i in 1:N){
 Focal_Age[,i] = standardize_strand(Colombia_Data$Individual$Age, center = FALSE)      
}

############# Create the STRAND data object
outcome = list(Friends = Colombia_Data$Friends)

dyad = list(Relatedness = standardize_strand(Colombia_Data$Relatedness, center = FALSE), 
            Distance = standardize_strand(Colombia_Data$Distance),
            Focal_Age = Focal_Age
            )

dat = make_strand_data(outcome = outcome,
                       block_covariates = NULL, 
                       individual_covariates = NULL, 
                       dyadic_covariates = dyad,
                       outcome_mode = "bernoulli",
                       link_mode = "logit")


# Model
fit = fit_block_plus_social_relations_model(data=dat,
                                            block_regression = ~ 1,
                                            focal_regression = ~ 1,
                                            target_regression = ~ 1,
                                            dyad_regression = ~ Distance + Relatedness*Focal_Age,
                                            mode="mcmc",
                                            mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 500, iter_sampling = 500,
                                                                          max_treedepth = 11, adapt_delta = 0.98)
)

# Summaries
res = summarize_strand_results(fit)

############# Simulate predictions from posterior
  K = 100
  i_samps = res$samples$srm_model_samples$block_parameters[[1]] # intercept
  d_samps = res$samples$srm_model_samples$dyadic_coeffs         # dyadic effects (order of effects is the same as in the res object. 
                                                                # Here 1 is for distance, and 2, 3, and 4 are for kin, age, and kin*age effects)

  # Note that we fit the model on standardized predictors, 
  # so we need to some wrangling to predict back on the natural scale

  # Simulate for a predicted sweep of ages
  foc_age = seq(min(standardize_strand(Colombia_Data$Individual$Age, center = FALSE)),
                max(standardize_strand(Colombia_Data$Individual$Age, center = FALSE)),
                 length.out=K)

  # and 4 levels of kinship
  kin_set = c() 
  kin_set[1] = 0.000/sd(Colombia_Data$Relatedness)
  kin_set[2] = 0.125/sd(Colombia_Data$Relatedness)
  kin_set[3] = 0.250/sd(Colombia_Data$Relatedness)
  kin_set[4] = 0.500/sd(Colombia_Data$Relatedness)
  
################################## Now generate the predictions for age interacting with each level of kinship
  viz = array(NA, c(length(i_samps), K, 4))
  for(q in 1:length(i_samps)){
   for(k in 1:K){
     viz[q,k,1] = (i_samps[q] + d_samps[q,2]*kin_set[1] + d_samps[q,3]*foc_age[k] + d_samps[q,4]*foc_age[k]*kin_set[1])   
     viz[q,k,2] = (i_samps[q] + d_samps[q,2]*kin_set[2] + d_samps[q,3]*foc_age[k] + d_samps[q,4]*foc_age[k]*kin_set[2])    
     viz[q,k,3] = (i_samps[q] + d_samps[q,2]*kin_set[3] + d_samps[q,3]*foc_age[k] + d_samps[q,4]*foc_age[k]*kin_set[3])    
     viz[q,k,4] = (i_samps[q] + d_samps[q,2]*kin_set[4] + d_samps[q,3]*foc_age[k] + d_samps[q,4]*foc_age[k]*kin_set[4])                                       
        }
  }

################################## Get mean and credible intervals of predictions
m_viz = apply(viz, 2:3, mean)
l_viz = apply(viz, 2:3, HPDI)[1,,]
h_viz = apply(viz, 2:3, HPDI)[2,,]

eth_viz = fa_viz = m_viz
for(i in 1:4) fa_viz[,i] = foc_age*sd(Colombia_Data$Individual$Age)
for(i in 1:4) eth_viz[,i] = c("r=0", "r=0.125", "r=0.25", "r=0.5")[i]

# Store and plot the results
dat_viz = data.frame(Focal_Age=c(fa_viz), Kinship=c(eth_viz), M=c(m_viz), L=c(l_viz), H=c(h_viz))

p1 = ggplot() +
  geom_ribbon(data = dat_viz, aes(x = Focal_Age, y = M, ymin = L, ymax = H), alpha=0.5) +
  geom_line(data = dat_viz, aes(x = Focal_Age, y = M)) +
  facet_wrap(.~Kinship) +
  theme(legend.position="bottom")
p1 


