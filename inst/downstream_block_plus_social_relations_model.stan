data{   
  //# Array dimension variables                                   
    int N_id;                                        //# Number of people                                                                                                   
    int N_params;                                    //# Number of focal, target, and dyadic predictors

  //# Predictor variables                                                                                                         
    matrix[N_id, N_params] focal_set;                //# Predictor variables    

  //# Nodal variables                                                                                                         
    matrix[N_id, 2] sr_mus;                          //# Predictor variables   
    matrix[N_id, 2] sr_sds;                          //# Predictor variables   
    array[2] real Z;

  //# Outcome and exposure data
    array[N_id] int outcomes;                        //# Outcome if int
    array[N_id] real outcomes_real;                  //# Outcome if real
    array[N_id] int exposure;                        //# Exposure for each outcome
    array[N_id] int mask;                            //# Mask for each outcome

  //# Accessory parameters 
    matrix[23, 2] priors;                            //# Priors in a matrix, see details in the make_priors() function
    int export_network;                              //# Controls export of predictions
    int outcome_mode;                                //# 
    int link_mode;                                   //# 
}

transformed data{
  //# Refactor to the first predictor slot, because it is unity
    matrix[N_id, N_params-1] focal_predictors;    //# Same as focal_set without first column

  //# Make pruned data for predictor variables, by dropping first column
    if(N_params>1){
     for(i in 2:N_params){
     focal_predictors[ , i-1] = focal_set[,i];  
     }}
}

parameters{
    //# Variation of sender-receiver effects
    matrix[N_id,2] sr_raw;

    //# Effects of covariate
    vector[N_params-1] focal_effects;

    //# Error in Gaussian, Beta, Neg Bin, and Gamma model
    real<lower=0> error_sigma;  

    //# Basic effects
    real alpha;
    vector[2] kappa;    
}

model{
  //# Local storage to make code more readable
    matrix[N_id,2] sr; //# Sender and receiver effects
    vector[N_id] Theta;
                 
    //# The first step, is to run measurement error model
     to_vector(sr_raw) ~ normal(0,1);

     for(i in 1:N_id){
      sr[i] = sr_mus[i] + (sr_sds[i] .* sr_raw[i]);
     }

    //# Priors on effects of covariates
     alpha ~ normal(0, 2.5);
     kappa ~ normal(priors[12,1], priors[12,2]);
     focal_effects ~ normal(priors[12,1], priors[12,2]);
     error_sigma ~ normal(priors[23,1], priors[23,2]);

    //# likelihood
    for(i in 1:N_id){
       if(mask[i]==0){

      //# Linear model
       Theta[i] = alpha + kappa[1]*sr[i,1]*Z[1] + kappa[2]*sr[i,2]*Z[2] + dot_product(focal_effects,  to_vector(focal_predictors[i]));   

      //# Bernoulli model of the outcomes
      if(outcome_mode==1){
        if(link_mode==1){
         outcomes[i] ~ bernoulli_logit(Theta[i]);  
        }
        if(link_mode==2){
         outcomes[i] ~ bernoulli(Phi(Theta[i]));  
        }
       }

      //# Binomial model of the outcomes
      if(outcome_mode==2){
        if(link_mode==1){
         outcomes[i] ~ binomial_logit(exposure[i], Theta[i]);  
        }
        if(link_mode==2){
         outcomes[i] ~ binomial(exposure[i], Phi(Theta[i])); 
        }
       }
      
      //# Poisson model of the outcomes
      if(outcome_mode==3){
       outcomes[i] ~ poisson_log(Theta[i]);  
       }

      //# Negative binomial model of the outcomes
      if(outcome_mode==4){
       outcomes[i] ~ neg_binomial(exp(Theta[i])*error_sigma, error_sigma);  
       }
      
      //# Gaussian model of the outcomes
      if(outcome_mode==5){
       outcomes_real[i] ~ normal(Theta[i], error_sigma);  //# Then model the outcomes
       }

      //# Beta model of the outcomes
      if(outcome_mode==6){
        if(link_mode==1){
         outcomes_real[i] ~ beta(inv_logit(Theta[i])*error_sigma, (1-inv_logit(Theta[i]))*error_sigma);  //# Then model the outcomes
        }
        if(link_mode==2){
         outcomes_real[i] ~ beta(Phi(Theta[i])*error_sigma, (1-Phi(Theta[i]))*error_sigma);  //# Then model the outcomes
        }
       }

      //# Gamma model of the outcomes
      if(outcome_mode==7){
       outcomes_real[i] ~ gamma(exp(Theta[i])*error_sigma, error_sigma);  
       }

     }}


 }

