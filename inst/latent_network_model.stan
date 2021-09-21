functions{
    //# probability of observed variables relevant to ij dyad
    real prob_sgij( int[] sij, int[] sji, int tie, vector fpri, vector fprj, vector rtti, vector rttj, real theta) {
       vector[2] y;
      
    //# model likelihood
    if(tie==0){
            //# prob i says i helps j
        y[1] = bernoulli_logit_lpmf(sij[1] | fpri[1]);

            //# prob j says i helps j
            if(sji[1] == 0){
        y[2] = bernoulli_logit_lpmf(sji[2] | fprj[2]);
                           }
              else{
        y[2] = log_mix(theta, 
                      bernoulli_logit_lpmf(sji[2] | 1),
                      bernoulli_logit_lpmf(sji[2] | fprj[2])); 
                           }
      }

    else{ //#  if(tie==1){
        //# prob i says i helps j
        y[1] = bernoulli_logit_lpmf(sij[1] | rtti[1]);

        //# prob j says i helps j
        if(sji[1] == 0){
        y[2] = bernoulli_logit_lpmf(sji[2] | rttj[2]);
                      }
              else{
        y[2] = log_mix(theta, 
                      bernoulli_logit_lpmf(sji[2] | 1),
                      bernoulli_logit_lpmf(sji[2] | rttj[2])); 
        }
      }
        return(sum(y));
    }
}

data{
    int N_networktypes;                                               
    int N_id;                                                         
    int N_groups;                                                     
    int N_responses;        

    int N_params [6];                                          

    int group_ids[N_id];                                              
    int outcomes[N_id,N_id,N_responses];                              

    matrix[N_id, N_params[1]] focal_set;
    matrix[N_id, N_params[2]] target_set;
    matrix[N_id, N_params[3]] fpr_set;
    matrix[N_id, N_params[4]] rtt_set;
    matrix[N_id, N_params[5]] theta_set;

    real dyad_set[N_id, N_id, N_params[6]];

    matrix[22, 2] priors;

    int export_network;
}

transformed data{
 real S;
 real penalty;
 int N_per_group [N_groups];

 matrix[N_id, N_params[1]-1] focal_individual_predictors; 
 matrix[N_id, N_params[2]-1] target_individual_predictors; 
 matrix[N_id, N_params[3]-1] fpr_individual_predictors; 
 matrix[N_id, N_params[4]-1] rtt_individual_predictors; 
 matrix[N_id, N_params[5]-1] theta_individual_predictors; 

 real dyad_individual_predictors[N_id, N_id, N_params[6]-1];

 //# By group Ns 
 for(k in 1: N_groups){
  N_per_group[k] = 0;
  }

 for(i in 1:N_id){
  N_per_group[group_ids[i]] += 1;
  }


 //# Make penalty terms
 S = 0;
 for(i in 1:N_id){
   for(j in 1:N_id){
     if(i != j)
   S += (outcomes[i,j,1] + outcomes[j,i,2])/2.0;
 }}

  penalty = S^(1/priors[19,1]);

 //# Make pruned data
  
  if(N_params[1]>1){
  for(i in 2:N_params[1]){
  focal_individual_predictors[ , i-1] = focal_set[,i];  
   }}

  if(N_params[2]>1){
  for(i in 2:N_params[2]){
  target_individual_predictors[ , i-1] = target_set[,i];  
   }}

  if(N_params[3]>1){
  for(i in 2:N_params[3]){
  fpr_individual_predictors[ , i-1] = fpr_set[,i];  
   }}

  if(N_params[4]>1){
  for(i in 2:N_params[4]){
  rtt_individual_predictors[ , i-1] = rtt_set[,i];  
   }}

  if(N_params[5]>1){
  for(i in 2:N_params[5]){
  theta_individual_predictors[ , i-1] = theta_set[,i];  
   }}

  if(N_params[6]>1){
  for(i in 2:N_params[6]){
  dyad_individual_predictors[ , , i-1] = dyad_set[,,i];  
   }}
}

parameters{
    //######################################################## Observation model
    //# Measurement model
    vector<lower=0, upper=1>[N_networktypes] false_positive_rate;
    vector<lower=0, upper=1>[N_networktypes] recall_of_true_ties;
    real<lower=0, upper=1> theta_mean;

    vector<lower=0>[N_networktypes] fpr_sigma;
    vector<lower=0>[N_networktypes] rtt_sigma;
    real<lower=0> theta_sigma;

    vector[N_networktypes] fpr_raw[N_id];
    vector[N_networktypes] rtt_raw[N_id];
    real theta_raw[N_id];

    vector[N_params[3]-1] fpr_effects[N_networktypes];
    vector[N_params[4]-1] rtt_effects[N_networktypes]; 
    vector[N_params[5]-1] theta_effects;   

    //########################################################### Latent Netowrk
    //# SBM + SRM model
    matrix[N_groups, N_groups] B;

    vector<lower=0>[2] sr_sigma;  //# Variation of sender-receiver effects
    cholesky_factor_corr[2] sr_L;
    vector[2] sr_raw[N_id];

    real<lower=0> dr_sigma;     //# Variation of dyadic effects
    cholesky_factor_corr[2] dr_L;
    matrix[N_id, N_id] dr_raw;

    //# Effects of covariate
    vector[N_params[1]-1] focal_effects;
    vector[N_params[2]-1] target_effects;
    vector[N_params[6]-1] dyad_effects;  
}

model{
  vector[2] sr[N_id];
  matrix[N_id, N_id] dr;
  vector[N_networktypes] fpr[N_id];
  vector[N_networktypes] rtt[N_id];
  real theta[N_id];

  vector[2] scrap;

  matrix[N_id, N_id] p;
  matrix[N_id, N_id] mixed_p;

    //# Priors on effects of covariates
     focal_effects ~ normal(priors[12,1], priors[12,2]);
     target_effects ~ normal(priors[13,1], priors[13,2]);
     dyad_effects ~ normal(priors[14,1], priors[14,2]);
     theta_effects ~ normal(priors[9,1], priors[9,2]);

   for(k in 1:N_networktypes){
    fpr_effects[k] ~ normal(priors[7,1], priors[7,2]);
    rtt_effects[k] ~ normal(priors[8,1], priors[8,2]); 
    }
    
    //# Priors for measurement model
    false_positive_rate ~ beta(priors[1,1], priors[1,2]);
    recall_of_true_ties ~ beta(priors[2,1], priors[2,2]);
    theta_mean ~ beta(priors[3,1], priors[3,2]);

    fpr_sigma ~ exponential(priors[4,1]);
    rtt_sigma ~ exponential(priors[5,1]);
    theta_sigma ~ exponential(priors[6,1]);

    for(i in 1:N_id){
    fpr_raw[i] ~ normal(0,1);
    rtt_raw[i] ~ normal(0,1);
    theta_raw[i] ~ normal(0,1);
    }

    for(i in 1:N_id){
    vector[N_networktypes] fpr_terms;
    vector[N_networktypes] rtt_terms;

     for(k in 1:N_networktypes){
      fpr_terms[k] = dot_product(fpr_effects[k],  to_vector(fpr_individual_predictors[i])); 
      rtt_terms[k] = dot_product(rtt_effects[k],  to_vector(rtt_individual_predictors[i])); 
     }

    fpr[i] = logit(false_positive_rate) + fpr_sigma .* fpr_raw[i] + fpr_terms;
    rtt[i] = logit(recall_of_true_ties) + rtt_sigma .* rtt_raw[i] + rtt_terms;
    theta[i] = inv_logit(logit(theta_mean) + theta_sigma * theta_raw[i] + dot_product(theta_effects,  to_vector(theta_individual_predictors[i])));
    }    

    //# Sender-receiver priors for social relations model
    for(i in 1:N_id)
    sr_raw[i] ~ normal(0,1);

    sr_sigma ~ exponential(priors[15,1]);
    sr_L ~ lkj_corr_cholesky(priors[17,1]);

    for(i in 1:N_id){
     vector[2] sr_terms;

     sr_terms[1] = dot_product(focal_effects,  to_vector(focal_individual_predictors[i]));
     sr_terms[2] = dot_product(target_effects,  to_vector(target_individual_predictors[i]));  

     sr[i] = diag_pre_multiply(sr_sigma, sr_L) * sr_raw[i] + sr_terms;
     }

    //# Dyadic priors for social relations model
    to_vector(dr_raw) ~ normal(0,1);
    dr_sigma ~ exponential(priors[16,1]);
    dr_L ~ lkj_corr_cholesky(priors[18,1]);

    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
     scrap[1] = dr_raw[i,j];
     scrap[2] = dr_raw[j,i];
     scrap = rep_vector(dr_sigma, 2) .* (dr_L*scrap);
     dr[i,j] = scrap[1] + dot_product(dyad_effects,  to_vector(dyad_individual_predictors[i, j, ]));
     dr[j,i] = scrap[2] + dot_product(dyad_effects,  to_vector(dyad_individual_predictors[j, i, ]));
     }}

    for(i in 1:N_id){
     dr[i,i] = -99; //# ignore this :)
    }

    //# priors for B
    for ( i in 1:N_groups ){
        for ( j in 1:N_groups ) {
            if ( i==j ) {
                B[i,j] ~ normal(logit(priors[10,1]/sqrt(N_per_group[i])), priors[10,2]);   //# transfers more likely within groups
            } else {
                B[i,j] ~ normal(logit(priors[11,1]/sqrt(N_per_group[i]*0.5 + N_per_group[j]*0.5)), priors[11,2]); //# transfers less likely between groups
            }
        }}

    //# likelihood
    for ( i in 1:N_id ) {
        for ( j in 1:N_id ) {
            if ( i != j ) {
                vector[2] terms;

                //# consider each possible state of true tie and compute prob of data
                for(tie in 0:1) {
                  terms[tie+1] = prob_sgij(outcomes[i,j,], outcomes[j,i,], tie, fpr[i], fpr[j], rtt[i], rtt[j], theta[j]);
                }
                
      p[i,j] = inv_logit( B[group_ids[i], group_ids[j]] + sr[i,1] + sr[j,2] + dr[i,j] );  //# Model as a mixture distribution
      mixed_p[i,j] = log_mix( p[i,j] , terms[2] , terms[1] );  //# Model as a mixture distribution
      }
   }
  }

  for ( i in 1:N_id ) {
   p[i,i] = 0; //# we only defined the off diagonal above. here we set the diagonal to zero, in order to call sum below
   mixed_p[i,i] = 0;
   }

 target += sum(mixed_p);

 target += normal_lpdf(sum(p) | S, penalty);
 }

generated quantities{
    // compute posterior prob of each network tie
    matrix[N_id*export_network, N_id*export_network] p_tie_out;
 
    if(export_network==1){                
                vector[2] terms;
                int tie;
                vector[N_networktypes] fpr[N_id];
                vector[N_networktypes] rtt[N_id];
                real theta[N_id];
                vector[2] sr[N_id];
                matrix[N_id, N_id] dr;
                vector[2] scrap;
                matrix[N_id, N_id] p;
            
    for(i in 1:N_id){
     vector[2] sr_terms;

     sr_terms[1] = dot_product(focal_effects,  to_vector(focal_individual_predictors[i]));
     sr_terms[2] = dot_product(target_effects,  to_vector(target_individual_predictors[i]));  

     sr[i] = diag_pre_multiply(sr_sigma, sr_L) * sr_raw[i] + sr_terms;
     }

    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
     scrap[1] = dr_raw[i,j];
     scrap[2] = dr_raw[j,i];
     scrap = rep_vector(dr_sigma, 2) .* (dr_L*scrap);
     dr[i,j] = scrap[1] + dot_product(dyad_effects,  to_vector(dyad_individual_predictors[i, j, ]));
     dr[j,i] = scrap[2] + dot_product(dyad_effects,  to_vector(dyad_individual_predictors[j, i, ]));
    }}

    for(i in 1:N_id){
     dr[i,i] = -99; //# ignore this :)
    }

    for(i in 1:N_id){
    vector[N_networktypes] fpr_terms;
    vector[N_networktypes] rtt_terms;

     for(k in 1:N_networktypes){
      fpr_terms[k] = dot_product(fpr_effects[k],  to_vector(fpr_individual_predictors[i])); 
      rtt_terms[k] = dot_product(rtt_effects[k],  to_vector(rtt_individual_predictors[i])); 
     }

    fpr[i] = logit(false_positive_rate) + fpr_sigma .* fpr_raw[i] + fpr_terms;
    rtt[i] = logit(recall_of_true_ties) + rtt_sigma .* rtt_raw[i] + rtt_terms;
    theta[i] = inv_logit(logit(theta_mean) + theta_sigma * theta_raw[i] + dot_product(theta_effects,  to_vector(theta_individual_predictors[i])));
    }  

    for ( i in 1:N_id ) {
        for ( j in 1:N_id ) {
            if ( i != j ) {
                // consider each possible state of true tie and compute prob of data
                p[i,j] = inv_logit( B[group_ids[i], group_ids[j]] + sr[i,1] + sr[j,2] + dr[i,j]);

                tie = 0;
                terms[1] = 
                    log1m( p[i,j] ) + 
                    prob_sgij(outcomes[i,j,], outcomes[j,i,], tie, fpr[i], fpr[j], rtt[i], rtt[j], theta[j]);

                tie = 1;
                terms[2] = 
                    log( p[i,j] ) + 
                    prob_sgij(outcomes[i,j,], outcomes[j,i,], tie, fpr[i], fpr[j], rtt[i], rtt[j], theta[j]);

                p_tie_out[i,j] = exp(terms[2] - log_sum_exp( terms ));
            }
        }//j
    }//i

  for ( i in 1:N_id ) {
   p[i,i] = 0; 
   p_tie_out[i,i] = 0;
   }
 }
}



