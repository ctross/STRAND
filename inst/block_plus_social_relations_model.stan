data{
    int N_networktypes;                                               
    int N_id;                                                         
    int N_groups;                                                     
    int N_responses;        

    int N_params [3];                                          

    int group_ids[N_id];                                              
    int outcomes[N_id,N_id,N_responses];                              

    matrix[N_id, N_params[1]] focal_set;
    matrix[N_id, N_params[2]] target_set;
    real dyad_set[N_id, N_id, N_params[3]];

    matrix[22, 2] priors;

    int export_network;
}

transformed data{
 matrix[N_id, N_params[1]-1] focal_individual_predictors; 
 matrix[N_id, N_params[2]-1] target_individual_predictors; 
 real dyad_individual_predictors[N_id, N_id, N_params[3]-1]; 

 int N_per_group [N_groups];

  //# By group Ns 
 for(k in 1: N_groups){
  N_per_group[k] = 0;
  }

 for(i in 1:N_id){
  N_per_group[group_ids[i]] += 1;
  }

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
  dyad_individual_predictors[ , , i-1] = dyad_set[,,i];  
   }}
}

parameters{
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
    vector[N_params[3]-1] dyad_effects;  
}

model{
  vector[2] sr[N_id];
  matrix[N_id, N_id] dr;

  vector[2] scrap;

    //# Priors on effects of covariates
     focal_effects ~ normal(priors[12,1], priors[12,2]);
     target_effects ~ normal(priors[13,1], priors[13,2]);
     dyad_effects ~ normal(priors[14,1], priors[14,2]);

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
     outcomes[i,j,1] ~ bernoulli_logit( B[group_ids[i], group_ids[j]] + sr[i,1] + sr[j,2] + dr[i,j] );  
       }
      }
     }


 }


generated quantities{
    //# compute posterior prob of each network tie
    matrix[N_id*export_network, N_id*export_network] p;
    vector[2*export_network] sr[N_id*export_network];
    matrix[N_id*export_network, N_id*export_network] dr;
 
    if(export_network==1){                
                vector[2] terms;
                int tie;
                vector[2] scrap;
            
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
     dr[i,j] = scrap[1] + dot_product(dyad_effects,  to_vector(dyad_individual_predictors[i, j, ])) + B[group_ids[i], group_ids[j]];
     dr[j,i] = scrap[2] + dot_product(dyad_effects,  to_vector(dyad_individual_predictors[j, i, ])) + B[group_ids[j], group_ids[i]];
    }}

    for(i in 1:N_id){
     dr[i,i] = -99; //# ignore this :)
    }


    for ( i in 1:N_id ) {
        for ( j in 1:N_id ) {
            if ( i != j ) {
                // consider each possible state of true tie and compute prob of data
                p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j]);
            }
        }//j
    }//i

  for ( i in 1:N_id ) {
   p[i,i] = 0; 
   }
 }
}



