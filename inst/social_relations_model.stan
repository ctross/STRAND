data{
    int N_networktypes;                                               
    int N_id;                                                                                                            
    int N_responses;        

    array[3] int N_params;                                          
                                             
    array[N_id,N_id,N_responses] int outcomes;  
    array[N_id,N_id,N_responses] int exposure; 
    array[N_id,N_id,N_responses] int mask;                                   

    matrix[N_id, N_params[1]] focal_set;
    matrix[N_id, N_params[2]] target_set;

    array[N_id, N_id, N_params[3]] real dyad_set;

    matrix [22, 2] priors;
    
    int export_network;
    int outcome_mode;
    int link_mode;                           
}

transformed data{
 matrix[N_id, N_params[1]-1] focal_individual_predictors; 
 matrix[N_id, N_params[2]-1] target_individual_predictors; 

 array[N_id, N_id, N_params[3]-1] real dyad_individual_predictors; 

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
    matrix[1,1] B;

    vector<lower=0>[2] sr_sigma;  //# Variation of sender-receiver effects
    cholesky_factor_corr[2] sr_L;
    array[N_id] vector[2] sr_raw;

    real<lower=0> dr_sigma;     //# Variation of dyadic effects
    cholesky_factor_corr[2] dr_L;
    matrix[N_id, N_id] dr_raw;

    //# Effects of covariate
    vector[N_params[1]-1] focal_effects;
    vector[N_params[2]-1] target_effects;
    vector[N_params[3]-1] dyad_effects;  
}

transformed parameters{
    matrix[2*N_responses, 2*N_responses] G_corr; 
    matrix[2*N_responses, 2*N_responses] D_corr; 

    G_corr = tcrossprod(sr_L); 
    D_corr = tcrossprod(dr_L);  
}

model{
  array[N_id] vector[2] sr;
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

    //# priors for 
    B[1,1] ~ normal(logit(priors[10,1]/sqrt(N_id)), priors[10,2]);


    //# likelihood
    for( i in 1:N_id ) {
     for( j in 1:N_id ) {
       if( i != j ) {
        if(mask[i,j,1]==0){
          //# Loops above run over relevant outcomes, below we split by outcome types and links

      //########################### Bernoulli    
      if(outcome_mode==1){
        if(link_mode==1){
         outcomes[i,j,1] ~ bernoulli_logit(B[1,1] + sr[i,1] + sr[j,2] + dr[i,j]);  //# Then model the outcomes
        }

        if(link_mode==2){
         outcomes[i,j,1] ~ bernoulli(Phi(B[1,1] + sr[i,1] + sr[j,2] + dr[i,j]));  //# Then model the outcomes
        }   
       }

      //########################### Binomial 
      if(outcome_mode==2){
        if(link_mode==1){
         outcomes[i,j,1] ~ binomial_logit(exposure[i,j,1], B[1,1] + sr[i,1] + sr[j,2] + dr[i,j]);  //# Then model the outcomes
        }

        if(link_mode==2){
         outcomes[i,j,1] ~ binomial(exposure[i,j,1], Phi(B[1,1] + sr[i,1] + sr[j,2] + dr[i,j]));  //# Then model the outcomes
        }
       }
      
      //########################### Poisson: link is always log, so skip the if-staement for now 
      if(outcome_mode==3){
      outcomes[i,j,1] ~ poisson_log(B[1,1] + sr[i,1] + sr[j,2] + dr[i,j]);  //# Then model the outcomes
       }

       }
      }
     }}


 }


generated quantities{
    //# compute posterior prob of each network tie
    matrix[N_id*export_network, N_id*export_network] p;
    array[N_id*export_network] vector[2*export_network] sr;
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
     dr[i,j] = scrap[1] + dot_product(dyad_effects,  to_vector(dyad_individual_predictors[i, j, ])) + B[1, 1];
     dr[j,i] = scrap[2] + dot_product(dyad_effects,  to_vector(dyad_individual_predictors[j, i, ])) + B[1, 1];
    }}

    for(i in 1:N_id){
     dr[i,i] = -99; //# ignore this :)
    }


    for ( i in 1:N_id ) {
        for ( j in 1:N_id ) {
            if ( i != j ) {
      // consider each possible state of true tie and compute prob of data
      if(outcome_mode==1){
       if(link_mode==1){
       p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j]);
       }

       if(link_mode==2){
       p[i,j] = Phi( sr[i,1] + sr[j,2] + dr[i,j]);
       }
      }

      if(outcome_mode==2){
       if(link_mode==1){
       p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j]);
       }

       if(link_mode==2){
       p[i,j] = Phi( sr[i,1] + sr[j,2] + dr[i,j]);
       }
     }

      if(outcome_mode==3){
       p[i,j] = exp(sr[i,1] + sr[j,2] + dr[i,j]);  
       }
            }
        }
    }

  for ( i in 1:N_id ) {
   p[i,i] = 0; 
   }
 }
}
