functions{
  //# Custom probability mass function: A Logit Model for Bivariate Binary Responses. Purhadi and Fathurahman.
     real bivariate_bernoulli_lpmf(int [] Y, real p, real q, real rho) {
       real a = 1 + (p + q)*(rho-1);
       real b = 4*rho*(1-rho)*p*q;
       real y11 = (a - sqrt(a^2 + b))*0.5*(rho-1)^(-1);
       //# real y01 = q - y11;
       //# real y10 = p - y11;
       //# real y00 = 1 - p - q + y11;

       if(Y[1] == 0 && Y[2] == 0)     
        return log(1 - p - q + y11);

       if(Y[1] == 1 && Y[2] == 1)  
        return log(y11);

       if(Y[1] == 0 && Y[2] == 1)  
        return log(q - y11);

       if(Y[1] == 1 && Y[2] == 0)  
        return log(p - y11);

        return log(0);
   }
}

data{   
  //# Array dimension variables                                   
    int N_id;                                  //# Number of people
    int N_hh;                                  //# Number of households                                                                                                   
    int N_responses;                           //# Number of outcome networks
    int N_params[3];                           //# Number of focal, target, and dyadic predictors
    int N_params_hh[3];                        //# Number of HH focal, target, and dyadic predictors

  //# Block predictor variables 
    int N_group_vars;                          //# Number of block structure variables
    int max_N_groups;                          //# Max number of group labels in any variable
    int N_groups_per_var[N_group_vars];        //# Number of group labels, per variable type
    int block_set[N_id, N_group_vars];         //# Dataframe holding the group ID codes for each person (rows) for each variable type (cols)

  //# Focal, target, and dyadic predictor variables                                                                                                      
    matrix[N_id, N_params[1]] focal_set;       //# Focal slash decider predictor variables    
    matrix[N_id, N_params[2]] target_set;      //# Target slash alter predictor variables
    real dyad_set[N_id, N_id, N_params[3]];    //# Dyadic predictor variables

  //# HH Focal, target, and dyadic predictor variables                                                                                                      
    matrix[N_hh, N_params_hh[1]] hh_focal_set;       //# Focal slash decider predictor variables    
    matrix[N_hh, N_params_hh[2]] hh_target_set;      //# Target slash alter predictor variables
    real hh_dyad_set[N_hh, N_hh, N_params_hh[3]];    //# Dyadic predictor variables
    int HH[N_id];                                    //# Househould ID

  //# Outcome and exposure data
    int outcomes[N_id,N_id,N_responses];             //# Outcome network of binary ties
    int exposure[N_id,N_id,N_responses];             //# Exposure for each outcome

  //# Accessory paramters 
    matrix[22, 2] priors;                            //# Priors in a matrix, see details in the make_priors() function
    int export_network;                              //# Controls export of predictions
    int outcome_mode;                                //# Are outcomes binomial
}

transformed data{
  //# Refactor to the first predictor slot, becuase it is unity
    matrix[N_id, N_params[1]-1] focal_predictors;     //# Same as focal_set without first column
    matrix[N_id, N_params[2]-1] target_predictors;    //# Same as target_set without first column
    real dyad_predictors[N_id, N_id, N_params[3]-1];  //# Same as dyad_set without first shelf

    //# Refactor to the first predictor slot, becuase it is unity
    matrix[N_hh, N_params_hh[1]-1] hh_focal_predictors;     //# Same as focal_set without first column
    matrix[N_hh, N_params_hh[2]-1] hh_target_predictors;    //# Same as target_set without first column
    real hh_dyad_predictors[N_hh, N_hh, N_params_hh[3]-1];  //# Same as dyad_set without first shelf

  //# Store some key indexes
    int N_per_group [max_N_groups, N_group_vars];     //# Number of people in each block-type for each group variable
    int block_indexes[N_group_vars+1];                //# The indexes of each block parameter when stored as a vector instead of ragged array
    int block_param_size;                             //# Total number of block-level parameters

  //# Get size of parameters for block model
    block_param_size = 0;                             //# Start at zero
    block_indexes[1] = 0;                             //# Start at zero for first index 
    
    for(q in 1: N_group_vars){
     block_param_size += N_groups_per_var[q]*N_groups_per_var[q];                      //# Count up number of parameters in each K by K block matrix and add to total
     block_indexes[1+q] = N_groups_per_var[q]*N_groups_per_var[q] + block_indexes[q];  //# Create cummulative sum of block indices, by adding new sum to old sum
     }

  //# First fill with scrap
    for(q in 1: N_group_vars){
    for(k in 1: max_N_groups){
     N_per_group[k, q] = 0;   
     }}

  //# Now fill in real values
    for(q in 1: N_group_vars){
    for(i in 1:N_id){
     N_per_group[block_set[i,q],q] += 1;
     }}

  //# Make pruned data for predictor variables, by dropping first column
    if(N_params[1]>1){
     for(i in 2:N_params[1]){
     focal_predictors[ , i-1] = focal_set[,i];  
     }}

    if(N_params[2]>1){
     for(i in 2:N_params[2]){
     target_predictors[ , i-1] = target_set[,i];  
     }}

    if(N_params[3]>1){
     for(i in 2:N_params[3]){
     dyad_predictors[ , , i-1] = dyad_set[,,i];  
     }}

  //# Make pruned data for predictor variables, by dropping first column
    if(N_params_hh[1]>1){
     for(i in 2:N_params_hh[1]){
     hh_focal_predictors[ , i-1] = hh_focal_set[,i];  
     }}

    if(N_params_hh[2]>1){
     for(i in 2:N_params_hh[2]){
     hh_target_predictors[ , i-1] = hh_target_set[,i];  
     }}

    if(N_params_hh[3]>1){
     for(i in 2:N_params_hh[3]){
     hh_dyad_predictors[ , , i-1] = hh_dyad_set[,,i];  
     }}
}

parameters{
    //# Block effects, stored as a vector to save space
    vector[block_param_size] block_effects;

    //# Variation of sender-receiver effects
    vector<lower=0>[2] sr_sigma;  
    cholesky_factor_corr[2] sr_L;
    vector[2] sr_raw[N_id];

    //# Dyadic reciprocity term on log odds scale
    real<lower=0> dr_sigma;       

    //# Effects of covariate
    vector[N_params[1]-1] focal_effects;
    vector[N_params[2]-1] target_effects;
    vector[N_params[3]-1] dyad_effects;   

    //# Variation of HH sender-receiver effects
    vector<lower=0>[2] hh_sr_sigma;  
    cholesky_factor_corr[2] hh_sr_L;
    vector[2] hh_sr_raw[N_hh];

    //# Variation of HH dyadic effects
    real hh_dr_within_hh_offset;    
    real<lower=0> hh_dr_sigma;       
    cholesky_factor_corr[2] hh_dr_L;
    matrix[N_hh, N_hh] hh_dr_raw; 

    //# Effects of HH covariate
    vector[N_params_hh[1]-1] hh_focal_effects;
    vector[N_params_hh[2]-1] hh_target_effects;
    vector[N_params_hh[3]-1] hh_dyad_effects;   
}

model{
  //# Local storage to make code more readable
    vector[2] sr[N_id];                                   //# Sender and receiver effects
    vector[2] hh_sr[N_hh];                                //# HH Sender and receiver effects
    matrix[N_id, N_id] dr;                                //# Dyadic effects
    matrix[N_hh, N_hh] hh_dr;                             //# HH Dyadic effects
    matrix[max_N_groups, max_N_groups] B [N_group_vars];  //# Block effects, in array form
    vector[N_group_vars] br_ij;                           //# Sum of block effects per dyad    
    vector[N_group_vars] br_ji;                           //# Sum of block effects per dyad    
    vector[2] scrap;                                      //# Local storage  
    int outcome_data [2];                                 //# Local storage  
    real pred_data [2];                                   //# Local storage                
    
    //# The first step, is to transform the vector of block effects into a list of matrices
    for(q in 1:N_group_vars){
      B[q,1:N_groups_per_var[q], 1:N_groups_per_var[q]] = to_matrix(block_effects[(block_indexes[q]+1):(block_indexes[q+1])], N_groups_per_var[q], N_groups_per_var[q]);
    }

    //# Then put priors on B, which scale loosely with the block size
    for ( q in 1:N_group_vars ){
    for ( i in 1:N_groups_per_var[q] ){
        for ( j in 1:N_groups_per_var[q] ) {
            if ( i==j ) {
                B[q,i,j] ~ normal(logit(priors[10,1]/sqrt(N_per_group[i,q])), priors[10,2]);   //# transfers more likely within groups
            } else {
                B[q,i,j] ~ normal(logit(priors[11,1]/sqrt(N_per_group[i,q]*0.5 + N_per_group[j,q]*0.5)), priors[11,2]); //# transfers less likely between groups
            }
        }}
    }

    //# Priors on effects of covariates
     focal_effects ~ normal(priors[12,1], priors[12,2]);
     target_effects ~ normal(priors[13,1], priors[13,2]);
     dyad_effects ~ normal(priors[14,1], priors[14,2]);

    //# Priors on effects of hh covariates
     hh_focal_effects ~ normal(priors[12,1], priors[12,2]);
     hh_target_effects ~ normal(priors[13,1], priors[13,2]);
     hh_dyad_effects ~ normal(priors[14,1], priors[14,2]);

    //# Sender-receiver priors for social relations model
    for(i in 1:N_id)
    sr_raw[i] ~ normal(0,1);
    sr_sigma ~ exponential(priors[15,1]);
    sr_L ~ lkj_corr_cholesky(priors[17,1]);

    for(i in 1:N_id){
     vector[2] sr_terms;

     sr_terms[1] = dot_product(focal_effects,  to_vector(focal_predictors[i]));
     sr_terms[2] = dot_product(target_effects,  to_vector(target_predictors[i]));  

     sr[i] = diag_pre_multiply(sr_sigma, sr_L) * sr_raw[i] + sr_terms;
     }

    //# HH Sender-receiver priors for social relations model
    for(i in 1:N_hh)
    hh_sr_raw[i] ~ normal(0,1);
    hh_sr_sigma ~ exponential(priors[15,1]);
    hh_sr_L ~ lkj_corr_cholesky(priors[17,1]);
    hh_dr_within_hh_offset ~ normal(0, 1);

    for(i in 1:N_hh){
     vector[2] hh_sr_terms;

     hh_sr_terms[1] = dot_product(hh_focal_effects,  to_vector(hh_focal_predictors[i]));
     hh_sr_terms[2] = dot_product(hh_target_effects,  to_vector(hh_target_predictors[i]));  

     hh_sr[i] = diag_pre_multiply(hh_sr_sigma, hh_sr_L) * hh_sr_raw[i] + hh_sr_terms;
     }

    //# Dyadic priors for social relations model
    dr_sigma ~ exponential(priors[16,1]);

    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
     dr[i,j] = dot_product(dyad_effects,  to_vector(dyad_predictors[i, j, ]));
     dr[j,i] = dot_product(dyad_effects,  to_vector(dyad_predictors[j, i, ]));
     }}

    for(i in 1:N_id){
     dr[i,i] = -999; 
    }

    //# HH Dyadic priors for social relations model
    to_vector(hh_dr_raw) ~ normal(0,1);
    hh_dr_sigma ~ exponential(priors[16,1]);
    hh_dr_L ~ lkj_corr_cholesky(priors[18,1]);

    for(i in 1:(N_hh-1)){
    for(j in (i+1):N_hh){
     scrap[1] = hh_dr_raw[i,j];
     scrap[2] = hh_dr_raw[j,i];
     scrap = rep_vector(hh_dr_sigma, 2) .* (hh_dr_L*scrap);
     hh_dr[i,j] = scrap[1] + dot_product(hh_dyad_effects,  to_vector(hh_dyad_predictors[i, j, ]));
     hh_dr[j,i] = scrap[2] + dot_product(hh_dyad_effects,  to_vector(hh_dyad_predictors[j, i, ]));
     }}

    for(i in 1:N_hh){
     hh_dr[i,i] = hh_dr_within_hh_offset + hh_dr_sigma*hh_dr_raw[i,i];
    }

    //# likelihood
    for(i in 1:(N_id-1)){
     for(j in (i+1):N_id){

         for(q in 1:N_group_vars){
          br_ij[q] = B[q,block_set[i,q], block_set[j,q]]; //# Extract all of the block components for this dyad
         }

         for(q in 1:N_group_vars){
          br_ji[q] = B[q,block_set[j,q], block_set[i,q]]; //# Extract all of the block components for this dyad
         }

        outcome_data[1] = outcomes[i,j,1];
        outcome_data[2] = outcomes[j,i,1];

        pred_data[1] = inv_logit(sum(br_ij) + sr[i,1] + sr[j,2] + dr[i,j] + hh_sr[HH[i],1] + hh_sr[HH[j],2] + hh_dr[HH[i],HH[j]]);  //# Then model the outcomes
        pred_data[2] = inv_logit(sum(br_ji) + sr[j,1] + sr[i,2] + dr[j,i] + hh_sr[HH[j],1] + hh_sr[HH[i],2] + hh_dr[HH[j],HH[i]]);  //# Then model the outcomes
       
        outcome_data ~ bivariate_bernoulli(pred_data[1], pred_data[2], dr_sigma);

      }}

 }


generated quantities{
    //# compute posterior prob of each network tie
    matrix[N_id*export_network, N_id*export_network] p;
    vector[2*export_network] sr[N_id*export_network];
    vector[2*export_network] hh_sr[N_hh*export_network];
    matrix[N_id*export_network, N_id*export_network] dr;
    matrix[N_hh*export_network, N_hh*export_network] hh_dr;
 
    if(export_network==1){                
     vector[2] terms;
     int tie;
     matrix[max_N_groups, max_N_groups] B[N_group_vars];

    for(i in 1:N_group_vars){
     B[i,1:N_groups_per_var[i], 1:N_groups_per_var[i]] = to_matrix(block_effects[(block_indexes[i]+1):(block_indexes[i+1])], N_groups_per_var[i], N_groups_per_var[i]);
    }
            
    for(i in 1:N_id){
     vector[2] sr_terms;

     sr_terms[1] = dot_product(focal_effects,  to_vector(focal_predictors[i]));
     sr_terms[2] = dot_product(target_effects,  to_vector(target_predictors[i]));  

     sr[i] = diag_pre_multiply(sr_sigma, sr_L) * sr_raw[i] + sr_terms;
     }

    for(i in 1:N_hh){
     vector[2] hh_sr_terms;

     hh_sr_terms[1] = dot_product(hh_focal_effects,  to_vector(hh_focal_predictors[i]));
     hh_sr_terms[2] = dot_product(hh_target_effects,  to_vector(hh_target_predictors[i]));  

     hh_sr[i] = diag_pre_multiply(hh_sr_sigma, hh_sr_L) * hh_sr_raw[i] + hh_sr_terms;
     }

    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
      vector[2] scrap;
      vector[N_group_vars] br1;
      vector[N_group_vars] br2;

       for(q in 1:N_group_vars){
        br1[q] = B[q,block_set[i,q], block_set[j,q]];
        br2[q] = B[q,block_set[j,q], block_set[i,q]];
         }

       dr[i,j] = dot_product(dyad_effects,  to_vector(dyad_predictors[i, j, ])) + sum(br1);
       dr[j,i] = dot_product(dyad_effects,  to_vector(dyad_predictors[j, i, ])) + sum(br2);
    }}

    for(i in 1:N_id){
      dr[i,i] = 0; 
    }

    for(i in 1:(N_hh-1)){
    for(j in (i+1):N_hh){
      vector[2] scrap;

      scrap[1] = hh_dr_raw[i,j];
      scrap[2] = hh_dr_raw[j,i];
      scrap = rep_vector(hh_dr_sigma, 2) .* (hh_dr_L*scrap);

      hh_dr[i,j] = scrap[1] + dot_product(hh_dyad_effects,  to_vector(hh_dyad_predictors[i, j, ]));
      hh_dr[j,i] = scrap[2] + dot_product(hh_dyad_effects,  to_vector(hh_dyad_predictors[j, i, ]));
    }}

    for(i in 1:N_hh){
      hh_dr[i,i] = hh_dr_within_hh_offset + hh_dr_sigma*hh_dr_raw[i,i];
    }


    for ( i in 1:N_id ) {
        for ( j in 1:N_id ) {
            if ( i != j ) {
       p[i,j] = inv_logit( sr[i,1] + sr[j,2] + dr[i,j] + hh_sr[HH[i],1] + hh_sr[HH[j],2] + hh_dr[HH[i],HH[j]]);
            }
        }
    }

  for ( i in 1:N_id ) {
   p[i,i] = 0; 
   }
 }
}



