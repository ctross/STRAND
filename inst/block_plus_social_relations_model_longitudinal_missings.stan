data{   
  //# Array dimension variables                                   
    int N_id;                                                //# Number of people                                                                                                   
    int N_responses;                                         //# Number of outcome networks
    array[3] int N_params;                                   //# Number of focal, target, and dyadic predictors

  //# Block predictor variables 
    int N_group_vars;                                        //# Number of block structure variables
    int max_N_groups;                                        //# Max number of group labels in any variable
    array[N_group_vars] int N_groups_per_var;                //# Number of group labels, per variable type
    array[N_id, N_group_vars, N_responses] int block_set;    //# Dataframe holding the group ID codes for each person (rows) for each variable type (cols) for each timepoint

  //# Focal, target, and dyadic predictor variables                                                                                                      
    array[N_id, N_params[1], N_responses] real focal_set;         //# Focal slash decider predictor variables for each timepoint  
    array[N_id, N_params[2], N_responses] real target_set;        //# Target slash alter predictor variables for each timepoint
    array[N_id, N_id, N_params[3], N_responses] real dyad_set;    //# Dyadic predictor variables for each timepoint

  //# Outcome and exposure data
    array[N_id,N_id,N_responses] int outcomes;       //# Outcome network of binary ties for each timepoint
    array[N_id,N_id,N_responses] int exposure;       //# Exposure for each outcome for each timepoint
    array[N_id,N_id,N_responses] int mask;           //# Censoring mask for each outcome for each timepoint

    int N_missing_focal_set;  
    int N_missing_target_set;  
    int N_missing_dyad_set;  

    array[N_missing_focal_set,3] int locations_missing_focal_set;  
    array[N_missing_target_set,3] int locations_missing_target_set;  
    array[N_missing_dyad_set,4] int locations_missing_dyad_set;  

    matrix[2, N_params[1]-1] focal_lims;  
    matrix[2, N_params[2]-1] target_lims;  
    matrix[2, N_params[3]-1] dyad_lims;  

  //# Accessory paramters 
    matrix[23, 2] priors;                            //# Priors in a matrix, see details in the make_priors() function
    int export_network;                              //# Controls export of predictions
    int outcome_mode;                                //# Are outcomes binomial
    int link_mode;                                   //# Link type
    real bandage_penalty;                            //# Stitching strength
    int random_effects_mode;                         //# Mode for random effects
    int coefficient_mode;                            //# Mode for other effects
}

transformed data{
  //# Refactor to the first predictor slot, becuase it is unity
    array[N_id, N_params[1]-1, N_responses] real focal_predictors;            //# Same as focal_set without first column
    array[N_id, N_params[2]-1, N_responses] real target_predictors;           //# Same as target_set without first column
    array[N_id, N_id, N_params[3]-1, N_responses] real dyad_predictors;       //# Same as dyad_set without first shelf

  //# Store some key indexes
    array[max_N_groups, N_group_vars, N_responses] int N_per_group;      //# Number of people in each block-type for each group variable
    array[N_group_vars+1] int block_indexes;                             //# The indexes of each block parameter when stored as a vector instead of ragged array
    int block_param_size;                                                //# Total number of block-level parameters

  //# Get size of parameters for block model
    block_param_size = 0;                             //# Start at zero
    block_indexes[1] = 0;                             //# Start at zero for first index 
    
    for(q in 1: N_group_vars){
     block_param_size += N_groups_per_var[q]*N_groups_per_var[q];                       //# Count up number of parameters in each K by K block matrix and add to total
     block_indexes[1+q] = N_groups_per_var[q]*N_groups_per_var[q] + block_indexes[q];   //# Create cummulative sum of block indices, by adding new sum to old sum
     }

  for(l in 1:N_responses){
  //# First fill with scrap
    for(q in 1: N_group_vars){
    for(k in 1: max_N_groups){
     N_per_group[k, q, l] = 0;   
     }}

  //# Now fill in real values
    for(q in 1: N_group_vars){
    for(i in 1:N_id){
     N_per_group[block_set[i,q,l],q,l] += 1;
     }}

  //# Make pruned data for predictor variables, by dropping first column
    if(N_params[1]>1){
     for(i in 2:N_params[1]){
     focal_predictors[ , i-1, l] = focal_set[ , i, l ];  
     }}

    if(N_params[2]>1){
     for(i in 2:N_params[2]){
     target_predictors[ , i-1, l] = target_set[ , i , l];  
     }}

    if(N_params[3]>1){
     for(i in 2:N_params[3]){
     dyad_predictors[ , , i-1, l] = dyad_set[ , , i, l];  
     }}
   }

    //# Missing data parameter limits
    vector[N_missing_focal_set] imp_focal_set_L;  
    vector[N_missing_target_set] imp_target_set_L;  
    vector[N_missing_dyad_set] imp_dyad_set_L;  

    vector[N_missing_focal_set] imp_focal_set_H;  
    vector[N_missing_target_set] imp_target_set_H;  
    vector[N_missing_dyad_set] imp_dyad_set_H;  

    //# Now map in missings
    for(q in 1:N_missing_focal_set){
     imp_focal_set_L[q] = focal_lims[1, locations_missing_focal_set[q, 2] - 1];
     imp_focal_set_H[q] = focal_lims[2, locations_missing_focal_set[q, 2] - 1];
    }

    for(q in 1:N_missing_target_set){
     imp_target_set_L[q] = target_lims[1, locations_missing_target_set[q, 2] - 1];
     imp_target_set_H[q] = target_lims[2, locations_missing_target_set[q, 2] - 1];
    }

    for(q in 1:N_missing_dyad_set){
     imp_dyad_set_L[q] = dyad_lims[1, locations_missing_dyad_set[q, 3] - 1];
     imp_dyad_set_H[q] = dyad_lims[2, locations_missing_dyad_set[q, 3] - 1];
    }
}

parameters{
    //# Block effects, stored as a vector to save space
    array[N_responses] vector[block_param_size] block_effects_raw;

    //# Effects of covariate
    array[N_responses] vector[N_params[1]-1] focal_effects_raw;
    array[N_responses] vector[N_params[2]-1] target_effects_raw;
    array[N_responses] vector[N_params[3]-1] dyad_effects_raw;   

    //# Variation of sender-receiver effects
    vector<lower=0>[2*N_responses] sr_sigma;            
    cholesky_factor_corr[2*N_responses] sr_L;  
    array[N_id] vector[2*N_responses] sr_raw;                

    //# Variation of dyadic effects
    vector<lower=0>[N_responses] dr_sigma;              
    cholesky_factor_corr[2*N_responses] dr_L;   
    array[N_responses] matrix[N_id, N_id] dr_raw; 

    //# Error in Gaussian model
    vector<lower=0>[N_responses] error_sigma;  

    //# Missing data parameters
    vector<lower=0, upper=1>[N_missing_focal_set] imp_focal_set;  
    vector<lower=0, upper=1>[N_missing_target_set] imp_target_set;  
    vector<lower=0, upper=1>[N_missing_dyad_set] imp_dyad_set;    
}

transformed parameters{
    //# Block effects, stored as a vector to save space
    array[N_responses] vector[block_param_size] block_effects;

    //# Effects of covariate
    array[N_responses] vector[N_params[1]-1] focal_effects;
    array[N_responses] vector[N_params[2]-1] target_effects;
    array[N_responses] vector[N_params[3]-1] dyad_effects;   

    matrix[2*N_responses, 2*N_responses] G_corr; 
    matrix[2*N_responses, 2*N_responses] D_corr; 

    G_corr = tcrossprod(sr_L); 
    D_corr = tcrossprod(dr_L);  

    //# If coefficient_mode = 1, for time-variant
    if(coefficient_mode == 1){
     for(m in 1:N_responses){
      block_effects[m] = block_effects_raw[m];
      focal_effects[m] = focal_effects_raw[m];
      target_effects[m] = target_effects_raw[m];
      dyad_effects[m] = dyad_effects_raw[m];
    }}

    //# If coefficient_mode = 2, for time-invariant
    if(coefficient_mode == 2){
     for(m in 1:N_responses){
      block_effects[m] = block_effects_raw[1];
      focal_effects[m] = focal_effects_raw[1];
      target_effects[m] = target_effects_raw[1];
      dyad_effects[m] = dyad_effects_raw[1];
    }}
}

model{
  //# Local storage to make code more readable
    array[N_id] vector[2] sr;                                  //# Sender and receiver effects
    array[N_id] vector[2*N_responses] sr_multi;                //# Sender and receiver effects long form
    matrix[N_id, N_id] dr;                                     //# Dyadic effects
    array[N_responses] matrix[N_id, N_id] dr_multi;            //# Dyadic effects long form
    array[N_group_vars] matrix[max_N_groups, max_N_groups] B;  //# Block effects, in array form
    vector[N_group_vars] br;                                   //# Sum of block effects per dyad    
    vector[2*N_responses] scrap;                               //# Local storage  

    array[N_id, N_params[1]-1, N_responses] real focal_predictors_mixed = focal_predictors;
    array[N_id, N_params[2]-1, N_responses] real target_predictors_mixed = target_predictors;
    array[N_id, N_id, N_params[3]-1, N_responses] real dyad_predictors_mixed = dyad_predictors;

    //# Priors on imputed values
     imp_focal_set ~ uniform(0, 1);  
     imp_target_set ~ uniform(0, 1);   
     imp_dyad_set ~ uniform(0, 1);  

    //# Now map in missings
    for(q in 1:N_missing_focal_set){
     focal_predictors_mixed[locations_missing_focal_set[q,1], locations_missing_focal_set[q,2]-1, locations_missing_focal_set[q,3]] = imp_focal_set_L[q] + (imp_focal_set_H[q] - imp_focal_set_L[q]) * imp_focal_set[q];
    }
    
    for(q in 1:N_missing_target_set){
     target_predictors_mixed[locations_missing_target_set[q,1], locations_missing_target_set[q,2]-1, locations_missing_target_set[q,3]] = imp_target_set_L[q] + (imp_target_set_H[q] - imp_target_set_L[q]) * imp_target_set[q];
    }

    for(q in 1:N_missing_dyad_set){
     dyad_predictors_mixed[locations_missing_dyad_set[q,1], locations_missing_dyad_set[q,2], locations_missing_dyad_set[q,3]-1, locations_missing_dyad_set[q,4]] = imp_dyad_set_L[q] + (imp_dyad_set_H[q] - imp_dyad_set_L[q]) * imp_dyad_set[q];
    }                     
    
    //# Stitch together the dyadic matrix
    for(m in 1:(N_responses-1)){
    for(n in (m+1):N_responses){
     target += normal_lpdf(D_corr[m+N_responses, n+N_responses] | D_corr[m, n], bandage_penalty);
     target += normal_lpdf(D_corr[m, n+N_responses]   | D_corr[n, m+N_responses], bandage_penalty);
    }}

    //# Add stiches to the dyadic matrix for longitudinal symmetries
    //# If random_effects_mode = 2, for time-invariant
    if(random_effects_mode==2){
     for(k in 1:(N_responses-1)){
     for(m in 1:(N_responses-k)){
      target += normal_lpdf(D_corr[m, m+k] | D_corr[1, k+1],   bandage_penalty);
      target += normal_lpdf(D_corr[m, N_responses+m+k] | D_corr[1, N_responses+k+1], bandage_penalty);
      }
     }

     for(m in 1:N_responses){
      target += normal_lpdf(D_corr[m, m+N_responses] | D_corr[1, N_responses+1], bandage_penalty);
      }
    }

    //# Add stiches to the generalized matrix for longitudinal symmetries
    //# If random_effects_mode = 2, for time-invariant
    if(random_effects_mode==2){
     for(k in 1:(N_responses-1)){
     for(m in 1:(N_responses-k)){
      target += normal_lpdf(G_corr[m, m+k] | G_corr[1, k+1], bandage_penalty);
      target += normal_lpdf(G_corr[m, N_responses+m+k] | G_corr[1, N_responses+k+1], bandage_penalty);
      target += normal_lpdf(G_corr[N_responses + m, N_responses + m+k] | G_corr[N_responses + 1, N_responses + k+1], bandage_penalty);
      target += normal_lpdf(G_corr[m + k, N_responses + m] | G_corr[k + 1, N_responses + 1], bandage_penalty);
      }}

     for(m in 1:N_responses){
      target += normal_lpdf(G_corr[m, m+N_responses] | G_corr[1, N_responses+1], bandage_penalty);
     }
    }

    //# Sender-receiver priors for social relations model
    for(i in 1:N_id)
    sr_raw[i] ~ normal(0,1);
    sr_sigma ~ normal(priors[15,1], priors[15,2]);
    sr_L ~ lkj_corr_cholesky(priors[17,1]);

    //# Dyadic priors for social relations model
    for(l in 1:N_responses)
    to_vector(dr_raw[l]) ~ normal(0,1);
    dr_sigma ~ normal(priors[16,1], priors[16,2]);
    dr_L ~ lkj_corr_cholesky(priors[18,1]);
    
    for(i in 1:N_id){
      sr_multi[i] = diag_pre_multiply(sr_sigma, sr_L) * sr_raw[i];
     }

    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){

         for(l in 1:N_responses){
     scrap[l] = dr_raw[l,i,j];
     scrap[l+N_responses] = dr_raw[l,j,i];
            }

     scrap = append_row(dr_sigma, dr_sigma) .* (dr_L*scrap);

          for(l in 1:N_responses){
     dr_multi[l,i,j] = scrap[l];
     dr_multi[l,j,i] = scrap[l+N_responses];
            }

     }}

    // This is a loop over years
    for(l in 1:N_responses){

    //# The first step, is to transform the vector of block effects into a list of matrices
     for(q in 1:N_group_vars){
      B[q,1:N_groups_per_var[q], 1:N_groups_per_var[q]] = to_matrix(block_effects[l,(block_indexes[q]+1):(block_indexes[q+1])], N_groups_per_var[q], N_groups_per_var[q]);
     }

    //# Then put priors on B, which scale loosely with the block size
    for ( q in 1:N_group_vars ){
    for ( i in 1:N_groups_per_var[q] ){
        for ( j in 1:N_groups_per_var[q] ) {
            if ( i==j ) {
                B[q,i,j] ~ normal(logit(priors[10,1]/sqrt(N_per_group[i,q,l])), priors[10,2]);                              //# transfers more likely within groups
            } else {
                B[q,i,j] ~ normal(logit(priors[11,1]/sqrt(N_per_group[i,q,l]*0.5 + N_per_group[j,q,l]*0.5)), priors[11,2]); //# transfers less likely between groups
            }
        }}
    }

    //# Priors on effects of covariates
     focal_effects[l] ~ normal(priors[12,1], priors[12,2]);
     target_effects[l] ~ normal(priors[13,1], priors[13,2]);
     dyad_effects[l] ~ normal(priors[14,1], priors[14,2]);

     error_sigma[l] ~ normal(priors[23,1], priors[23,2]);

     for(i in 1:N_id){                                                                
     sr[i,1] = sr_multi[i, l] + dot_product(focal_effects[l],  to_vector(focal_predictors_mixed[i, ,l]));
     sr[i,2] = sr_multi[i, l + N_responses] + dot_product(target_effects[l],  to_vector(target_predictors_mixed[i, ,l]));  
     }

    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
     dr[i,j] = dr_multi[l,i,j] + dot_product(dyad_effects[l],  to_vector(dyad_predictors_mixed[i, j, , l]));   
     dr[j,i] = dr_multi[l,j,i] + dot_product(dyad_effects[l],  to_vector(dyad_predictors_mixed[j, i, , l]));
     }}

    for(i in 1:N_id){
     dr[i,i] = -99; //# ignore this :)
    }

    //# likelihood
    for(i in 1:N_id){
     for(j in 1:N_id){
       if(i != j){
         if(mask[i,j,l]==0){
        for(q in 1:N_group_vars){
          br[q] = B[q,block_set[i,q,l], block_set[j,q,l]]; //# Extract all of the block components for this dyad
         }

      if(outcome_mode==1){
        if(link_mode==1){
         outcomes[i,j,l] ~ bernoulli_logit(sum(br) + sr[i,1] + sr[j,2] + dr[i,j]);  //# Then model the outcomes
        }
        if(link_mode==2){
         outcomes[i,j,l] ~ bernoulli(Phi(sum(br) + sr[i,1] + sr[j,2] + dr[i,j]));  //# Then model the outcomes
        }
       }

      if(outcome_mode==2){
        if(link_mode==1){
         outcomes[i,j,l] ~ binomial_logit(exposure[i,j,l], sum(br) + sr[i,1] + sr[j,2] + dr[i,j]);  //# Then model the outcomes
        }
        if(link_mode==2){
         outcomes[i,j,l] ~ binomial(exposure[i,j,l], Phi(sum(br) + sr[i,1] + sr[j,2] + dr[i,j]));  //# Then model the outcomes
        }
       }

      if(outcome_mode==3){
         outcomes[i,j,l] ~ poisson_log(sum(br) + sr[i,1] + sr[j,2] + dr[i,j]);  //# Then model the outcomes
       }

      if(outcome_mode==4){
         outcomes[i,j,l] ~ normal(sum(br) + sr[i,1] + sr[j,2] + dr[i,j], error_sigma[l]);  //# Then model the outcomes
       }

       }
      }
     }
    }
    }
 }
