functions{
    //# Shifted and scaled inv_logit link
    real inv_logit2(real x) {
    return((inv_logit(x) - 0.5)*2);
    }
}

data{   
  //# Array dimension variables                                   
    int N_id;                                                     //# Number of people                                                                                                   
    int N_responses;                                              //# Number of network time-points
    array[4] int N_params;                                        //# Number of focal, target, and dyadic predictors

  //# Block predictor variables   
    int N_group_vars;                                             //# Number of block structure variables
    int max_N_groups;                                             //# Max number of group labels in any variable
    array[N_group_vars] int N_groups_per_var;                     //# Number of group labels, per variable type
    array[N_id, N_group_vars, N_responses] int block_set;         //# Dataframe holding the group ID codes for each person (rows) for each variable type (cols) for each timepoint

  //# Focal, target, and dyadic predictor variables                                                                                                      
    array[N_id, N_params[1], N_responses] real focal_set;         //# Focal slash decider predictor variables for each timepoint  
    array[N_id, N_params[2], N_responses] real target_set;        //# Target slash alter predictor variables for each timepoint
    array[N_id, N_id, N_params[3], N_responses] real dyad_set;    //# Dyadic predictor variables for each timepoint

    array[N_id, N_params[4], N_responses] real ind_focal_set;     //# Focal slash decider predictor variables for each timepoint  

  //# Network and exposure data
    array[N_id,N_id,N_responses] int outcomes;                    //# Network of ties for each timepoint
    array[N_id,N_id,N_responses] int exposure;                    //# Exposure for each outcome for each timepoint
    array[N_id,N_id,N_responses] int mask;                        //# Censoring mask for each outcome for each timepoint

  //# Diffusing trait data and exposure
    array[N_id,N_responses] int diffusion_outcomes;  //# Network of ties for each timepoint
    array[N_id,N_responses] int diffusion_exposure;  //# Exposure for each outcome for each timepoint
    array[N_id,N_responses] int diffusion_mask;      //# Censoring mask for each outcome for each timepoint

  //# Accessory parameters 
    matrix[22, 2] priors;                            //# Priors in a matrix, see details in the make_priors() function
    int export_network;                              //# Controls export of predictions
    int outcome_mode;                                //# Are outcomes binomial
    int link_mode;                                   //# Link type
}

transformed data{
  //# Refactor to the first predictor slot, becuase it is unity
    array[N_id, N_params[1]-1, N_responses] real focal_predictors;            //# Same as focal_set without first column
    array[N_id, N_params[2]-1, N_responses] real target_predictors;           //# Same as target_set without first column
    array[N_id, N_id, N_params[3]-1, N_responses] real dyad_predictors;       //# Same as dyad_set without first shelf
    array[N_id, N_params[4]-1, N_responses] real ind_focal_predictors;        //# Same as ind_focal_set without first column

  //# Store some key indexes
    array[max_N_groups, N_group_vars, N_responses] int N_per_group;           //# Number of people in each block-type for each group variable
    array[N_group_vars+1] int block_indexes;                                  //# The indexes of each block parameter when stored as a vector instead of ragged array
    int block_param_size;                                                     //# Total number of block-level parameters

  //# Get size of parameters for block model
    block_param_size = 0;                                                     //# Start at zero
    block_indexes[1] = 0;                                                     //# Start at zero for first index 
    
    for(q in 1: N_group_vars){
     block_param_size += N_groups_per_var[q]*N_groups_per_var[q];                       //# Count up number of parameters in each K by K block matrix and add to total
     block_indexes[1+q] = N_groups_per_var[q]*N_groups_per_var[q] + block_indexes[q];   //# Create cummulative sum of block indices, by adding new sum to old sum
     }

  for(t in 1:N_responses){
  //# First fill with scrap
    for(q in 1: N_group_vars){
    for(k in 1: max_N_groups){
     N_per_group[k, q, t] = 0;   
     }}

  //# Now fill in real values
    for(q in 1: N_group_vars){
    for(i in 1:N_id){
     N_per_group[block_set[i,q,t],q,t] += 1;
     }}

  //# Make pruned data for predictor variables, by dropping first column
    if(N_params[1]>1){
     for(i in 2:N_params[1]){
     focal_predictors[ , i-1, t] = focal_set[ , i, t];  
     }}

    if(N_params[2]>1){
     for(i in 2:N_params[2]){
     target_predictors[ , i-1, t] = target_set[ , i , t];  
     }}

    if(N_params[3]>1){
     for(i in 2:N_params[3]){
     dyad_predictors[ , , i-1, t] = dyad_set[ , , i, t];  
     }}

    if(N_params[4]>1){
     for(i in 2:N_params[4]){
     ind_focal_predictors[ , i-1, t] = ind_focal_set[ , i, t];  
     }}

   }
}

parameters{
    //# Block effects, stored as a vector to save space
    vector[block_param_size] block_effects;

    //# Effects of covariate
    vector[N_params[1]-1] focal_effects;
    vector[N_params[2]-1] target_effects;
    vector[N_params[3]-1] dyad_effects; 
    vector[N_params[4]-1] ind_focal_effects;   

    //# Individual learning base-rate
    real lambda;
}

model{
  //# Local storage to make code more readable
    array[N_group_vars] matrix[max_N_groups, max_N_groups] B;  //# Block effects, in array form
    vector[N_group_vars] br;                                   //# Sum of block effects per dyad    
    vector[2*N_responses] scrap;                               //# Local storage         
    array[N_id, N_id, N_responses] real br_network;            //# Block strucuture  
    vector[N_id] weighted_tie_vector;                          //# Tie vector
    vector[N_id] latent_prediction;
    real theta;
    vector[N_id] psi;


    //# The first step, is to transform the vector of block effects into a list of matrices
     for(q in 1:N_group_vars){
      B[q, 1:N_groups_per_var[q], 1:N_groups_per_var[q]] = to_matrix(block_effects[(block_indexes[q]+1):(block_indexes[q+1])], N_groups_per_var[q], N_groups_per_var[q]);
     }

    //# Then put priors on B, which scale loosely with the block size
    for ( q in 1:N_group_vars ){
    for ( i in 1:N_groups_per_var[q] ){
        for ( j in 1:N_groups_per_var[q] ) {
            if ( i==j ) {
                B[q,i,j] ~ normal(logit(priors[10,1]/sqrt(N_per_group[i,q,t])), priors[10,2]);                              //# transfers more likely within groups
            } else {
                B[q,i,j] ~ normal(logit(priors[11,1]/sqrt(N_per_group[i,q,t]*0.5 + N_per_group[j,q,t]*0.5)), priors[11,2]); //# transfers less likely between groups
            }
        }}
    }

    //# Priors on effects of covariates
     focal_effects ~ normal(priors[12,1], priors[12,2]);
     ind_focal_effects ~ normal(priors[12,1], priors[12,2]);
     target_effects ~ normal(priors[13,1], priors[13,2]);
     dyad_effects ~ normal(priors[14,1], priors[14,2]);

     
   for(t in 1:N_responses){
     for(i in 1:N_id){
      for(j in 1:N_id){
       if(i != j){
        for(q in 1:N_group_vars){
          br[q] = B[q,block_set[i,q,t], block_set[j,q,t]]; //# Extract all of the block components for this dyad
         }
         br_network[i,j,t] = sum(br);
         }}}
    }

        lambda ~ normal(0, 2.5);

    //# Likelihood
    for(t in 2:N_responses){
     for(i in 1:N_id){
        if(diffusion_mask[i,1]==0){
          if(diffusion_outcomes[i,t-1]==0){

      //################################################################ Model of individual learning
          theta = inv_logit(lambda + dot_product(ind_focal_effects,  to_vector(ind_focal_predictors[i, ,t])) );
         
      //################################################################ Model of social learning
          psi = inv_logit(to_vector(br_network[i,,t]) + 
                                    dot_product(dyad_effects[t],  to_vector(dyad_predictors[i, j, , t])) +
                                    dot_product(focal_effects[t],  to_vector(focal_predictors[i, ,t])) +
                                    dot_product(target_effects[t],  to_vector(target_predictors[i, ,t]))
                                             );


      //################################################################ Integration of weights
      latent_prediction[i] = theta[i] + (1-theta[i])*inv_logit2(sum(network_filtered[i,]));

      //# Build links
      //# social_information_weights gives the effective strength of social information conditional on covariate data: focal, target, and dyadic
      //# social_information_weights = inv_logit(to_vector(br_network[i,,t]) + 
      //#                                       dot_product(dyad_effects[t],  to_vector(dyad_predictors[i, j, , t])) +
      //#                                       dot_product(focal_effects[t],  to_vector(focal_predictors[i, ,t])) +
      //#                                       dot_product(target_effects[t],  to_vector(target_predictors[i, ,t]))
      //#                                       );

      //# weighted_tie_vector = to_vector(A[i,]) .* social_information_weights;

     

      latent_prediction[i,t] = lambda; //#+ dot_product(weighted_tie_vector, to_vector(diffusion_outcomes[,t-1]));
  
      
      //################################################################ Outcome models // DONE
      //# Bernoulli outcomes, two link options
      if(outcome_mode==1){
        if(link_mode==1){
         diffusion_outcomes[i,t] ~ bernoulli_logit(latent_prediction[i,t]);  
        }
        if(link_mode==2){
         diffusion_outcomes[i,t] ~ bernoulli(Phi(latent_prediction[i,t]));  
        }
       }

      //# Binomial outcomes, two link options
      if(outcome_mode==2){
        if(link_mode==1){
         diffusion_outcomes[i,t] ~ binomial_logit(diffusion_exposure[i,t], latent_prediction[i,t]);  
        }
        if(link_mode==2){
         diffusion_outcomes[i,t] ~ binomial(diffusion_exposure[i,t], Phi(latent_prediction[i,t]));  
        }
       }

      //# Poisson outcomes, one link options
      if(outcome_mode==3){
         diffusion_outcomes[i,t] ~ poisson_log(latent_prediction[i,t]);  
        }

         }}  
   }}
 }
