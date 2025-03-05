functions{
    //# Shifted and scaled inv_logit link
     real inv_logit_shifted(real x) {
     return((inv_logit(x) - 0.5)*2);
    }

    //# CES function
    vector CES(vector K, vector L, real alpha, real sigma, real eta) {
      real rho = (sigma - 1)/sigma;
      return(pow(alpha*pow(K, rho) + (1-alpha)*pow(L, rho), (eta/rho)));
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
    array[N_id, N_params[1]-1, N_responses] real focal_predictors;                      //# Same as focal_set without first column
    array[N_id, N_params[2]-1, N_responses] real target_predictors;                     //# Same as target_set without first column
    array[N_id, N_id, N_params[3]-1, N_responses] real dyad_predictors;                 //# Same as dyad_set without first shelf
    array[N_id, N_params[4]-1, N_responses] real ind_focal_predictors;                  //# Same as ind_focal_set without first column
    array[N_id,N_id,N_responses] real A;                                                //# Network of ties for each timepoint

  //# Store some key indexes
    array[max_N_groups, N_group_vars] real N_per_group;                                 //# Number of people in each block-type for each group variable
    array[N_group_vars+1] int block_indexes;                                            //# The indexes of each block parameter when stored as a vector instead of ragged array
    int block_param_size;                                                               //# Total number of block-level parameters

  //# Get size of parameters for block model
    block_param_size = 0;                                                               //# Start at zero
    block_indexes[1] = 0;                                                               //# Start at zero for first index 
    
    for(q in 1: N_group_vars){
     block_param_size += N_groups_per_var[q]*N_groups_per_var[q];                       //# Count up number of parameters in each K by K block matrix and add to total
     block_indexes[1+q] = N_groups_per_var[q]*N_groups_per_var[q] + block_indexes[q];   //# Create cummulative sum of block indices, by adding new sum to old sum
     }

  //# First fill with scrap
    for(q in 1: N_group_vars){
    for(k in 1: max_N_groups){
     N_per_group[k, q] = 0;   
     }}

  for(t in 1:N_responses){
  //# Now fill in real values
    for(q in 1: N_group_vars){
    for(i in 1:N_id){
     N_per_group[block_set[i,q,t],q] += 1.0 / N_responses;
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

 // Temp SRI style
     for(t in 1:N_responses){
     for(i in 1:N_id){
      for(j in 1:N_id){
        if(i == j){
          A[i,j,t] = 0;
        } else{
          A[i,j,t] = outcomes[i,j,t]*1.0 / exposure[i,j,t]; 
         }
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

    //# Directed information. Defines how social information from directed networks is integrated.
    real<lower=0, upper=1> alpha;
    real<lower=0> sigma;
    real<lower=0> eta;
}

model{
  //# Local storage to make code more readable
    array[N_group_vars] matrix[max_N_groups, max_N_groups] B;  //# Block effects, in array form
    vector[N_group_vars] br;                                   //# Sum of block effects per dyad        
    vector[N_id] attention_weighted_network;                   
    real latent_prediction;
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
            if ( i == j ) {
                B[q,i,j] ~ normal(logit(priors[10,1]/sqrt(N_per_group[i,q])), priors[10,2]);                              //# transfers more likely within groups
            } else {
                B[q,i,j] ~ normal(logit(priors[11,1]/sqrt(N_per_group[i,q]*0.5 + N_per_group[j,q]*0.5)), priors[11,2]);   //# transfers less likely between groups
            }
        }}
    }

    //# Priors on effects of covariates
     lambda ~ normal(0, 2.5);

     alpha ~ beta(2, 2);
     sigma ~ normal(2, 2);
     eta ~ normal(1, 0.33);

     focal_effects ~ normal(priors[12,1], priors[12,2]);
     ind_focal_effects ~ normal(priors[12,1], priors[12,2]);
     target_effects ~ normal(priors[13,1], priors[13,2]);
     dyad_effects ~ normal(priors[14,1], priors[14,2]);
  
    //# Likelihood
    for(t in 2:N_responses){
     for(i in 1:N_id){
      if(diffusion_outcomes[i, t-1]==0){
      //################################################################ Model of individual learning
          //# theta gives the weight of individual learning as a function of focal characteristics
          theta = inv_logit(lambda + dot_product(ind_focal_effects,  to_vector(ind_focal_predictors[i, ,t])) );
         
      //################################################################ Model of social learning
        for(j in 1:N_id){
          //# Loop over block variables and add to sum
          for(q in 1:N_group_vars){
           br[q] = B[q, block_set[i,q,t], block_set[j,q,t]]; //# Extract all of the block components for dyad i j
           }

          //# psi gives the weight of social learning as a function of block, focal, target, and dyadic characteristics
          psi[j] = inv_logit(sum(br) + 
                             dot_product(focal_effects,  to_vector(focal_predictors[i, ,t])) +
                             dot_product(target_effects,  to_vector(target_predictors[j, ,t])) +
                             dot_product(dyad_effects,  to_vector(dyad_predictors[i, j, , t])) 
                                       );
          }
         
         //# to get to social learning from weights, we need to multply psi by the effective exposure to social information  
         attention_weighted_network =  psi .* CES(to_vector(A[i,,t]) .* to_vector(diffusion_outcomes[, t-1]), 
                                                  to_vector(A[,i,t]) .* to_vector(diffusion_outcomes[, t-1]), 
                                                  alpha, sigma, eta);
           //# The CES function integrates information on outgoing ties and incoming ties, and the status of alters as having the binary trait of interest.
           //# The first two terms are ties weights times indicators for the binary trait.
           //# The next three terms are parameters. Alpha gives the relative importance of i to j ties, relative to j to i ties, in i learning the diffusing outcome from j.
           //# sigma gives the elasticity of substitution, and eta the returns to scale.

      //######################## Integration of individual and social learning weights
        latent_prediction = theta + (1-theta)*inv_logit_shifted(sum(attention_weighted_network));
        //# Note, we map the sum social learning weight in (0, infinity) to (0,1) with inv_logit_shifted()
      
      //################################################################ Outcome model
      //# Bernoulli outcomes
        diffusion_outcomes[i,t] ~ bernoulli(latent_prediction);  
     }}  
   }
 }
