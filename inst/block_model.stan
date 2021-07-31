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
}

transformed data{
 matrix[N_id, N_params[1]-1] focal_individual_predictors; 
 matrix[N_id, N_params[2]-1] target_individual_predictors; 
 real dyad_individual_predictors[N_id, N_id, N_params[3]-1]; 

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

    //# Effects of covariate
    vector[N_params[1]-1] focal_effects;
    vector[N_params[2]-1] target_effects;
    vector[N_params[3]-1] dyad_effects;  
}

model{
  vector[2] sr[N_id];
  matrix[N_id, N_id] dr;

    //# Priors on effects of covariates
    focal_effects ~ normal(0,1);
    target_effects ~ normal(0,1);
    dyad_effects ~ normal(0,1);

    for(i in 1:N_id){
     vector[2] sr_terms;

     sr_terms[1] = dot_product(focal_effects,  to_vector(focal_individual_predictors[i]));
     sr_terms[2] = dot_product(target_effects,  to_vector(target_individual_predictors[i]));  

     sr[i] = sr_terms;
     }

    //# Dyadic priors 
    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
     dr[i,j] =  dot_product(dyad_effects,  to_vector(dyad_individual_predictors[i, j, ]));
     dr[j,i] =  dot_product(dyad_effects,  to_vector(dyad_individual_predictors[j, i, ]));
     }}

    for(i in 1:N_id){
     dr[i,i] = -99; //# ignore this :)
    }

    //# priors for B
    for ( i in 1:N_groups ){
        for ( j in 1:N_groups ) {
            if ( i==j ) {
                B[i,j] ~ normal(logit(1/sqrt(N_id)), 0.5);   //# transfers more likely with groups
            } else {
                B[i,j] ~ normal(logit(0.1/sqrt(N_id)), 0.5); //# transfers less likely between groups
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


