functions{
    //# probability of observed variables relevant to ij dyad
    real prob_sgij( int[] sij, int[] sji, int tie, vector fpri, vector fprj, vector rtti, vector rttj) {
       vector[2] y;
      
    //# model likelihood
    if(tie==0){
            //# prob i says i helps j
        y[1] = bernoulli_logit_lpmf(sij[1] | fpri[1]);

            //# prob j says i helps j
        y[2] = bernoulli_logit_lpmf(sji[2] | fprj[2]);

      }

    else{ //#  if(tie==1){
        //# prob i says i helps j
        y[1] = bernoulli_logit_lpmf(sij[1] | rtti[1]);

        //# prob j says i helps j
        y[2] = bernoulli_logit_lpmf(sji[2] | rttj[2]);
      }
        return(sum(y));
    }
}

data{
    int N_networktypes;           //# Number of network types
    int N_id;                     //# Number of respondents
    int N_groups;                 //# Number of observed groups in block strucutre
    int N_responses;              //# Number of responses in the self-report network, 2 is double sampled
    int group[N_id];              //# Group ID code of individual i
    int s[N_id,N_id,N_responses]; //# Networks of responses
    vector[N_id] predictor_1;     //# Individual level covariate
}

transformed data{
 real S;
 real penalty;

 S = 0;
 for(i in 1:N_id){
   for(j in 1:N_id){
     if(i != j)
   S += (s[i,j,1] + s[j,i,2] + 0.0)/2.0;
 }}

  penalty = S^(1/1.5);
}

parameters{
    //######################################################## Observation model
    //# Measurement model
    vector<lower=0, upper=1>[N_networktypes] false_positive_rate;
    vector<lower=0, upper=1>[N_networktypes] recall_of_true_ties;

    vector<lower=0>[N_networktypes] fpr_sigma;
    vector<lower=0>[N_networktypes] rtt_sigma;

    vector[N_networktypes] fpr_raw[N_id];
    vector[N_networktypes] rtt_raw[N_id];

    vector[N_networktypes] fpr_effects_1;
    vector[N_networktypes] rtt_effects_1;   


    //########################################################### Latent Netowrk
    //# SBM + SRM model
    real out_block;
    real in_block;

    vector<lower=0>[2] sr_sigma;  //# Variation of sender-receiver effects
    cholesky_factor_corr[2] sr_L;
    vector[2] sr_raw[N_id];

    real<lower=0> dr_sigma;     //# Variation of dyadic effects
    cholesky_factor_corr[2] dr_L;
    matrix[N_id, N_id] dr_raw;

    //# Effects of covariate
    vector[N_responses] sr_effects_1;
}

model{
  vector[2] sr[N_id];
  matrix[N_id, N_id] dr;
  vector[N_networktypes] fpr[N_id];
  vector[N_networktypes] rtt[N_id];

  vector[2] scrap;
  matrix[N_groups,N_groups] B;

  matrix[N_id, N_id] p;
  matrix[N_id, N_id] mixed_p;

    //# Priors on effects of covariates
    sr_effects_1 ~ normal(0,1);
    fpr_effects_1 ~ normal(0,1);
    rtt_effects_1 ~ normal(0,1);

    //# Priors for measurement model
    false_positive_rate ~ beta(1,10);
    recall_of_true_ties ~ beta(10,1);

    fpr_sigma ~ exponential(1);
    rtt_sigma ~ exponential(1);

    for(i in 1:N_id){
    fpr_raw[i] ~ normal(0,1);
    rtt_raw[i] ~ normal(0,1);
    }

    for(i in 1:N_id){
    fpr[i] = logit(false_positive_rate) + predictor_1[i]*fpr_effects_1 + fpr_sigma .* fpr_raw[i];
    rtt[i] = logit(recall_of_true_ties) + predictor_1[i]*rtt_effects_1 + rtt_sigma .* rtt_raw[i];
    }

    //# Sender-receiver priors for social relations model
    for(i in 1:N_id)
    sr_raw[i] ~ normal(0,1);
    sr_sigma ~ exponential(1);
    sr_L ~ lkj_corr_cholesky(2.5);

    for(i in 1:N_id)
    sr[i] = predictor_1[i]*sr_effects_1 + diag_pre_multiply(sr_sigma, sr_L) * sr_raw[i];

    //# Dyadic priors for social relations model
    to_vector(dr_raw) ~ normal(0,1);
    dr_sigma ~ exponential(1);
    dr_L ~ lkj_corr_cholesky(2.5);

    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
     scrap[1] = dr_raw[i,j];
     scrap[2] = dr_raw[j,i];
     scrap = rep_vector(dr_sigma, 2) .* (dr_L*scrap);
     dr[i,j] = scrap[1];
     dr[j,i] = scrap[2];
     }}

    for(i in 1:N_id){
     dr[i,i] = -99; //# ignore this :)
    }

    //# priors for B
    in_block ~ normal(logit(0.1/sqrt(N_id)), 0.5);
    out_block ~ normal(logit(0.01/sqrt(N_id)), 0.5);

    for ( i in 1:N_groups ){
        for ( j in 1:N_groups ) {
            if ( i==j ) {
                B[i,j] = in_block; //# transfers more likely within groups
            } else {
                B[i,j] = out_block; //# transfers less likely between groups
            }
        }}

    //# likelihood
    for ( i in 1:N_id ) {
        for ( j in 1:N_id ) {
            if ( i != j ) {
                vector[2] terms;

                //# consider each possible state of true tie and compute prob of data
                for(tie in 0:1) {
                  terms[tie+1] = prob_sgij(s[i,j,], s[j,i,], tie, fpr[i], fpr[j], rtt[i], rtt[j]);
                }
      p[i,j] = inv_logit( B[group[i], group[j]] + sr[i,1] + sr[j,2] + dr[i,j] );  //# Model as a mixture distribution
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
    matrix[N_id, N_id] p_tie_out;
    matrix[N_id, N_id] p;

    {                
                vector[2] terms;
                int tie;
                vector[N_networktypes] alpha;
                vector[N_networktypes] beta;
                vector[N_networktypes] fpr[N_id];
                vector[N_networktypes] rtt[N_id];
                matrix[N_groups,N_groups] B;
                vector[2] sr[N_id];
                matrix[N_id, N_id] dr;
                vector[2] scrap;
            

    for(i in 1:N_id){
     sr[i] = predictor_1[i]*sr_effects_1 + diag_pre_multiply(sr_sigma, sr_L) * sr_raw[i];
    }

    for(i in 1:(N_id-1)){
    for(j in (i+1):N_id){
     scrap[1] = dr_raw[i,j];
     scrap[2] = dr_raw[j,i];
     scrap = rep_vector(dr_sigma, 2) .* (dr_L*scrap);
     dr[i,j] = scrap[1];
     dr[j,i] = scrap[2];
    }}

    for(i in 1:N_id){
     dr[i,i] = -99; //# ignore this :)
    }

    for(i in 1:N_id){
    fpr[i] = logit(false_positive_rate) + predictor_1[i]*fpr_effects_1 + fpr_sigma .* fpr_raw[i];
    rtt[i] = logit(recall_of_true_ties) + predictor_1[i]*rtt_effects_1 + rtt_sigma .* rtt_raw[i];
    }

    for ( i in 1:N_groups ){
    for ( j in 1:N_groups ) {
       if ( i==j ) {
         B[i,j] = in_block; //# transfers more likely within groups
            } else {
         B[i,j] = out_block; //# transfers less likely between groups
      }
    }}


    for ( i in 1:N_id ) {
        for ( j in 1:N_id ) {
            if ( i != j ) {
                // consider each possible state of true tie and compute prob of data
                p[i,j] = inv_logit( B[group[i], group[j]] + sr[i,1] + sr[j,2] + dr[i,j]);

                tie = 0;
                terms[1] = 
                    log1m( p[i,j] ) + 
                    prob_sgij(s[i,j,], s[j,i,], tie, fpr[i], fpr[j], rtt[i], rtt[j]);

                tie = 1;
                terms[2] = 
                    log( p[i,j] ) + 
                    prob_sgij(s[i,j,], s[j,i,], tie, fpr[i], fpr[j], rtt[i], rtt[j]);

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






  