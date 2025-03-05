data{                                                                                                                    
    int N_responses;                           //# number of outcome layers
    int setting;                               //# setting = 1 for multiplex dr, setting = 2 for longitudinal dr, setting = 3 for longitudinal gr
    real bandage_penalty;                      //# stitching strength
    real<lower=0> eta;                         //# prior
}

parameters{         
    cholesky_factor_corr[2*N_responses] dr_L;   
}

transformed parameters{
    matrix[2*N_responses, 2*N_responses] D_corr; 
    D_corr = tcrossprod(dr_L);  
}

model{          
    //# Stitch together the dyadic matrix
    if(setting<3){
    for(m in 1:(N_responses-1)){
    for(n in (m+1):N_responses){
     target += normal_lpdf(D_corr[m+N_responses, n+N_responses] | D_corr[m, n],   bandage_penalty);
     target += normal_lpdf(D_corr[m, n+N_responses]   | D_corr[n, m+N_responses], bandage_penalty);
     }}
    }

    //# Add stitches to the dyadic matrix for longitudinal symmetries
    if(setting==2){
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

    //# Add stitches to the generalized matrix for longitudinal symmetries
    if(setting==3){
     for(k in 1:(N_responses-1)){
     for(m in 1:(N_responses-k)){
      target += normal_lpdf(D_corr[m, m+k] | D_corr[1, k+1], bandage_penalty);
      target += normal_lpdf(D_corr[m, N_responses+m+k] | D_corr[1, N_responses+k+1], bandage_penalty);
      target += normal_lpdf(D_corr[N_responses + m, N_responses + m+k] | D_corr[N_responses + 1, N_responses + k+1], bandage_penalty);
      target += normal_lpdf(D_corr[m + k, N_responses + m] | D_corr[k + 1, N_responses + 1], bandage_penalty);
      }}

     for(m in 1:N_responses){
      target += normal_lpdf(D_corr[m, m+N_responses] | D_corr[1, N_responses+1], bandage_penalty);
     }
    }

    dr_L ~ lkj_corr_cholesky(eta); 
 }
