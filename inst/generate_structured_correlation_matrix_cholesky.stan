functions {
  vector lb_ub_lp(vector y, real lb, real ub) {
    target += log(ub - lb) + log_inv_logit(y) + log1m_inv_logit(y);
    
    return lb + (ub - lb) * inv_logit(y);
  }
  
  real lb_ub_lp(real y, real lb, real ub) {
    target += log(ub - lb) + log_inv_logit(y) + log1m_inv_logit(y);
    
    return lb + (ub - lb) * inv_logit(y);
  }
  
  matrix cholesky_corr_constrain_lp(int K, vector raw, int N_blocks,
                                    array[,] int res_index,
                                    array[] int res_id, vector lb, vector ub) {
    matrix[K, K] L = rep_matrix(0, K, K);
    int cnt = 1;
    int N_res = num_elements(res_id);
    vector[N_blocks] x_cache;
    array[N_blocks] int res_id_cnt = ones_int_array(N_blocks);
    int res_row = 1;
    
    L[1, 1] = 1;
    L[2, 1] = lb_ub_lp(raw[cnt], lb[cnt], ub[cnt]);
    L[2, 2] = sqrt(1 - L[2, 1] ^ 2);
    cnt += 1;
    
    if (res_index[res_row, 1] == 2) {
      x_cache[res_id[res_row]] = L[2, 1];
      res_id_cnt[res_id[res_row]] += 1;
      res_row += 1;
    }
    
    for (i in 3 : K) {
      if (res_index[res_row, 1] == i && res_index[res_row, 2] == 1) {
        if (res_id_cnt[res_id[res_row]] == 1) {
          L[i, 1] = lb_ub_lp(raw[cnt], lb[cnt], ub[cnt]);
          x_cache[res_id[res_row]] = L[i, 1];
          res_id_cnt[res_id[res_row]] += 1;
          cnt += 1;
        } else {
          L[i, 1] = x_cache[res_id[res_row]];
        }
        res_row += 1;
      } else {
        L[i, 1] = lb_ub_lp(raw[cnt], lb[cnt], ub[cnt]);
        cnt += 1;
      }
      
      L[i, 2] = sqrt(1 - L[i, 1] ^ 2);
      real l_ij_old = log1m(L[i, 1] ^ 2);
      for (j in 2 : i - 1) {
        real b1 = dot_product(L[j, 1 : (j - 1)], L[i, 1 : (j - 1)]);
        real stick_length = exp(0.5 * l_ij_old);
        real low = max({-stick_length, (lb[cnt] - b1) / L[j, j]});
        real up = min({stick_length, (ub[cnt] - b1) / L[j, j]});
        
        if (res_index[res_row, 1] == i && res_index[res_row, 2] == j) {
          if (res_id_cnt[res_id[res_row]] == 1) {
            L[i, j] = lb_ub_lp(raw[cnt], low, up);
            x_cache[res_id[res_row]] = L[i, j] * L[j, j] + b1;
            res_id_cnt[res_id[res_row]] += 1;
            cnt += 1;
          } else {
            L[i, j] = (x_cache[res_id[res_row]] - b1) / L[j, j];
            target += -log(L[j, j]);
          }
          res_row = res_row == N_res ? N_res : res_row + 1;
        } else {
          L[i, j] = lb_ub_lp(raw[cnt], low, up);
          cnt += 1;
        }
        l_ij_old = log_diff_exp(l_ij_old, 2 * log(abs(L[i, j])));
      }
      L[i, i] = exp(0.5 * l_ij_old);
    }
    return L;
  }
}

data{   
  //# Array dimension variables                                                                                                                            
    int N_responses;                           //# Number of outcome networks

  //# Dyadic reciprocity control parameters
    real<lower=0> eta;
    int<lower=0> N_dr_params;
    int<lower=0> N_dr_indices;
    int<lower=0> N_dr_bindings;
    array[N_dr_indices, 2] int dr_indices;
    array[N_dr_indices] int dr_id;
}

transformed data{
  //# Dyadic reciprocity control parameters
   int<lower=0> N_off_diag = ((2*N_responses) * ((2*N_responses) - 1)) %/% 2;
   vector[N_off_diag] lb;
   vector[N_off_diag] ub;

  //# Dyadic reciprocity control parameters
    lb = rep_vector(-1, N_off_diag);
    ub = rep_vector(1, N_off_diag);
}

parameters{         
    vector[N_dr_params] dr_par_set;                 
}

transformed parameters{
    matrix[2*N_responses, 2*N_responses] dr_L = cholesky_corr_constrain_lp(2*N_responses, dr_par_set, N_dr_bindings, dr_indices, dr_id, lb, ub);                                                      
}

model{
  //# Dyadic priors for social relations model
    dr_L ~ lkj_corr_cholesky(eta);
 }

generated quantities {
  matrix[2*N_responses, 2*N_responses] D_corr;
  D_corr = multiply_lower_tri_self_transpose(dr_L);
}
