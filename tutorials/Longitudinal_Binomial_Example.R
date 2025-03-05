############################################## Baboon example
library(STRAND)
library(stringr)
library(ggplot2)
library(psych)

# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors


############################################## 
# Load data
 data(Baboon_Longitudinal_Data)
 d = Baboon_Longitudinal_Data

############################################## 
# Process data into longitudinal form
dat_long = NULL
# Loop over y days of data, making a day-specific network data-set
# Some days are missing outcome data. This is dealt with via the mask layer
# if outcome data is missing, mask[i,j]=1
# currently, missings in the predictors aren't supported in STRAND, but will be eventually
for(y in 1:14){
 # Merge data
 dat_long[[y]] = make_strand_data(
  outcome = list("Affiliative" = d$Affiliative[[y]]),
  exposure = list("Affiliative" = d$Exposure[[y]]),
  mask = list("Affiliative" = d$Mask[[y]]),
  block_covariates = NULL, 
  individual_covariates = d$Individual, 
  dyadic_covariates = list("Presenting" = t(d$Presenting[[y]])),
  longitudinal = TRUE,
  outcome_mode="binomial",
  link_mode="logit"
  )
 }
names(dat_long) = paste("Time", c(1:14))

############################################## 
# Fit model with time-varying slopes
fit_2a = fit_longitudinal_model(
 long_data=dat_long,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ Presenting,
 coefficient_mode="varying",
 random_effects_mode="varying",
 bandage_penalty = -1,
  mode="mcmc",
 stan_mcmc_parameters = list(seed = 1, chains = 1, init=0,
    parallel_chains = 1, refresh = 1, iter_warmup = 1000,
    iter_sampling = 1000, max_treedepth = 12)
  )

res_2a = summarize_longitudinal_bsrm_results(fit_2a)

############################################## 
# Visualize results
pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,2))
longitudinal_plot(fit_2a, type="dyadic", save_plot="Baboon_dyadic_long_1.pdf", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,2,9,3,4))
longitudinal_plot(fit_2a, type="generalized", save_plot="Baboon_generalized_long_1.pdf", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,7,5,10,8,6))
longitudinal_plot(fit_2a,type="coefficient", 
    parameter_set = list(
    focal="Age", target="Age", 
    focal="SexMale", target="SexMale",
    dyadic="Presenting"),
    palette=pal,
    normalized=TRUE,
    height=4, width=9,
    save_plot="Slopes_Baboon.pdf")


############################################## 
# Fit model with time-invariant slopes
fit_2b = fit_longitudinal_model(
 long_data=dat_long,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ Presenting,
 coefficient_mode="fixed",
 random_effects_mode="fixed",
  mode="mcmc",
 stan_mcmc_parameters = list(seed = 1, chains = 1, 
    parallel_chains = 1, refresh = 1, iter_warmup = 100,
    iter_sampling = 1000, max_treedepth = 12),
 priors=NULL
  )

res_2b = summarize_longitudinal_bsrm_results(fit_2b)

############################################## 
# Visualize results
pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,2))
longitudinal_plot(fit_2b, type="dyadic", save_plot="Baboon_dyadic_long_2.pdf", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,2,9,3,4))
longitudinal_plot(fit_2b, type="generalized", save_plot="Baboon_generalized_long_2.pdf", palette = pal, height=6, width=6.5)


############################################## 
# Make a merged plot
bab2b_set = c("focal effects coeffs (out-degree), Time 1 - Age", 
"focal effects coeffs (out-degree), Time 1 - SexMale", 
"target effects coeffs (in-degree), Time 1 - Age", 
"target effects coeffs (in-degree), Time 1 - SexMale",
"dyadic effects coeffs, Time 1 - Presenting")

to_add = res_2b$summary[which(res_2b$summary$Variable %in% bab2b_set),]

    parameter_set = list(
    focal="Age", target="Age", 
    focal="SexMale", target="SexMale",
    dyadic="Presenting")

base_set = longitudinal_plot_c(fit_2a, parameters=as.vector(unlist(parameter_set)), type=names(parameter_set),
    plot=FALSE,
    normalized=FALSE,
    export_as_table=TRUE)

to_add$short_names = c("Focal - Time 1 - Age", "Focal - Time 1 - SexMale",
                "Target - Time 1 - Age", "Target - Time 1 - SexMale",
                 "Dyadic - Time 1 - Presenting")
to_add$time_point = rep("Time 0", 5)
to_add$time_point_int = rep(0, 5)
to_add$extra_short_names = c("Focal - Age", "Focal - SexMale",
                "Target - Age", "Target - SexMale",
                 "Dyadic - Presenting")

to_add$type_set = c("Focal", "Focal", "Target", "Target", "Dyadic")

to_add$Model = "Fixed"
base_set$Model = "Varying"

full_set = rbind(base_set,to_add)


full_set2 = full_set[,c("Variable", "time_point", "extra_short_names", "extra_short_names", "Median", "HPDI:0.05", "HPDI:0.95", "Mean", "SD", "time_point_int", "type_set")]
colnames(full_set2) = c("Variable", "Layer", "Target" , "Base" , "Median", "L", "H", "Mean","SD","LayerNumeric","Type")
      
     Diff = as.numeric(full_set2$H)-as.numeric(full_set2$L)   
     full_set2$Median = as.numeric(full_set2$Median)/Diff
     full_set2$L = as.numeric(full_set2$L)/Diff
     full_set2$H =  as.numeric(full_set2$H)/Diff

     full_set2$Model = ifelse(full_set2$"LayerNumeric"==0, "Fixed", "Varying")
     full_set2$Model2 = ifelse(full_set2$"Model"=="Fixed", "Fixed", full_set2$Type)

     full_set2$Model2 = factor(full_set2$Model2)
     full_set2$Model2 = factor(full_set2$Model2, levels=c("Fixed", "Dyadic", "Focal", "Target"))

p = ggplot(full_set2, aes(x=LayerNumeric, y=as.numeric(Median), ymin=as.numeric(L), ymax=as.numeric(H), group=Target, color=Target))+ 
     geom_linerange(size=1, position = position_dodge(width = 0.3)) + facet_grid(.~Model2, scales="free", space="free") + 
     geom_point(size=2, position = position_dodge(width = 0.3))+
     geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
     labs(y="Effect size", x="Time step") + 
     theme(strip.text.x = element_text(size=12,face="bold"), 
      strip.text.y = element_text(size=12,face="bold"),
      axis.text = element_text(size=12),
      axis.title = element_text(size=14, face="bold"))+
     theme(strip.text.y = element_text(angle = 360)) + 
    # coord_flip() + 
     theme(panel.spacing = grid::unit(1, "lines")) + scale_color_manual(values = pal) + 
     theme(legend.position="bottom") + theme(legend.title = element_blank()) + scale_x_continuous(breaks=1:14,expand = c(0, 0.95))

p

ggsave("Slopes_Baboon_merged.pdf", p, height=6, width=13.5)


res_2a = summarize_longitudinal_bsrm_results(fit_2a)


strand_VPCs(fit_2a, n_partitions = 3, HPDI=0.9, include_reciprocity=TRUE, mode="adj")
multiplex_plot(fit_2a, type="dyadic", mode="cor", HPDI=0.9)
multiplex_plot(fit_2a, type="dyadic", mode="adj", HPDI=0.9)




rlkjcorr = function (n, K, eta = 1) 
{
    stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
    stopifnot(eta > 0)
    f <- function() {
        alpha <- eta + (K - 2)/2
        r12 <- 2 * rbeta(1, alpha, alpha) - 1
        R <- matrix(0, K, K)
        R[1, 1] <- 1
        R[1, 2] <- r12
        R[2, 2] <- sqrt(1 - r12^2)
        if (K > 2) 
            for (m in 2:(K - 1)) {
                alpha <- alpha - 0.5
                y <- rbeta(1, m/2, alpha)
                z <- rnorm(m, 0, 1)
                z <- z/sqrt(crossprod(z)[1])
                R[1:m, m + 1] <- sqrt(y) * z
                R[m + 1, m + 1] <- sqrt(1 - y)
            }
        return(crossprod(R))
    }
    R <- replicate(n, f())
    if (dim(R)[3] == 1) {
        R <- R[, , 1]
    }
    else {
        R <- aperm(R, c(3, 1, 2))
    }
    return(R)
}

silly_prod = function(X,v){
    Y=X
    for(i in 1:length(v)) Y[i,i]=X[i,i]*v[i]
    return(Y)
}

merge = function(A,D){
X1 = cbind(A,D)
X2 = cbind(D,A)
Y = rbind(X1,X2)
return(Y)
} 

K = 5
A = rlkjcorr(n=1, K=K, eta=2.5)
B = rlkjcorr(n=1, K=K, eta=2.5)
C = runif(K,-1, 1)
D = silly_prod(B,C)

M = merge(A,D)
L = chol(M)







to_corr = function(K){
    N = choose(K,2)
    M = diag(rep(0,K))
    M[lower.tri(M)] = rnorm(N, 0, 2)

    Z = tanh(M)

    X = matrix(0, nrow=K, ncol=K)
    X[1,] = c(1, rep(0, K-1))

    for(i in 2:K){
       for(j in 1:K){
         if(i == j) X[i,j] = sqrt(1 - sum(X[i,1:(j-1)]^2));
          if(i > j) X[i,j] = Z[i,j] * sqrt(1 - sum(X[i,1:(j-1)]^2));
       } 
    }

   R = X %*% t(X)

   return(X)
}

return_pars = function(X){
 K = nrow(X)
    Z1 = matrix(0, nrow=K, ncol=K)
    Z1[1,] = rep(0, K)

    for(i in 2:K){
       for(j in K:1){
         if(i == j) Z1[i,j] = 0;
         if(j>1){
          if(i > j) Z1[i,j] = X[i,j] / sqrt(1 - sum(X[i,1:(j-1)]^2));
          }
         if(j==1){
            Z1[i,j] = X[i,j]
         }
       } 
    }
   

  M1 =  0.5*(log(1 + Z1) - log(1 - Z1))
  return(M1)
}

bob = return_pars(X)
bob2 = return_pars(rlkjcorr(n=1,K=8,eta=4))


to_corr2 = function(M){
    Z = tanh(M)

    X = matrix(0, nrow=K, ncol=K)
    X[1,] = c(1, rep(0, K-1))

    for(i in 2:K){
       for(j in 1:K){
         if(i == j) X[i,j] = sqrt(1 - sum(X[i,1:(j-1)]^2));
          if(i > j) X[i,j] = Z[i,j] * sqrt(1 - sum(X[i,1:(j-1)]^2));
       } 
    }

   R = X %*% t(X)

   return(X)
}





A = rlkjcorr(n=1,K=4,eta=4)
B = rlkjcorr(n=1,K=4,eta=4)
C = runif(4,-0.01, 0.01)
D = silly_prod(B,C)

L = chol(merge(A,D))

X = t(L)






L = matrix(0, nrow=4, ncol=4)
H = rlkjcorr(n=1,K=2,eta=1.5)

L[1:2,1:2] = t(chol(H))
L[3,1] = runif(1, -1, 1)
L[3,2] = runif(1, -1, 1)*sqrt(1 - L[3,1]^2)
L[4,1] = sum(L[2,1:2]*L[3,1:2])
L[3,3] = sqrt(1 - sum(L[3,1:2]^2))

  thresh = (L[2,1] - L[3,3] - L[3,1]*L[4,1])/L[3,2]

  Q = L[2,1]/L[3,3] - (L[3,1]/L[3,3])*L[4,1]
  G = L[3,2]/L[3,3]
  X = L[4,2]

  A = (1 + G^2)
  B = 2*Q*G
  C = -(1 - L[4,1]^2 - Q^2)

  S1 = (-B - sqrt(B^2 - 4*A*C))/(2*A)
  S2 = (-B + sqrt(B^2 - 4*A*C))/(2*A)
  
  thresh_0 = c(S1, S2)

  if(L[3,2] > 0) L[4,2] = runif(1, max(thresh, thresh_0[1]), thresh_0[2])  
  if(L[3,2] < 0) L[4,2] = runif(1, thresh_0[1], min(thresh, thresh_0[2])) 

L[4,3] = (L[2,1] - sum(L[3,1:2]*L[4,1:2]))/L[3,3]
L[4,4] = sqrt(1 - sum(L[4,1:3]^2))

L %*% t(L)








Q = L[2,1]/L[3,3] - (L[3,1]/L[3,3])*L[4,1]
G = L[3,2]/L[3,3]
X = L[4,2]

(1 + G^2)*X^2 + 2*Q*G*X  < 1 - L[4,1]^2 - Q^2

A = (1 + G^2)
B = 2*Q*G
C = -(1 - L[4,1]^2 - Q^2)

S1 = (-B + sqrt(B^2 - 4*A*C))/(2*A)
S2 =  (B + sqrt(B^2 - 4*A*C))/(2*A)

X_test = mean(S1,S2)
R = (1 + G^2)*X_test^2 + 2*Q*G*X_test  < 1 - L[4,1]^2 - Q^2
R

X_test = seq(-1,1, by=0.01)
R = (1 + G^2)*X_test^2 + 2*Q*G*X_test  < 1 - L[4,1]^2 - Q^2
R









L = matrix(0, nrow=4, ncol=4)
H = diag(rep(1,2))
lim = 0.5
H[1,2] = H[2,1] = runif(1, -lim, lim)

L[1:2,1:2] = t(chol(H))
L[3,1] = runif(1, -lim, lim)
L[3,2] = runif(1, -lim, lim)*sqrt(1 - L[3,1]^2)
L[4,1] = sum(L[2,1:2]*L[3,1:2])
L[3,3] = sqrt(1 - sum(L[3,1:2]^2))
L[4,2] = runif(1, -lim, lim)  
L[4,3] = (L[2,1] - sum(L[3,1:2]*L[4,1:2]))/L[3,3]
L[4,4] = sqrt(1 - sum(L[4,1:3]^2))

L %*% t(L)



inverse = function(x){return(solve(x))}

K = 4
A = rlkjcorr(n=1,K=K,eta=5)
B = rlkjcorr(n=1,K=K,eta=5)
C = runif(K,0, 0)
D = silly_prod(B,C)


I = diag(rep(1,K))
Z = diag(rep(0,K))

Schur_thing = A - D %*% inverse(A) %*% D

F1 = rbind(cbind(I,Z),cbind(D %*% inverse(A),I))

F2 = rbind(cbind(A,Z), cbind(Z, Schur_thing))

F3 = rbind(cbind(I,inverse(A) %*% D), cbind(Z,I))

chol(F1 %*% F2 %*% F3)



functions{
  //# Function to build a dyadic reciprocity matrix by hand
  matrix multiply_diag(matrix B, vector C){
    matrix[rows(B),cols(B)] D = B;

    for(i in 1:size(C))
     D[i,i] = B[i,i]*C[i];
    return D;
  }

  matrix build_dr_matrix(matrix A, matrix B, vector C){
    matrix[rows(B),cols(B)] D = multiply_diag(B, C);
    return append_row(append_col(A, multiply_diag(B, C)), append_col(multiply_diag(B, C), A));
  }
}
