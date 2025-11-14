#' A data-parser function to translate a matrix of Block ID variables into a tensor of indicators
#'
#' This is an internal function.
#'
#' @param x A data frame of indicators.
#' @param priors A data frame of priors.
#' @return block_set_to_dyadic_block_set(x)
#' @export

block_set_to_dyadic_block_set = function(x, priors = NULL){
    N_vars = ncol(x)
    N_id = nrow(x)
    N_groups_per_var = rep(NA, N_vars)

    for(v in 1:N_vars){
     N_groups_per_var[v] = max(x[,v])
    }

    max_N_groups_per_var = max(N_groups_per_var)

    N_per_group = matrix(1, nrow=N_vars, ncol=max_N_groups_per_var)

    for(v in 1:N_vars){
     N_per_group[v,1:N_groups_per_var[v]] = table(x[,v])
    }

    N_pars = sum(N_groups_per_var^2)

    N_groups_per_var2 = N_groups_per_var^2

    ################################### Parameter structure
    B_V = B_I = B_J = B_In = B_Base = B_SS = array(0, c(N_vars, max_N_groups_per_var, max_N_groups_per_var))

    if(is.null(priors)){
     block_priors =  make_priors()
      } else{
     block_priors = priors
      }

    for(v in 1:N_vars){
    for(i in 1:N_groups_per_var[v]){
    for(j in 1:N_groups_per_var[v]){
     B_V[v,i,j] = v
     B_I[v,i,j] = i
     B_J[v,i,j] = j
     B_In[v,i,j] = 1
     B_Base[v,i,j] = ifelse(i==j, block_priors[10,1], block_priors[11,1])
     B_SS[v,i,j] = sqrt(N_per_group[v,i]*0.5 + N_per_group[v,j]*0.5)
      }
     }
    }

    dat = data.frame(V = c(B_V), I = c(B_I), J = c(B_J), In = c(B_In), Base = c(B_Base), SS = c(B_SS))

    dat2 = dat[which(dat$In==1),]

    dat2$Mu = logit(dat2$Base/dat2$SS)

    if(block_priors[10,1] == block_priors[11,1]){
        stop("block_priors[10,1] and block_priors[11,1] may not be identical. Add at least a tiny offset.")
    }
    
    dat2$Sigma = ifelse(dat2$Base == block_priors[10,1], block_priors[10,2], block_priors[11,2])

    dat2 = dat2[order(dat2$V),]

    ########## Build binary tensor

     Y = array(0, c(N_id, N_id, N_pars))

     for(i in 1:N_id){
        for(j in 1:N_id){
             scrap_i = scrap_j = rep(0, N_pars)
          for(v in 1:N_vars){
              scrap_i[which(dat2$V==v & dat2$I==x[i,v])] = 1
              scrap_j[which(dat2$V==v & dat2$J==x[j,v])] = 1
              }
            Y[i,j,] = scrap_i*scrap_j
            }
        }
     
    return(list(Y=Y, Mu=c(dat2$Mu), Sigma=c(dat2$Sigma)))

    }
