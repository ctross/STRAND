#' A function to add accesssory data to a STRAND data object for use in longitudinal analysis
#'
#' This is an internal function.
#'
#' @param data A STRAND data object.
#' @param verbose If TRUE, then print structure of the bindings.
#' @return A STRAND data object with extra information about the structure of the dyadic reciprocity matrix.
#' @export
#'

build_multiplex_bindings_gr_longitudinal_numpyro = function(data, verbose = FALSE){
  # Compute banding structure
   K = data$N_responses
  fill_set = c()
   for(k in 1:(K-1)){
    fill_set = c(fill_set, seq(k, 1))
   }
  
  # Build A matrix
  bindings_mat_A = diag(NA, K)
  bindings_mat_A[upper.tri(bindings_mat_A)] = fill_set
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0


  # Build gr bindings 1
  out_list = in_list = NULL

  for(k in 1:(K-1)){
   out_set = which(bindings_mat_A==k,arr.ind=TRUE)
   out_set_rep = out_set[-1,, drop=FALSE]
   in_set_rep = matrix(rep(out_set[1,], dim(out_set_rep)[1]), nrow=dim(out_set_rep)[1], ncol=2, byrow=TRUE)
 
   out_list[[k]] = out_set_rep
   in_list[[k]] = in_set_rep
  }

  out_full_bindings1 = do.call(rbind,out_list)
  in_full_bindings1 = do.call(rbind,in_list)


  # Build gr bindings 2
  out_full_bindings2 = out_full_bindings1
  out_full_bindings2[,2] = out_full_bindings2[,2] + K

  in_full_bindings2 = in_full_bindings1
  in_full_bindings2[,2] = in_full_bindings2[,2] + K
 
  out_dr_bindings = cbind(rep(1, K-1), rep(1+K, K-1))
  in_dr_bindings = cbind(c(2:K),c(2:K) + K )


  # Build gr bindings 3
  out_full_bindings3 = out_full_bindings2
  out_full_bindings3[,1] = out_full_bindings3[,1] + K

  in_full_bindings3 = in_full_bindings2
  in_full_bindings3[,1] = in_full_bindings3[,1] + K

 
  # Build gr bindings 4
  out_full_bindings4 = out_full_bindings1
  out_full_bindings4[,1] = out_full_bindings4[,1] + K

  in_full_bindings4 = in_full_bindings1
  in_full_bindings4[,1] = in_full_bindings4[,1] + K
 

  out_all = rbind(out_full_bindings1, out_full_bindings2, out_dr_bindings, out_full_bindings3, out_full_bindings4)
  in_all = rbind(in_full_bindings1, in_full_bindings2, in_dr_bindings, in_full_bindings3, in_full_bindings4)

  if(verbose == TRUE){
   print(out_all)
   print(in_all)
  }
  
  data$numpyro_gr_long_bindings_out = out_all
  data$numpyro_gr_long_bindings_in = in_all
   return(data)
}
