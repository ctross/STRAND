#' A function to add accesssory data to a STRAND data object for use in multiplex analysis
#'
#' This is an internal function.
#'
#' @param data A STRAND data object.
#' @param verbose If TRUE, then print structure of the bindings.
#' @return A STRAND data object with extra information about the structure of the dyadic reciprocity matrix.
#' @export
#'

build_multiplex_bindings_dr_numpyro = function(data, verbose = FALSE){
  # Compute banding structure
  K = data$N_responses
  fill_set = c(1:choose(K,2))

  # Build A matrix
  bindings_mat_A = diag(NA, K)
  bindings_mat_A[upper.tri(bindings_mat_A)] = fill_set
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0

  bind_set = which(bindings_mat_A>0, arr.ind = TRUE)
  bind_set = bind_set[order(bind_set[,1]),]

  bind_set_1 = bind_set
  bind_set_2 = bind_set
  bind_set_3 = bind_set
  bind_set_4 = bind_set

  bind_set_2[,1] = bind_set[,1] + K
  bind_set_2[,2] = bind_set[,2] + K

  bind_set_3[,1] = bind_set[,2] 
  bind_set_3[,2] = bind_set[,1] + K

  bind_set_4[,1] = bind_set[,1] 
  bind_set_4[,2] = bind_set[,2] + K


  out_all = rbind(bind_set_1, bind_set_3)
  in_all = rbind(bind_set_2, bind_set_4)


  if(verbose == TRUE){
   print(out_all)
   print(in_all)
  }
  
  data$numpyro_dr_bindings_out = out_all
  data$numpyro_dr_bindings_in = in_all
   return(data)
}

