#' A function to add accesssory data to a STRAND data object for use in longitudinal analysis
#'
#' This is an internal function.
#'
#' @param data A STRAND data object.
#' @param verbose If TRUE, then print structure of generalized reciprocity matrix.
#' @return A STRAND data object with extra information about the structure of the generalized reciprocity matrix.
#' @export
#'

build_multiplex_bindings_sr_longitudinal = function(data, verbose = FALSE){
  # Compute banding structure
  fill_set = c()
   for(k in 1:(data$N_responses-1)){
    fill_set = c(fill_set, seq(k, 1))
   }
  
  # Build A matrix
  bindings_mat_A = diag(NA, data$N_responses)
  bindings_mat_A[upper.tri(bindings_mat_A)] = fill_set
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0

  # Build C matrix
  bindings_mat_C = diag(NA, data$N_responses)
  bindings_mat_C[upper.tri(bindings_mat_C)] = fill_set + (data$N_responses - 1)
  bindings_mat_C[lower.tri(bindings_mat_C)] = 0
  
  # Build B matrice
  bindings_mat_B = diag(0, data$N_responses)
  bindings_mat_B[upper.tri(bindings_mat_B)] = fill_set + 2*(data$N_responses - 1)
  bindings_mat_B[lower.tri(bindings_mat_B)] = 0

  bindings_mat_Bt = diag(0, data$N_responses)
  bindings_mat_Bt[upper.tri(bindings_mat_Bt)] = fill_set + 3*(data$N_responses - 1)
  bindings_mat_Bt[lower.tri(bindings_mat_Bt)] = 0

  bindings_mat_B = bindings_mat_B + t(bindings_mat_Bt)
  diag(bindings_mat_B) = 4*(data$N_responses - 1) + 1

  # Assemble full matrix
  bindings_mat = diag(NA, 2*data$N_responses)
  bindings_mat[upper.tri(bindings_mat)] = 0

  bindings_mat[1:data$N_responses, 1:data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, (1:data$N_responses)+data$N_responses] = t(bindings_mat_C)
  bindings_mat[(1:data$N_responses)+data$N_responses, 1:data$N_responses] = bindings_mat_B

  # Compute lookup table
  sr_indices = which(bindings_mat > 0, arr.ind=TRUE)
  sr_indices = sr_indices[order(sr_indices[,1]),]
  sr_id = bindings_mat[sr_indices]

  # Compute other info 
  data$N_sr_bindings_l = choose(2*data$N_responses,2) # Not sure about this one
  data$N_sr_params_l = 4*(data$N_responses - 1) + 1
  data$N_sr_indices_l = nrow(sr_indices)
  data$sr_indices_l = sr_indices
  data$sr_id_l = bindings_mat[sr_indices]

  if(verbose == TRUE){
   print(bindings_mat)
  }

  return(data)
 }
