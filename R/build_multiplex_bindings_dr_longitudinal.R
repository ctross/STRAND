#' A function to add accesssory data to a STRAND data object for use in longitudinal analysis
#'
#' This is an internal function.
#'
#' @param data A STRAND data object.
#' @param verbose If TRUE, then print structure of dyadic reciprocity matrix.
#' @return A STRAND data object with extra information about the structure of the dyadic reciprocity matrix.
#' @export
#'

build_multiplex_bindings_dr_longitudinal = function(data, verbose = FALSE){
  # Compute banding structure
  fill_set = c()
   for(k in 1:(data$N_responses-1)){
    fill_set = c(fill_set, seq(k, 1))
   }
  
  # Build A matrix
  bindings_mat_A = diag(NA, data$N_responses)
  bindings_mat_A[upper.tri(bindings_mat_A)] = fill_set
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0
  
  # Build B matrix
  bindings_mat_B = diag(0, data$N_responses)
  bindings_mat_B[upper.tri(bindings_mat_B)] = fill_set + data$N_responses - 1
  bindings_mat_B[lower.tri(bindings_mat_B)] = 0
  bindings_mat_B = bindings_mat_B + t(bindings_mat_B)
  diag(bindings_mat_B) = 2*data$N_responses - 1

  # Assemble full matrix
  bindings_mat = diag(NA, 2*data$N_responses)
  bindings_mat[upper.tri(bindings_mat)] = 0

  bindings_mat[1:data$N_responses, 1:data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, (1:data$N_responses)+data$N_responses] = t(bindings_mat_A)
  bindings_mat[(1:data$N_responses)+data$N_responses, 1:data$N_responses] = bindings_mat_B

  # Compute lookup table
  dr_indices = which(bindings_mat > 0, arr.ind=TRUE)
  dr_indices = dr_indices[order(dr_indices[,1]),]
  dr_id = bindings_mat[dr_indices]

  # Compute other info 
  data$N_dr_bindings_l = choose(2*data$N_responses,2) # Not sure about this one
  data$N_dr_params_l = 2*data$N_responses - 1
  data$N_dr_indices_l = nrow(dr_indices)
  data$dr_indices_l = dr_indices
  data$dr_id_l = bindings_mat[dr_indices]

  if(verbose == TRUE){
   print(bindings_mat)
  }

  return(data)
 }
 