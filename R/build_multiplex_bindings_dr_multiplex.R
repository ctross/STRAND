#' A function to add accesssory data to a STRAND data object for use in multiplex analysis
#'
#' This is an internal function.
#'
#' @param data A STRAND data object.
#' @param verbose If TRUE, then print structure of dyadic reciprocity matrix.
#' @return A STRAND data object with extra information about the structure of the dyadic reciprocity matrix.
#' @export
#'

build_multiplex_bindings_dr_multiplex = function(data, verbose = FALSE){
  # Build A matrix
  bindings_mat_A = diag(NA, data$N_responses)
  bindings_mat_A[upper.tri(bindings_mat_A)] = seq(1,choose(data$N_responses,2))
  bindings_mat_A[lower.tri(bindings_mat_A)] = 0

  # Build B matrix
  bindings_mat_B = diag(0, data$N_responses)
  bindings_mat_B[lower.tri(bindings_mat_B)] = seq(1,choose(data$N_responses,2)) + choose(data$N_responses,2)
  bindings_mat_B[upper.tri(bindings_mat_B)] = 0
  bindings_mat_B = bindings_mat_B + t(bindings_mat_B)
  
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
  data$N_dr_bindings = 2*choose(data$N_responses,2)
  data$N_dr_params = data$N_dr_bindings + data$N_responses
  data$N_dr_indices = nrow(dr_indices)
  data$dr_indices = dr_indices
  data$dr_id = bindings_mat[dr_indices]

  if(verbose == TRUE){
   print(bindings_mat)
  }

  return(data)
 }
 