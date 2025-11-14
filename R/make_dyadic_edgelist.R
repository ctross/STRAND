#' A function to create a long-form dyad list.
#'
#' This is an internal function.
#'
#' @param N_id A count of nodes.
#' @return make_dyadic_edgelist(N_id)
#' @export

make_dyadic_edgelist = function(N_id){
  Q = expand.grid(1:N_id, 1:N_id)
  Q = Q[,2:1]
  Q = Q[which(Q[,1]<Q[,2]),]
  return(Q)
}
