#' A function to map a long-form dyad list to a matrix style.
#'
#' This is an internal function.
#'
#' @param X A set of S samples from an array of S x Ndyad x 2.
#' @param N_id Rows in the dyad list.
#' @return long_to_dyadic_set(X, N_id)
#' @export

long_to_dyadic_set = function(X, N_id){
  S = dim(X)[1]
  Y = array(0, c(S, N_id, N_id))   

  locs = make_dyadic_edgelist(N_id) 
  L = nrow(locs)

 for(s in 1:S){
  Y[cbind(rep(s,L), locs[,1], locs[,2])] = X[s, , 1]
  Y[cbind(rep(s,L), locs[,2], locs[,1])] = X[s, , 2]
 }

  return(Y)
}


