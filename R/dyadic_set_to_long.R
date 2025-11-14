#' A function to map dyadic predictors to a long-form dyad list.
#'
#' This is an internal function.
#'
#' @param X A 3d-array.
#' @return dyadic_set_to_long(X)
#' @export

dyadic_set_to_long = function(X){
  N_id = nrow(X)
  V = dim(X)[3]
  Y = array(NA, c((N_id*(N_id-1)/2), 2, V))   

  locs = make_dyadic_edgelist(N_id) 

  for(v in 1:V){
    for(q in 1:nrow(locs)){
        Y[q,1,v] = X[locs[q,1],locs[q,2],v]
        Y[q,2,v] = X[locs[q,2],locs[q,1],v]
    }
  }

  return(Y)
}
