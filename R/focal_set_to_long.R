#' A function to map focal or sender predictors to a long-form dyad list.
#'
#' This is an internal function.
#'
#' @param X A data.frame.
#' @return focal_set_to_long(X)
#' @export

focal_set_to_long = function(X){
  N_id = nrow(X)
  V = ncol(X)
  Y = array(NA, c( (N_id*(N_id-1)/2), 2, V))   

  locs = make_dyadic_edgelist(N_id) 

  for(v in 1:V){
    for(q in 1:nrow(locs)){
        Y[q,1,v] = X[locs[q,1],v]
        Y[q,2,v] = X[locs[q,2],v]
    }
  }

  return(Y)
}
