#' A function to calculate prob of sign error, or Bayesian p value
#'
#' @param x value
#' @return bayesian_p(x)
#' @export

bayesian_p = function(x){
 M_x = median(x)

 if(M_x<0){
   N_x = length(x)
   P_x = length(which(x>0))
  } else{
   N_x = length(x)
   P_x = length(which(x<0))
  }
 return(P_x/N_x)
    }
