#' A function to take inverse logit shifted
#'
#' Map x from (0, inf) to (0, 1)
#'
#' @param x value
#' @return inv_logit_shifted(x)
#' @export

inv_logit_shifted = function(x){
 return((inv_logit(x) - 0.5)*2)
}
