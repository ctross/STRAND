#' A function to take inverse logit
#'
#' See McElreaths Rethinking package for details
#'
#' @param x value
#' @return inverse_logit(x)
#' @export

inv_logit = function(x){
    p = 1/(1 + exp(-x))
    p = ifelse(x == Inf, 1, p)
    p
}
