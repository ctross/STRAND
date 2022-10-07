#' A function to take logit
#'
#' See McElreaths Rethinking package for details
#'
#' @param 
#' x value
#' @return logit(x)
#' @export

logit = function (x) 
{
    log(x) - log(1 - x)
}
