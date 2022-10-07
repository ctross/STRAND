#' A function to simulate Bernoulli data
#'
#' See McElreaths Rethinking package for details
#'
#' @param 
#' n Number of samples
#' @param 
#' prob Probability
#' @return Binary samples
#' @export

rbern = function (n, prob = 0.5) 
{
    rbinom(n, size = 1, prob = prob)
}
