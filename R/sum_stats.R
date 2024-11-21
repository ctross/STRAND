#' A function to summarize samples
#'
#' This is a simple helper function to summarize samples. 
#'
#' @param y A string.
#' @param x Samples.
#' @param z HPDI level
#' @return A vector of summaries.
#' @export
#'

sum_stats = function(y, x, z){
      bob = rep(NA, 7)
      dig = 3
      bob[1] = y
      bob[2] = round(median(x),dig)
      bob[3] = round(HPDI(x, z)[1],dig)
      bob[4] = round(HPDI(x, z)[2],dig)
      bob[5] = round(mean(x),dig)
      bob[6] = round(sd(x),dig)
      bob[7] = round(bayesian_p(x),dig)

      return(bob)
      }
