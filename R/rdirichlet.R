#' A function to simulate Dirichlet data
#'
#' See https://github.com/dkahle/dirichlet/blob/master/R/dirichlet.R for details
#'
#' @param n number of samples
#' @param mu mean vector
#' @return samples
#' @export

rdirichlet = function(n, mu){
  normalize = function(.) . / sum(.)
  samps = vapply(mu, function(al) rgamma(n, al, 1), numeric(n))
  if (n == 1) {
    matrix(normalize(samps), nrow = 1, ncol = length(samps))
  } else {
    t(apply(samps, 1, normalize))
  }  
}