#' A function to get Highest Posterior Density Interval
#'
#' This is a simple helper function to get density intervals. See McElreaths Rethinking for details.
#'
#' @param samples A vector of samples to summarize.
#' @param prob Interval range.
#' @return A vector of intervals.
#' @export
#'

HPDI = function(samples, prob = 0.89){
  concat = function (...){
    paste(..., collapse = "", sep = "")
    }

    coerce.list = c("numeric", "matrix", "data.frame", "integer", 
        "array")
    if (inherits(samples, coerce.list)) {
        samples = coda::as.mcmc(samples)
    }
    x = sapply(prob, function(p) coda::HPDinterval(samples, 
        prob = p))
    n = length(prob)
    result = rep(0, n * 2)
    for (i in 1:n) {
        low_idx = n + 1 - i
        up_idx = n + i
        result[low_idx] = x[1, i]
        result[up_idx] = x[2, i]
        names(result)[low_idx] = concat("|", prob[i])
        names(result)[up_idx] = concat(prob[i], "|")
    }
    return(result)
}
