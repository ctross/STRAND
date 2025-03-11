#' A function to compute the min by col
#'
#' Simpler helper function.
#'
#' @param x Input array.
#' @param margin Margin.
#' @return colMins(x, margin)
#' @export

colMins = function(x, margin=2){
              m = apply(x, MARGIN=c(margin), min, na.rm=TRUE)
              return(m)
             }  
             