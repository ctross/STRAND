#' A function to compute the max by col
#'
#' Simpler helper function.
#'
#' @param x Input array.
#' @param margin Margin.
#' @return colMaxs(x, margin)
#' @export

colMaxs = function(x, margin=2){
              m = apply(x, MARGIN=c(margin), max, na.rm=TRUE)
              return(m)
             }  
