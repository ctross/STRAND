#' A function to mean center data
#'
#' This is a simple helper function to mean center data. We recommmend to mean center covariates as it makes parameter interpretation easier.
#'
#' @param 
#' input A vector of data to center.
#' @return A vector of centered data.
#' @export
#' @examples
#' \dontrun{
#' age_centered = center(input=age)
#' }
#'

center = function(input){
  y = input - mean(input, na.rm=TRUE)
  return(y)
}
