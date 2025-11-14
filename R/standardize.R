#' A function to standardize data
#'
#' This is a simple helper function to standardize data. We recommmend to standardize covariates as it makes parameter interpretation easier.
#'
#' @param input A vector of data to standardize.
#' @param type Divide by "max" or "sd"?
#' @param center Subtract mean?
#' @return A vector of standardize data.
#' @export
#' @examples
#' \dontrun{
#' age_standardize = standardize_strand(input=age)
#' }
#'

standardize_strand = function(input, center=TRUE, type="sd"){
  if(center==TRUE){
   y = input - mean(input, na.rm=TRUE)
  }

  if(center==FALSE){
   y = input
  }

  if(type=="sd"){
   z = y/sd(y, na.rm=TRUE)
  }

  if(type=="max"){
   z = y/max(y, na.rm=TRUE)
  }

  if(!type %in% c("sd","max")){
   stop("type must be sd or max.")
  }

  return(z)
}
