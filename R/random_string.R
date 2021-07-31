#' A function to make random ID codes
#' 
#' This is a small helper function to create random ID codes.
#'
#' @param 
#' N Number of ID codes to make
#' @param 
#' length Length of ID codes
#' @return A vector of unique ID codes
#' @export
#' @examples
#' \dontrun{
#' IDs = random_string(N=100, length=5)
#' }
#'

random_string = function(N=1, length=12){
    randomString = rep(NA, N)               
    for (i in 1:N){
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS), length, replace=TRUE), collapse="")
                  }

     if(sum(duplicated(randomString))>0){
        stop("At least one ID code string was duplicated. Consider a larger string length.")
                  }

    return(randomString)
    }

