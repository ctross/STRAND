#' A function to check the data type
#'
#' @param x value
#' @return is_numeric_not_binary(x)
#' @export

is_numeric_not_binary = function(x){
  if(is.numeric(x)==FALSE){
    return(FALSE)
    } else{
         if(all(x %in% c(0,1,NA))==TRUE){
           return(FALSE)} 
             else{
              return(TRUE)  
            }
          }
        }