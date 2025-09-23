#' A function to check for errors
#'
#' Taken from berryFunctions by Berry Boessenkool.
#' See https://github.com/brry/berryFunctions for details.
#'
#' @param expr expression
#' @param tell Logical: Should the error message be printed?
#' @param force Logical: Should an error be returned if the expression is not an error? 
#' @return Boolean for error being thrown.
#' @export

is_error = function (expr, tell = FALSE, force = FALSE){
    expr_name = deparse(substitute(expr))
    test = try(expr, silent = TRUE)
    iserror = inherits(test, "try-error")
    if(tell) 
        if(iserror) 
            message("Note in is.error: ", test)
    if(force) 
        if(!iserror) 
            stop(expr_name, " is not returning an error.", call. = FALSE)
    return(iserror)
}
