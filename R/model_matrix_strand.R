#' A wrapper function for model.matrix that reports more helpful errors.
#'
#'
#' @param x Model equation.
#' @param d Data object.
#' @param name Name of equation for error message.
#' @return Boolean for error being thrown.
#' @export

model_matrix_strand = function(x, d, name){
    if(is_error(model.matrix(x, d))){
        stop(paste0("It appears that you have included a predictor variable in your regression equation that does not appear as (i.e., exactly match) a named variable in the make_strand_data object. \n Check your equations and data declarations.",
                    " This error occured in the: ", name, " equation."))
    }else{
        return(model.matrix(x, d))
    }
}
