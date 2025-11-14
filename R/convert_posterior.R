#' A function to map NumPyro samples to STRAND format
#'
#' This is an internal function, written by Sebastian Sosa.
#'
#' @param posteriors A NumPyro posterior.
#' @return convert_posterior(posteriors)
#' @export

convert_posterior = function(posteriors){
  np = reticulate::import('numpy')
  R_list = reticulate::py_to_r(posteriors)
  for(a in 1:length(R_list)){
    R_list[[a]] = reticulate::py_to_r(np$array(R_list[[a]]))
  }
  return(R_list)
}
