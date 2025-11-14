#' A function to compute summary stats from a JAX object
#' 
#' This is a small helper function for JAX ESS and Rhat checking.
#'
#' @param X A set of sample to summarize.
#' @return A table of MCMC summaries.
#' @export

jax_summary = function(X){
 np = reticulate::import("numpy")
 draws = posterior::as_draws_array(reticulate::py_to_r(np$array(X)))
 posterior::summarise_draws(draws)
}
