#' A function to compute summary stats from a JAX object
#' 
#' This is a small helper function for JAX ESS and Rhat checking.
#'
#' @param X A set of sample to summarize.
#' @param n_iter Iterations.
#' @param n_chains Chains.
#' @return A table of MCMC summaries.
#' @export

jax_summary = function(X, n_iter, n_chains){
 np = reticulate::import("numpy")
 draws = posterior::as_draws_array(reshape_jax_arrays(reticulate::py_to_r(np$array(X)), n_iter, n_chains))
 posterior::summarise_draws(draws)
}

