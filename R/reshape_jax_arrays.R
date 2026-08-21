#' A function to reshape a JAX samples object
#' 
#' This is a small helper function for JAX ESS and Rhat checking.
#'
#' @param samps A set of samples to summarize.
#' @param n_iter Iterations.
#' @param n_chains Chains.
#' @return A reshaped array that posterior can work with.
#' @export

reshape_jax_arrays = function(samps, n_iter, n_chains){
  expected_n = n_iter * n_chains

  # Require an array/vector with the expected number of draws
  if(length(samps) < expected_n || dim(samps)[1] != expected_n){
      stop(
        "First dimension must equal n_iter * n_chains.",
        "Expected: ", expected_n, 
        "Found: ", dim(samps)[1]
      )
    }

  # Preserve all dimensions after the first
    old_dim = dim(samps)

    new_dim = c(
      n_iter,
      n_chains,
      old_dim[-1]
    )

    # Preserve dimension names where possible
    old_dimnames = dimnames(samps)

    new_dimnames = c(
      list(
        iteration = as.character(seq_len(n_iter)),
        chain = as.character(seq_len(n_chains))
      ),
      old_dimnames[-1]
    )

   samps_new = array(
      samps,
      dim = new_dim,
      dimnames = new_dimnames
    )

   dim_sizes = dim(samps_new)
   dim(samps_new) = c(dim_sizes[1], dim_sizes[2], prod(dim_sizes[-(1:2)]))
 return(samps_new)
}
