#' A function to update mcmc control parameters
#'
#' @param user_params User-supplied updates
#' @return a list of parameters
#' @export

merge_mcmc_parameters = function(user_params = list()){
  default_mcmc_parameters = list(
    seed = 1, 
    chains = 1, 
    parallel_chains = 1, 
    refresh = 1, 
    iter_warmup = 500, 
    iter_sampling = 500, 
    max_treedepth = 12, 
    adapt_delta = 0.95, 
    chain_method = "vectorized", 
    cores=1, 
    init = 2
    )
  modifyList(default_mcmc_parameters, user_params)
}
