#' A function to generate correlation matrices with special inner symmetry
#'
#' This is an internal function.
#'
#' @param data A data list.
#' @param eta Prior on the LKJ Cholesky factor.
#' @param mode Which method to use? "cholesky" or "l2norm".
#' @param setting Which matrix to create? "multiplex_dyadic_reciprocity", "longitudinal_dyadic_reciprocity", or "longitudinal_generalized_reciprocity".
#' @param mcmc_parameters A list of Stan parameters that often need to be tuned.
#' @return A Stan object.
#' @export
#'

generate_structured_correlation_matrix = function(data, eta, mode, setting, mcmc_parameters) {
  data$eta = eta

  if (setting == "multiplex_dyadic_reciprocity") {
    data$setting = 1
    data = build_multiplex_bindings_dr_multiplex(data)
  }

  if (setting == "longitudinal_dyadic_reciprocity") {
    data$setting = 2
    data = build_multiplex_bindings_dr_longitudinal(data)
  }

  if (setting == "longitudinal_generalized_reciprocity") {
    data$setting = 3
    data = build_multiplex_bindings_sr_longitudinal(data)

    # Have to do a little renaming, so that the Stan file works for both dyadic and generalized matrices
    data$N_dr_bindings = data$N_sr_bindings
    data$N_dr_params = data$N_sr_params
    data$N_dr_indices = data$N_sr_indices
    data$dr_indices = data$sr_indices
    data$dr_id = data$sr_id
  }


  if (mode == "cholesky") {
    # Correlation matrix methods contributed to STRAND by Sean Pinkney
    # Copyright 2025 Sean Pinkney <sean.pinkney@gmail.com>
    # Subject to the BSD 3-Clause License
    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","generate_structured_correlation_matrix_cholesky.stan"))
  }

  if (mode == "l2norm") {
    model = cmdstanr::cmdstan_model(paste0(path.package("STRAND"),"/","generate_structured_correlation_matrix_l2norm.stan"))
  }

  fit = model$sample(
    data = unclass(data),
    seed = mcmc_parameters$seed,
    chains = mcmc_parameters$chains,
    parallel_chains = mcmc_parameters$parallel_chains,
    refresh = mcmc_parameters$refresh,
    iter_warmup = mcmc_parameters$iter_warmup,
    iter_sampling = mcmc_parameters$iter_sampling,
    max_treedepth = mcmc_parameters$max_treedepth,
    adapt_delta = mcmc_parameters$adapt_delta,
    init = mcmc_parameters$init
  )

  bob = list(data = data, fit = fit)
  return(bob)
}
