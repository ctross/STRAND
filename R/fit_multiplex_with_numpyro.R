#' A NumPyro implementation of the STRAND multiplex latent-network SRM
#'
#' @param data A STRAND data object.
#' @param bandage_penalty How tight are the dyadic bindings? 
#' @param mcmc_parameters A parameter list to control NUTS sampling. NumPyro uses the same control set as Stan.
#' @return A NumPyro results object.
#' @export
#' @examples
#' \dontrun{
#' res = fit_multiplex_with_numpyro(data, mcmc_parameters)
#' }
#'

fit_multiplex_with_numpyro = function(
    data,
    bandage_penalty,
    mcmc_parameters
){

########################### Input checks   
 if(data$link_mode == 2){stop("NumPyro back-end only supports logit links in the SRM. Fit with Stan using 'mcmc' if you want a probit model.")}  
 if(sum(data$mask) != 0){stop("NumPyro back-end does not support outcome masking. Fit with Stan using 'mcmc'.")}               
 
# Import numpy, jax, numpyro, and numpyro.distributions
  reticulate::py_require(c("numpy>=1.27","numpyro>=0.18", "jax>=0.7", "jaxlib>=0.7"))
  np = reticulate::import("numpy")
  jax = reticulate::import("jax")
  numpyro = reticulate::import("numpyro")
  dist = reticulate::import("numpyro.distributions")
  jnp = reticulate::import("jax.numpy")
  numpyro_init = reticulate::import("numpyro.infer.initialization")
  
  numpyro_multiplex = NULL
  reticulate::source_python(paste0(path.package("STRAND"),"/","numpyro_multiplex.py"))       

########################### NumPyro backend prep
  dyadic_block_set = block_set_to_dyadic_block_set(data$block_set, data$priors)
  block_mu = dyadic_block_set[[2]]
  block_sigma = dyadic_block_set[[3]]
  long_focal_set = focal_set_to_long(data$focal_set)
  long_target_set = target_set_to_long(data$target_set)
  long_dyad_set = dyadic_set_to_long(data$dyad_set)
  long_block_set = dyadic_set_to_long(dyadic_block_set[[1]])
  long_outcome_set = dyadic_set_to_long(data$outcomes)
  long_exposure_set = dyadic_set_to_long(data$exposure)
  long_ids = make_dyadic_edgelist(nrow(data$focal_set))

  # Convert to numpy arrays                         
  long_focal_set_np = np$array(long_focal_set)
  long_target_set_np = np$array(long_target_set)
  long_dyad_set_np = np$array(long_dyad_set)
  long_block_set_np = np$array(long_block_set)
  long_outcome_set_np = np$array(aperm(long_outcome_set, perm = c(3,1,2)))
  long_exposure_set_np = np$array(aperm(long_exposure_set, perm = c(3,1,2)))
  block_mu_np = np$array(block_mu)
  block_sigma_np = np$array(block_sigma)
  long_ids_np = np$array(
   cbind(as.integer(long_ids[,1]), as.integer(long_ids[,2]))-1,
   dtype = np$int32
   )

  long_ids_int = np$array(long_ids_np, np$int32)
  long_ids_int = jnp$astype(long_ids_int, jnp$int32)

  N_dyads = as.integer(dim(long_focal_set)[1])
  N_id = as.integer(data$N_id)
  N_layers = as.integer(dim(long_outcome_set_np)[1])

  N_var_focal = as.integer(dim(long_focal_set)[3])
  N_var_target = as.integer(dim(long_target_set)[3])
  N_var_dyad = as.integer(dim(long_dyad_set)[3])
  N_var_block = as.integer(dim(long_block_set_np)[3])

  A_F = rep(1, N_var_focal)
  A_T = rep(1, N_var_target)
  A_D = rep(1, N_var_dyad)

  A_D[1] = A_F[1] = A_T[1] = 0 
  A_F = np$array(A_F)
  A_T = np$array(A_T)
  A_D = np$array(A_D)

  priors_np = np$array(data$priors, dtype = np$float32)

  prior_12_1 = as.numeric(priors_np[12L,1L])
  prior_12_2 = as.numeric(priors_np[12L,2L])

  prior_13_1 = as.numeric(priors_np[13L,1L])
  prior_13_2 = as.numeric(priors_np[13L,2L])

  prior_14_1 = as.numeric(priors_np[14L,1L])
  prior_14_2 = as.numeric(priors_np[14L,2L])

  prior_15_1 = as.numeric(priors_np[15L,1L])
  prior_15_2 = as.numeric(priors_np[15L,2L])

  prior_16_1 = as.numeric(priors_np[16L,1L])
  prior_16_2 = as.numeric(priors_np[16L,2L])

  prior_17_1 = as.numeric(priors_np[17L,1L])

  prior_18_1 = as.numeric(priors_np[18L,1L])

  prior_23_1 = as.numeric(priors_np[23L,1L])
  prior_23_2 = as.numeric(priors_np[23L,2L])

  outcome_mode = as.integer(data$outcome_mode)
  export_network = as.integer(data$export_network)
  bandage_penalty = as.numeric(bandage_penalty)


  numpyro_dr_bindings_out_1 = np$array(as.integer(data$numpyro_dr_bindings_out[,1])-1L, np$int32)
  numpyro_dr_bindings_out_2 = np$array(as.integer(data$numpyro_dr_bindings_out[,2])-1L, np$int32)

  numpyro_dr_bindings_in_1 = np$array(as.integer(data$numpyro_dr_bindings_in[,1])-1L, np$int32)
  numpyro_dr_bindings_in_2 = np$array(as.integer(data$numpyro_dr_bindings_in[,2])-1L, np$int32)


  mcmc_key = jax$random$PRNGKey(as.integer(mcmc_parameters$seed))
  numpyro$set_host_device_count(n=as.integer(mcmc_parameters$cores))

  MCMC = reticulate::import("numpyro.infer")$MCMC
  NUTS = reticulate::import("numpyro.infer")$NUTS

  kernel = NUTS(numpyro_multiplex,
                target_accept_prob = mcmc_parameters$adapt_delta,
                max_tree_depth = as.integer(mcmc_parameters$max_treedepth),
                init_strategy=numpyro_init$init_to_uniform(radius=mcmc_parameters$init))

  mcmc = MCMC(kernel,
          num_warmup = as.integer(mcmc_parameters$iter_warmup),
          num_samples = as.integer(mcmc_parameters$iter_sampling),
          num_chains = as.integer(mcmc_parameters$chains),
          chain_method = mcmc_parameters$chain_method)

  mcmc$run(mcmc_key, np, jax, numpyro, dist, jnp,
           outcome_mode, export_network, 
           jnp$array(bandage_penalty),
           jnp$array(long_outcome_set_np), 
           jnp$array(long_exposure_set_np), 
           jnp$array(long_ids_int), 
           jnp$array(block_mu_np), 
           jnp$array(block_sigma_np), 
           jnp$array(prior_12_1), jnp$array(prior_12_2), 
           jnp$array(prior_13_1), jnp$array(prior_13_2), 
           jnp$array(prior_14_1), jnp$array(prior_14_2),
           jnp$array(prior_15_1), jnp$array(prior_15_2), 
           jnp$array(prior_16_1), jnp$array(prior_16_2), 
           jnp$array(prior_17_1), jnp$array(prior_18_1), 
           jnp$array(prior_23_1), jnp$array(prior_23_2),
           jnp$array(long_focal_set_np), 
           jnp$array(long_target_set_np), 
           jnp$array(long_dyad_set_np), 
           jnp$array(long_block_set_np),  
           jnp$array(jnp$astype(numpyro_dr_bindings_out_1, jnp$int32)), 
           jnp$array(jnp$astype(numpyro_dr_bindings_out_2, jnp$int32)), 
           jnp$array(jnp$astype(numpyro_dr_bindings_in_1, jnp$int32)), 
           jnp$array(jnp$astype(numpyro_dr_bindings_in_2, jnp$int32)),
           N_dyads, N_id, N_var_focal, N_var_target, N_var_dyad, N_var_block, N_layers,
           jnp$array(A_F), jnp$array(A_T), jnp$array(A_D))

  return(mcmc)
}



