#' A NumPyro implementation of the STRAND longitudinal latent-network SRM
#'
#' @param data A STRAND data object.
#' @param bandage_penalty How tight are the dyadic and longitudinal bindings? 
#' @param mcmc_parameters A parameter list to control NUTS sampling. NumPyro uses the same control set as Stan.
#' @return A NumPyro results object.
#' @export
#' @examples
#' \dontrun{
#' res = fit_longitudinal_with_numpyro(data, bandage_penalty, mcmc_parameters)
#' }
#'

fit_longitudinal_with_numpyro = function(
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
  
  numpyro_longitudinal = NULL
  reticulate::source_python(paste0(path.package("STRAND"),"/","numpyro_longitudinal.py"))   

########################### NumPyro backend prep 
  N_layers = data$N_responses
  long_outcome_set = dyadic_set_to_long(data$outcomes)
  long_exposure_set = dyadic_set_to_long(data$exposure)
  long_ids = make_dyadic_edgelist(dim(data$outcomes)[1])
  
  # First, pull time-step 1 to get dims
  dyad_set_scrap = data$dyad_set[,,,1,drop=FALSE]
  block_set_scrap = data$block_set[,,1,drop=FALSE]
  focal_set_scrap = data$focal_set[,,1,drop=FALSE]
  target_set_scrap = data$target_set[,,1,drop=FALSE]

  dyadic_block_set_1 = block_set_to_dyadic_block_set(array(block_set_scrap, dim = dim(block_set_scrap)[1:2]), data$priors)
  block_mu_1 = dyadic_block_set_1[[2]]
  block_sigma_1 = dyadic_block_set_1[[3]]
  long_block_set_1 = dyadic_set_to_long(dyadic_block_set_1[[1]])
  long_focal_set_1 = focal_set_to_long(array(focal_set_scrap, dim = dim(focal_set_scrap)[1:2]))
  long_target_set_1 = target_set_to_long(array(target_set_scrap, dim = dim(target_set_scrap)[1:2]))
  long_dyad_set_1 = dyadic_set_to_long(array(dyad_set_scrap, dim = dim(dyad_set_scrap)[1:3]))

  # Now fill across time-steps
  block_mu = matrix(NA, nrow=N_layers, ncol=length(block_mu_1))
  block_sigma = matrix(NA, nrow=N_layers, ncol=length(block_sigma_1))

  long_block_set = array(NA,c(N_layers, dim(long_block_set_1)))
  long_focal_set = array(NA,c(N_layers, dim(long_focal_set_1)))
  long_target_set = array(NA,c(N_layers, dim(long_target_set_1)))
  long_dyad_set = array(NA,c(N_layers, dim(long_dyad_set_1)))

  for(q in 1:N_layers){
   dyad_set_scrap = data$dyad_set[,,,q,drop=FALSE]
   block_set_scrap = data$block_set[,,q,drop=FALSE]
   focal_set_scrap = data$focal_set[,,q,drop=FALSE]
   target_set_scrap = data$target_set[,,q,drop=FALSE]

   temp = block_set_to_dyadic_block_set(array(block_set_scrap, dim = dim(block_set_scrap)[1:2]), data$priors)
   block_mu[q,] = temp[[2]]
   block_sigma[q,] = temp[[3]]
   long_block_set[q,,,] = dyadic_set_to_long(temp[[1]])
   long_focal_set[q,,,] = focal_set_to_long(array(focal_set_scrap, dim = dim(focal_set_scrap)[1:2]))
   long_target_set[q,,,] = target_set_to_long(array(target_set_scrap, dim = dim(target_set_scrap)[1:2]))
   long_dyad_set[q,,,] = dyadic_set_to_long(array(dyad_set_scrap, dim = dim(dyad_set_scrap)[1:3]))
  }


  # Convert to numpy arrays                         
  long_focal_set_np = np$array(long_focal_set)
  long_target_set_np = np$array(long_target_set)
  long_dyad_set_np = np$array(long_dyad_set)
  long_block_set_np = np$array(long_block_set)
  long_outcome_set_np = np$array(aperm(long_outcome_set, perm = c(3,1,2)))
  long_exposure_set_np = np$array(aperm(long_exposure_set, perm = c(3,1,2)))
  block_mu_np = np$array(block_mu)
  block_sigma_np = np$array(block_sigma)
  mean_block_mu_np = np$array(colMeans(block_mu))
  mean_block_sigma_np = np$array(colMeans(block_sigma))
  long_ids_np = np$array(
   cbind(as.integer(long_ids[,1]), as.integer(long_ids[,2]))-1,
   dtype = np$int32
   )

  long_ids_int = np$array(long_ids_np, np$int32)
  long_ids_int = jnp$astype(long_ids_int, jnp$int32)

  N_dyads = as.integer(dim(long_focal_set)[2])
  N_id = as.integer(data$N_id)
  N_layers = as.integer(dim(long_outcome_set_np)[1])

  N_var_focal = as.integer(dim(long_focal_set)[4])
  N_var_target = as.integer(dim(long_target_set)[4])
  N_var_dyad = as.integer(dim(long_dyad_set)[4])
  N_var_block = as.integer(dim(long_block_set_np)[4])

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
  bandage_penalty = as.numeric(data$bandage_penalty)

  random_effects_mode = as.integer(data$random_effects_mode)
  coefficient_mode = as.integer(data$coefficient_mode)


  numpyro_dr_bindings_out_1 = np$array(as.integer(data$numpyro_dr_bindings_out[,1])-1L, np$int32)
  numpyro_dr_bindings_out_2 = np$array(as.integer(data$numpyro_dr_bindings_out[,2])-1L, np$int32)

  numpyro_dr_bindings_in_1 = np$array(as.integer(data$numpyro_dr_bindings_in[,1])-1L, np$int32)
  numpyro_dr_bindings_in_2 = np$array(as.integer(data$numpyro_dr_bindings_in[,2])-1L, np$int32)


  numpyro_dr_long_bindings_out_1 = np$array(as.integer(data$numpyro_dr_long_bindings_out[,1])-1L, np$int32)
  numpyro_dr_long_bindings_out_2 = np$array(as.integer(data$numpyro_dr_long_bindings_out[,2])-1L, np$int32)

  numpyro_dr_long_bindings_in_1 = np$array(as.integer(data$numpyro_dr_long_bindings_in[,1])-1L, np$int32)
  numpyro_dr_long_bindings_in_2 = np$array(as.integer(data$numpyro_dr_long_bindings_in[,2])-1L, np$int32)


  numpyro_gr_long_bindings_out_1 = np$array(as.integer(data$numpyro_gr_long_bindings_out[,1])-1L, np$int32)
  numpyro_gr_long_bindings_out_2 = np$array(as.integer(data$numpyro_gr_long_bindings_out[,2])-1L, np$int32)

  numpyro_gr_long_bindings_in_1 = np$array(as.integer(data$numpyro_gr_long_bindings_in[,1])-1L, np$int32)
  numpyro_gr_long_bindings_in_2 = np$array(as.integer(data$numpyro_gr_long_bindings_in[,2])-1L, np$int32)


  mcmc_key = jax$random$PRNGKey(as.integer(mcmc_parameters$seed))
  numpyro$set_host_device_count(n=as.integer(mcmc_parameters$cores))

  MCMC = reticulate::import("numpyro.infer")$MCMC
  NUTS = reticulate::import("numpyro.infer")$NUTS

  kernel = NUTS(numpyro_longitudinal,
                target_accept_prob = mcmc_parameters$adapt_delta,
                max_tree_depth = as.integer(mcmc_parameters$max_treedepth),
                init_strategy=numpyro_init$init_to_uniform(radius=mcmc_parameters$init))

  mcmc = MCMC(kernel,
          num_warmup = as.integer(mcmc_parameters$iter_warmup),
          num_samples = as.integer(mcmc_parameters$iter_sampling),
          num_chains = as.integer(mcmc_parameters$chains),
          chain_method = mcmc_parameters$chain_method)

  mcmc$run(mcmc_key, 
         np, jax, numpyro, dist, jnp,
         outcome_mode, export_network, 
         jnp$array(bandage_penalty),
         jnp$array(long_outcome_set_np), 
         jnp$array(long_exposure_set_np), 
         jnp$array(long_ids_int), 
         jnp$array(block_mu_np), 
         jnp$array(block_sigma_np), 
         jnp$array(mean_block_mu_np), 
         jnp$array(mean_block_sigma_np), 
         random_effects_mode, 
         coefficient_mode, 
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
         N_dyads, N_id, N_var_focal, N_var_target, N_var_dyad, N_var_block, N_layers,
         jnp$array(jnp$astype(numpyro_dr_bindings_out_1, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_dr_bindings_out_2, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_dr_bindings_in_1, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_dr_bindings_in_2, jnp$int32)),
         jnp$array(jnp$astype(numpyro_dr_long_bindings_out_1, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_dr_long_bindings_out_2, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_dr_long_bindings_in_1, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_dr_long_bindings_in_2, jnp$int32)),
         jnp$array(jnp$astype(numpyro_gr_long_bindings_out_1, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_gr_long_bindings_out_2, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_gr_long_bindings_in_1, jnp$int32)), 
         jnp$array(jnp$astype(numpyro_gr_long_bindings_in_2, jnp$int32)),
         jnp$array(A_F), jnp$array(A_T), jnp$array(A_D))

  return(mcmc)
}



