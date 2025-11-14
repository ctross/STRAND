def numpyro_bsrm_model(
    np, jax, numpyro, dist, jnp, 
    outcome_mode,
    export_network,
    long_outcome_set_np,
    long_exposure_set_np,
    long_ids_int,
    block_mu_np,
    block_sigma_np,
    prior_12_1, prior_12_2,
    prior_13_1, prior_13_2,
    prior_14_1, prior_14_2,
    prior_15_1, prior_15_2,
    prior_16_1, prior_16_2,
    prior_17_1,
    prior_18_1,
    prior_23_1, prior_23_2,
    long_focal_set_np,
    long_target_set_np,
    long_dyad_set_np,
    long_block_set_np,
    N_dyads,
    N_id,
    N_var_focal,
    N_var_target,
    N_var_dyad,
    N_var_block,
    A_F, A_T, A_D):

    # Observation noise (for Gaussian outcome)
    error_sigma = numpyro.sample("error_sigma", dist.TruncatedNormal(loc=prior_23_1, scale=prior_23_2, low=0.0))

    # Block, focal, target, dyad effects
    block_effects = numpyro.sample("block_effects", dist.Normal(block_mu_np, block_sigma_np))
    focal_effects = numpyro.sample("focal_effects", dist.Normal(jnp.full((N_var_focal,), prior_12_1), jnp.full((N_var_focal,), prior_12_2)))
    target_effects = numpyro.sample("target_effects", dist.Normal(jnp.full((N_var_target,), prior_13_1), jnp.full((N_var_target,), prior_13_2)))
    dyad_effects = numpyro.sample("dyad_effects", dist.Normal(jnp.full((N_var_dyad,), prior_14_1), jnp.full((N_var_dyad,), prior_14_2)))

    # Linear predictor for each observation
    mu = (
        jnp.tensordot(long_block_set_np, block_effects, axes=1) +
        jnp.tensordot(long_focal_set_np, focal_effects * A_F, axes=1) +
        jnp.tensordot(long_target_set_np, target_effects * A_T, axes=1) +
        jnp.tensordot(long_dyad_set_np, dyad_effects * A_D, axes=1)
    )

    # Dyad random effects
    dr_L = numpyro.sample("dr_L", dist.LKJCholesky(2, prior_18_1))
    dr_sigma = numpyro.sample("dr_sigma", dist.TruncatedNormal(loc=prior_16_1, scale=prior_16_2, low=0.0))
    dr_sigma_temp = jnp.expand_dims(jnp.repeat(dr_sigma, 2), 1)
    dr_sigma_scaled = dr_sigma_temp * dr_L

    if export_network == 1:
        dr_raw = numpyro.sample("dr_raw", dist.Normal(jnp.zeros((2, N_dyads)), 1))
    else:
        dr_raw = numpyro.sample("dr_raw", dist.Normal(jnp.zeros((2, N_dyads)), 1), infer={"num_samples": 0})

    dr = jnp.transpose(jnp.matmul(dr_sigma_scaled, dr_raw))

    # Sender–receiver effects
    sr_L = numpyro.sample("sr_L", dist.LKJCholesky(2, prior_17_1))
    sr_sigma = numpyro.sample("sr_sigma", dist.TruncatedNormal(loc=jnp.full((2,), prior_15_1), scale=jnp.full((2,), prior_15_2), low=jnp.full((2,), 0.0)))
    sr_sigma_temp = jnp.expand_dims(sr_sigma, 1)
    sr_sigma_scaled = sr_sigma_temp * sr_L
    sr_raw = numpyro.sample("sr_raw", dist.Normal(jnp.zeros((2, N_id)), 1))

    gr = jnp.transpose(jnp.matmul(sr_sigma_scaled, sr_raw))

    # Expand sender/receiver effects into observation-level data
    i_idx = long_ids_int[:, 0]
    j_idx = long_ids_int[:, 1]

    S_i = gr[i_idx, 0]
    S_j = gr[j_idx, 0]
    R_i = gr[i_idx, 1]
    R_j = gr[j_idx, 1]
    gr_long = jnp.stack([S_i + R_j, S_j + R_i], axis=1)

    # Outcome models
    linear_model = mu + dr + gr_long
    if outcome_mode == 1:
        # Bernoulli (binary)
        numpyro.sample("obs", dist.Bernoulli(logits = linear_model), obs = long_outcome_set_np)
        if export_network == 1:
            numpyro.deterministic("p", jax.nn.sigmoid(linear_model))

    elif outcome_mode == 2:
        # Binomial
        numpyro.sample("obs", dist.Binomial(total_count = long_exposure_set_np, logits = linear_model), obs = long_outcome_set_np)
        if export_network == 1:
            numpyro.deterministic("p", jax.nn.sigmoid(linear_model))

    elif outcome_mode == 3:
        # Poisson
        numpyro.sample("obs", dist.Poisson(rate = jnp.exp(linear_model)), obs = long_outcome_set_np)
        if export_network == 1:
            numpyro.deterministic("p", jnp.exp(linear_model))

    elif outcome_mode == 4:
        # Gaussian
        numpyro.sample("obs", dist.Normal(loc = linear_model, scale = error_sigma), obs = long_outcome_set_np)
        if export_network == 1:
            numpyro.deterministic("p", linear_model)

    # Store correlation matrices for diagnostics
    numpyro.deterministic("G_corr", jnp.matmul(sr_L, sr_L.T))
    numpyro.deterministic("D_corr", jnp.matmul(dr_L, dr_L.T))
