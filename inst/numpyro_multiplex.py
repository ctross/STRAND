def numpyro_multiplex(
    np, jax, numpyro, dist, jnp,
    outcome_mode, 
    export_network, 
    bandage_penalty,
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
    prior_17_1, prior_18_1, 
    prior_23_1, prior_23_2,
    long_focal_set_np, 
    long_target_set_np, 
    long_dyad_set_np, 
    long_block_set_np, 
    numpyro_dr_bindings_out_1, 
    numpyro_dr_bindings_out_2, 
    numpyro_dr_bindings_in_1, 
    numpyro_dr_bindings_in_2, 
    N_dyads, N_id, 
    N_var_focal, N_var_target, N_var_dyad, 
    N_var_block, N_layers,
    A_F, A_T, A_D):

    # Observation noise (for Gaussian outcome)
    error_sigma = numpyro.sample("error_sigma", dist.TruncatedNormal(loc=jnp.full((N_layers, 1, 1), prior_23_1), scale=jnp.full((N_layers, 1, 1), prior_23_2), low=jnp.full((N_layers, 1, 1), 0.0)))

    # Block, focal, target, dyad effects
    focal_effects = numpyro.sample("focal_effects", dist.Normal(loc=jnp.full((N_layers, N_var_focal), prior_12_1), scale=jnp.full((N_layers, N_var_focal), prior_12_2)))
    target_effects = numpyro.sample("target_effects", dist.Normal(loc=jnp.full((N_layers, N_var_target), prior_13_1), scale=jnp.full((N_layers, N_var_target), prior_13_2)))
    dyad_effects = numpyro.sample("dyad_effects",dist.Normal(loc=jnp.full((N_layers, N_var_dyad), prior_14_1), scale=jnp.full((N_layers, N_var_dyad), prior_14_2)))
    block_effects = numpyro.sample("block_effects", dist.Normal(loc=jnp.broadcast_to(block_mu_np, (N_layers, N_var_block)), scale=jnp.broadcast_to(block_sigma_np, (N_layers, N_var_block))))

    # Linear predictor for each observation
    mu = (
        jnp.einsum("dxf,lf->ldx", long_focal_set_np, focal_effects * A_F) +
        jnp.einsum("dxf,lf->ldx", long_target_set_np, target_effects * A_T) +
        jnp.einsum("dxf,lf->ldx", long_dyad_set_np, dyad_effects * A_D) +
        jnp.einsum("dxf,lf->ldx", long_block_set_np, block_effects)
    )

    # Dyad random effects
    if export_network == 1:
        dr_raw = numpyro.sample("dr_raw", dist.Normal(jnp.zeros((int(N_layers)*2, N_dyads)), 1))
    else:
        dr_raw = numpyro.sample("dr_raw", dist.Normal(jnp.zeros((int(N_layers)*2, N_dyads)), 1), infer={"num_samples": 0})

    dr_L = numpyro.sample("dr_L", dist.LKJCholesky(int(2)*int(N_layers), prior_18_1))
    dr_sigma = numpyro.sample("dr_sigma", dist.TruncatedNormal(loc=jnp.full((int(N_layers),), prior_16_1), scale=jnp.full((int(N_layers),), prior_16_2), low=0))
    dr_sigma_temp = jnp.expand_dims(jnp.repeat(dr_sigma, 2), 1)
    dr_sigma_scaled = dr_sigma_temp * dr_L
    dr = jnp.matmul(dr_sigma_scaled, dr_raw)
    dr_long = jnp.stack(jnp.split(dr, 2), axis=2)

    # Sender and receiver effects
    sr_raw = numpyro.sample("sr_raw", dist.Normal(jnp.zeros((int(N_layers)*2, N_id)), 1))
    sr_L = numpyro.sample("sr_L", dist.LKJCholesky(int(N_layers)*2, prior_17_1))
    sr_sigma = numpyro.sample("sr_sigma", dist.TruncatedNormal(loc=jnp.full((int(N_layers)*2,), prior_15_1), scale=jnp.full((int(N_layers)*2,), prior_15_2),low=0))
    sr_sigma_temp = jnp.expand_dims(sr_sigma, 1)
    sr_sigma_scaled = sr_sigma_temp * sr_L
    gr = jnp.matmul(sr_sigma_scaled, sr_raw)

    # Expand sender/receiver effects into observation-level data
    gr_list = jnp.split(gr, 2)
    gr_sender = gr_list[0]
    gr_receiver = gr_list[1]

    S_i = gr_sender[:, long_ids_int[:, 0]]
    S_j = gr_sender[:, long_ids_int[:, 1]]
    R_i = gr_receiver[:, long_ids_int[:, 0]]
    R_j = gr_receiver[:, long_ids_int[:, 1]]

    gr_long = jnp.stack([S_i + R_j, S_j + R_i], axis=2)

    # Outcome models
    linear_model = mu + dr_long + gr_long

    if outcome_mode == 1:
        # Bernoulli logits
        numpyro.sample("obs", dist.Bernoulli(logits=linear_model), obs=long_outcome_set_np)

    elif outcome_mode == 2:
        # Binomial with provided total_count and logits
        numpyro.sample("obs", dist.Binomial(total_count=long_exposure_set_np, logits=linear_model),
                       obs=long_outcome_set_np)

    elif outcome_mode == 3:
        # Poisson using exp of linear predictor as rate
        numpyro.sample("obs", dist.Poisson(rate=jnp.exp(linear_model)), obs=long_outcome_set_np)

    elif outcome_mode == 4:
        # Gaussian observation
        numpyro.sample("obs", dist.Normal(loc=linear_model, scale=error_sigma), obs=long_outcome_set_np)

    G_corr = numpyro.deterministic("G_corr", jnp.matmul(sr_L, sr_L.T))
    D_corr = numpyro.deterministic("D_corr", jnp.matmul(dr_L, dr_L.T))

    numpyro.factor("dyadic_bandage_constraint", -0.5 * jnp.sum(((D_corr[numpyro_dr_bindings_in_1, numpyro_dr_bindings_in_2] - D_corr[numpyro_dr_bindings_out_1, numpyro_dr_bindings_out_2])/bandage_penalty) ** 2))
