#' Organize Stan output and provide summaries of model parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters
#'
#' @param 
#' input A STRAND model object, obtained by fitting a combined stochastic block and social relations model.
#' @param 
#' include_samples An indicator for the user to specify where raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param 
#' HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A STRAND results object including summary table, a summary list, and samples.
#' @export
#' @examples
#' \dontrun{
#' res = summarize_bsrm_results(input=fit)
#' }
#'

summarize_bsrm_results_with_measurement_bias = function(input, include_samples=TRUE, HPDI=0.9){
    input=fit_cens; include_samples=TRUE; HPDI=0.9
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_bsrm_results_with_measurement_bias() requires a fitted object of class: STRAND Model Object. Please use fit_block_plus_social_relations_model_with_measurement_bias() to run your model.")
    }

    if(attributes(input)$fit_type != "mcmc"){
        stop("Fitted results can only be reorganized for STRAND model objects fit using MCMC. Variational inference or optimization can be used in Stan
              during experimental model runs, but final inferences should be based on MCMC sampling.")   
    }

    ###################################################### Create samples 
    stanfit = posterior::as_draws_rvars(input$fit$draws())

    ################### Network model parameters
    sr_sigma = posterior::draws_of(stanfit$"sr_sigma")
    sr_L = posterior::draws_of(stanfit$"sr_L") 
    sr_raw = posterior::draws_of(stanfit$"sr_raw")

    s_mu = posterior::draws_of(stanfit$"s_mu") 
    s_sigma = posterior::draws_of(stanfit$"s_sigma")
    s_raw = posterior::draws_of(stanfit$"s_raw")

    c_mu = posterior::draws_of(stanfit$"c_mu") 
    c_sigma = posterior::draws_of(stanfit$"c_sigma")
    c_raw = posterior::draws_of(stanfit$"c_raw")

    dr_L = posterior::draws_of(stanfit$"dr_L")
    dr_raw = posterior::draws_of(stanfit$"dr_raw") 
    dr_sigma = posterior::draws_of(stanfit$"dr_sigma")
    
    if(dim(input$data$block_set)[2]>0)
    block_effects = posterior::draws_of(stanfit$"block_effects")

    if(dim(input$data$focal_set)[2]>1)
    focal_effects = posterior::draws_of(stanfit$"focal_effects") 

    if(dim(input$data$target_set)[2]>1)
    target_effects = posterior::draws_of(stanfit$"target_effects")  

    if(dim(input$data$sampling_set)[2]>1)
    sampling_effects = posterior::draws_of(stanfit$"sampling_effects")  

    if(dim(input$data$censoring_set)[2]>1)
    censoring_effects = posterior::draws_of(stanfit$"censoring_effects")  

    if(dim(input$data$dyad_set)[3]>1)
    dyad_effects = posterior::draws_of(stanfit$"dyad_effects")

    ################### Get index data for block-model samples
    block_indexes = c()
    block_indexes[1] = 0
    for(q in 1:input$data$N_group_vars){ 
    block_indexes[1+q] = input$data$N_groups_per_var[q]*input$data$N_groups_per_var[q] + block_indexes[q]
    }

    ################### Convert the block-model effects into an array form
    B = list()
    for(q in 1:input$data$N_group_vars){
      B[[q]] = array(NA, c(dim(block_effects)[1], input$data$N_groups_per_var[q], input$data$N_groups_per_var[q]  ))

      for(s in 1:dim(block_effects)[1]){
       B[[q]][s,,] =  array(block_effects[s,(block_indexes[q]+1):(block_indexes[q+1])], c(input$data$N_groups_per_var[q], input$data$N_groups_per_var[q]))
      }
    }   

    srm_samples = list(
            block_parameters=B,

            focal_target_sd=sr_sigma,
            focal_target_L=sr_L,
            focal_target_random_effects=sr_raw,

            dyadic_sd = dr_sigma,
            dyadic_L = dr_L,
            dyadic_random_effects=dr_raw,

            sampling_mu = s_mu,
            sampling_sd = s_sigma,
            sampling_random_effects = s_raw,

            censoring_mu = c_mu,
            censoring_sd = c_sigma,
            censoring_random_effects = c_raw
        )

    if(dim(input$data$focal_set)[2]>1)
    srm_samples$focal_coeffs = focal_effects

    if(dim(input$data$target_set)[2]>1)
    srm_samples$target_coeffs = target_effects

    if(dim(input$data$sampling_set)[2]>1)
    srm_samples$sampling_coeffs = sampling_effects

    if(dim(input$data$censoring_set)[2]>1)
    srm_samples$censoring_coeffs = censoring_effects

    if(dim(input$data$dyad_set)[3]>1)
    srm_samples$dyadic_coeffs = dyad_effects

    samples = list(srm_model_samples=srm_samples)

    if(input$return_predicted_network == TRUE){
        samples$predicted_network_sample = posterior::draws_of(stanfit$"p")
        }


    ###################################################### Create summary stats 
     sum_stats = function(y, x, z){
      bob = rep(NA, 6)
       dig = 3
      bob[1] = y
      bob[2] = round(median(x),dig)
      bob[3] = round(HPDI(x, z)[1],dig)
      bob[4] = round(HPDI(x, z)[2],dig)
      bob[5] = round(mean(x),dig)
      bob[6] = round(sd(x),dig)

      return(bob)
      }
     
     results_list = list()

    ################### SRM model
     Q1 = dim(input$data$focal_set)[2]-1
     Q2 = dim(input$data$target_set)[2]-1
     Q3 = dim(input$data$dyad_set)[3]-1
     Q4 = dim(input$data$sampling_set)[2]-1
     Q5 = dim(input$data$censoring_set)[2]-1

     results_srm_focal = matrix(NA, nrow=(1+Q1) , ncol=6)
     results_srm_target = matrix(NA, nrow=(1+Q2) , ncol=6)
     results_srm_sampling = matrix(NA, nrow=(2+Q4) , ncol=6)
     results_srm_censoring = matrix(NA, nrow=(2+Q5) , ncol=6)
     results_srm_dyadic = matrix(NA, nrow=(1+Q3) , ncol=6)

    ######### Calculate all focal effects
     results_srm_focal[1,] = sum_stats("focal effects sd", samples$srm_model_samples$focal_target_sd[,1], HPDI)
     if(Q1>0){
     coeff_names = colnames(input$data$focal_set)[-1]
        for(i in 1:Q1){
     results_srm_focal[1+i,] = sum_stats(paste0("focal effects coeffs (out-degree), ", coeff_names[i] ), samples$srm_model_samples$focal_coeffs[,i], HPDI)
        }
      }

      results_list[[1]] = results_srm_focal

    ######### Calculate all target effects
     results_srm_target[1,] = sum_stats("target effects sd", samples$srm_model_samples$focal_target_sd[,2], HPDI)
     if(Q2>0){
     coeff_names = colnames(input$data$target_set)[-1]
        for(i in 1:Q2){
     results_srm_target[1+i,] = sum_stats(paste0("target effects coeffs (in-degree), ", coeff_names[i] ), samples$srm_model_samples$target_coeffs[,i], HPDI)
        }
      }

      results_list[[2]] = results_srm_target

    ######### Calculate all sampling effects
     results_srm_sampling[1,] = sum_stats("sampling effects mu", samples$srm_model_samples$sampling_mu, HPDI)
     results_srm_sampling[2,] = sum_stats("sampling effects sd", samples$srm_model_samples$sampling_sd, HPDI)
     if(Q4>0){
     coeff_names = colnames(input$data$sampling_set)[-1]
        for(i in 1:Q4){
     results_srm_sampling[2+i,] = sum_stats(paste0("sampling effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$sampling_coeffs[,i], HPDI)
        }
      }

      results_list[[5]] = results_srm_sampling

    ######### Calculate all censoring effects
     results_srm_censoring[1,] = sum_stats("censoring effects mu", samples$srm_model_samples$censoring_mu, HPDI)
     results_srm_censoring[2,] = sum_stats("censoring effects sd", samples$srm_model_samples$censoring_sd, HPDI)
     if(Q4>0){
     coeff_names = colnames(input$data$censoring_set)[-1]
        for(i in 1:Q4){
     results_srm_censoring[2+i,] = sum_stats(paste0("censoring effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$censoring_coeffs[,i], HPDI)
        }
      }

      results_list[[6]] = results_srm_censoring

    ######### Calculate all dyad effects
     results_srm_dyadic[1,] = sum_stats("dyadic effects sd", c(samples$srm_model_samples$dyadic_sd), HPDI)
     if(Q3>0){
     coeff_names = dimnames(input$data$dyad_set)[[3]][-1]
        for(i in 1:Q3){
     results_srm_dyadic[1+i,] = sum_stats(paste0("dyadic effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$dyadic_coeffs[,i], HPDI)
        }
      }
     results_list[[3]] = results_srm_dyadic

    ######### Calculate all block effects
     results_srm_base = matrix(NA, nrow=2 + dim(block_effects)[2], ncol=6)
     results_srm_base[1,] = sum_stats("focal-target effects rho (generalized recipocity)", samples$srm_model_samples$focal_target_L[,2,1], HPDI)
     results_srm_base[2,] = sum_stats("dyadic effects rho (dyadic recipocity)", samples$srm_model_samples$dyadic_L[,2,1], HPDI)
 
     group_ids_character_df = cbind(rep("Any",input$data$N_id),attr(input$data, "group_ids_character"))

     if(is.null(colnames(group_ids_character_df))){
        colnames(group_ids_character_df) = paste0("(NoBlockingVars)", 1:ncol(group_ids_character_df))
     }
     
     colnames(group_ids_character_df)[1] = "(Intercept)"
     in_IDs = colnames(input$data$block_set)
     all_IDs = colnames(group_ids_character_df)
     group_ids_character_df = group_ids_character_df[,match(in_IDs, all_IDs), drop = FALSE]

     group_id_levels = append("Any", attr(input$data, "group_ids_levels"), 1)
     names(group_id_levels)[1]= "(Intercept)"
     
     ticker = 0
     for(q in 1:input$data$N_group_vars){
      group_ids_character = levels(factor(group_ids_character_df[,q]))
      test_sorting = group_id_levels[[which(names(group_id_levels) == colnames(group_ids_character_df)[q])]]
      if(all(group_ids_character==test_sorting)==FALSE){
        stop("Factors not sorted correctly.")
      }


      for(b1 in 1:input$data$N_groups_per_var[q]){
      for(b2 in 1:input$data$N_groups_per_var[q]){
       ticker = ticker + 1  
      results_srm_base[ 2+ ticker,] = sum_stats(paste0("offset, ", group_ids_character[b1], " to ", group_ids_character[b2]), 
                                                                         samples$srm_model_samples$block_parameters[[q]][,b1,b2], HPDI)
     }}

     }
     
     results_list[[4]] = results_srm_base

   ############# Finally, merge all effects into a list
     for(i in 1:6)
     colnames(results_list[[i]]) = c("Variable", "Median", "HPDI:L","HPDI:H","Mean","SD") 

     names(results_list) = c( "Focal effects: Out-degree", "Target effects: In-degree", "Dyadic effects", "Other estimates", "Sampling estimates", "Censoring estimates")
          
   results_out = rbind( results_srm_focal, results_srm_target,results_srm_dyadic, results_srm_base, results_srm_sampling, results_srm_censoring)
   
   df = data.frame(results_out)
   colnames(df) = c("Variable", "Median", "HPDI:L","HPDI:H","Mean","SD") 

   res_final = list(summary=df, summary_list=results_list)

  if(include_samples==TRUE){
    res_final$samples=samples
   }

   print(results_list)

    attr(res_final, "class") = "STRAND Results Object"
    return(res_final)
}



