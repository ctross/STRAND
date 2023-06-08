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




summarize_bsrm_hh_results2 = function(input, include_samples=TRUE, HPDI=0.9){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_bsrm_hh_results() requires a fitted object of class: STRAND Model Object. Please use fit_block_plus_social_relations_hh_model() to run your model.")
    }

    if(attributes(input)$fit_type != "mcmc"){
        stop("Fitted results can only be reorganized for STRAND model objects fit using MCMC. Variational inference or optimization can be used in Stan
              during experimental model runs, but final inferences should be based on MCMC sampling.")   
    }

    ###################################################### Create samples 
    stanfit = posterior::as_draws_rvars(input$fit$draws())
    
    ################### Network model parameters - Indiv
    sr_sigma = draws_of(stanfit$"sr_sigma")
    sr_L = draws_of(stanfit$"sr_L") 
    sr_raw = draws_of(stanfit$"sr_raw")

    dr_L = draws_of(stanfit$"dr_L")
    dr_raw = draws_of(stanfit$"dr_raw") 
    dr_sigma = draws_of(stanfit$"dr_sigma")

    if(dim(input$data$block_set)[2]>0)
    block_effects = draws_of(stanfit$"block_effects")

    if(dim(input$data$focal_set)[2]>1)
    focal_effects = draws_of(stanfit$"focal_effects")

    if(dim(input$data$target_set)[2]>1)
    target_effects = draws_of(stanfit$"target_effects") 

    if(dim(input$data$dyad_set)[3]>1)
    dyad_effects = draws_of(stanfit$"dyad_effects")

    ################### Network model parameters - HH
    hh_sr_sigma = draws_of(stanfit$"hh_sr_sigma")
    hh_sr_L = draws_of(stanfit$"hh_sr_L") 
    hh_sr_raw = draws_of(stanfit$"hh_sr_raw") 
    hh_within_mu = draws_of(stanfit$"hh_within_mu")

    hh_dr_sigma = draws_of(stanfit$"hh_dr_sigma")
    hh_dr_L = draws_of(stanfit$"hh_dr_L") 
    hh_dr_raw = draws_of(stanfit$"hh_dr_raw") 
    
    if(dim(input$data$hh_focal_set)[2]>1)
    hh_focal_effects = draws_of(stanfit$"hh_focal_effects")

    if(dim(input$data$hh_target_set)[2]>1)
    hh_target_effects = draws_of(stanfit$"hh_target_effects")

    if(dim(input$data$hh_within_set)[2]>1)
    hh_within_effects = draws_of(stanfit$"hh_within_effects")

    if(dim(input$data$hh_between_set)[3]>1)
    hh_between_effects = draws_of(stanfit$"hh_between_effects")


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

            hh_focal_target_sd=hh_sr_sigma,
            hh_focal_target_L=hh_sr_L,
            hh_focal_target_random_effects=hh_sr_raw,
            hh_within_mu = hh_within_mu,

            hh_dyadic_sd = hh_dr_sigma,
            hh_dyadic_L = hh_dr_L,
            hh_dyadic_random_effects = hh_dr_raw
        )

    if(attributes(input)$model_version == "ulre"){
      srm_samples$dyadic_L = dr_L
      srm_samples$dyadic_random_effects=dr_raw
      srm_samples$dyadic_sd = dr_sigma
     } 

    #### indiv covars
    if(dim(input$data$focal_set)[2]>1)
    srm_samples$focal_coeffs = focal_effects

    if(dim(input$data$target_set)[2]>1)
    srm_samples$target_coeffs = target_effects

    if(dim(input$data$dyad_set)[3]>1)
    srm_samples$dyadic_coeffs = dyad_effects

    #### hh covars
    if(dim(input$data$hh_focal_set)[2]>1)
    srm_samples$hh_focal_coeffs = hh_focal_effects

    if(dim(input$data$hh_target_set)[2]>1)
    srm_samples$hh_target_coeffs = hh_target_effects

    if(dim(input$data$hh_within_set)[2]>1)
    srm_samples$hh_within_coeffs = hh_within_effects

    if(dim(input$data$hh_between_set)[3]>1)
    srm_samples$hh_between_coeffs = hh_between_effects

    samples = list(srm_model_samples=srm_samples)

    if(input$return_predicted_network == TRUE){
        samples$predicted_network_sample = draws_of(stanfit$"p")
        samples$predicted_hh_sample = draws_of(stanfit$"hh_dr")
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

     sum_stats_miss = function(y, x, z){
      bob = rep(NA, 6)
       dig = 3
      bob[1] = y
   
      return(bob)
      }
     
     results_list = list()
    
    ############################################################# Indiv
    ################### SRM model
     Q1 = dim(input$data$focal_set)[2]-1
     Q2 = dim(input$data$target_set)[2]-1
     Q3 = dim(input$data$dyad_set)[3]-1

     results_srm_focal = matrix(NA, nrow=(1+Q1) , ncol=6)
     results_srm_target = matrix(NA, nrow=(1+Q2) , ncol=6)
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

    ######### Calculate all dyad effects
     if(attributes(input)$model_version == "ulre"){
       results_srm_dyadic[1,] = sum_stats("dyadic effects sd", c(samples$srm_model_samples$dyadic_sd), HPDI)
     } 
     if(attributes(input)$model_version == "fast_bb"){
       results_srm_dyadic[1,] = sum_stats("dyadic effects cross-ratio", c(samples$srm_model_samples$dyadic_cross_ratio), HPDI)
     }
     if(attributes(input)$model_version == "no_dr"){
       results_srm_dyadic[1,] = sum_stats_miss("dyadic effects sd", c(0), HPDI)
     }  

     
     if(Q3>0){
     coeff_names = dimnames(input$data$dyad_set)[[3]][-1]
        for(i in 1:Q3){
     results_srm_dyadic[1+i,] = sum_stats(paste0("dyadic effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$dyadic_coeffs[,i], HPDI)
        }
      }
     results_list[[3]] = results_srm_dyadic

    ######### Calculate all block effects
     results_srm_base = matrix(NA, nrow=2 + dim(block_effects)[2], ncol=6)


     if(attributes(input)$model_version == "ulre"){
       results_srm_base[1,] = sum_stats("focal-target effects rho (generalized reciprocity)", samples$srm_model_samples$focal_target_L[,2,1], HPDI)
       results_srm_base[2,] = sum_stats("dyadic effects rho (dyadic reciprocity)", samples$srm_model_samples$dyadic_L[,2,1], HPDI)
     } 

     if(attributes(input)$model_version == "no_dr"){
       results_srm_base[1,] = sum_stats("focal-target effects rho (generalized reciprocity)", samples$srm_model_samples$focal_target_L[,2,1], HPDI)
       results_srm_base[2,] = sum_stats_miss("dyadic effects rho (dyadic reciprocity)", c(0), HPDI)
     } 
 

     group_ids_character_df = cbind(rep("Any",input$data$N_id),attr(input$data, "group_ids_character"))
     
     colnames(group_ids_character_df)[1] = "(Intercept)"
     in_IDs = colnames(input$data$block_set)
     all_IDs = colnames(group_ids_character_df)
     group_ids_character_df = group_ids_character_df[,match(in_IDs, all_IDs)]

     group_id_levels = append("Any", attr(input$data, "group_ids_levels"), 1)
     
     ticker = 0
     for(q in 1:input$data$N_group_vars){
      group_ids_character = group_id_levels[[q]]

      for(b1 in 1:input$data$N_groups_per_var[q]){
      for(b2 in 1:input$data$N_groups_per_var[q]){
       ticker = ticker + 1  
      results_srm_base[ 2+ ticker,] = sum_stats(paste0("offset, ", group_ids_character[b1], " to ", group_ids_character[b2]), 
                                                                         samples$srm_model_samples$block_parameters[[q]][,b1,b2], HPDI)
     }}

     }
     
     results_list[[4]] = results_srm_base

    ############################################################# HH
    ################### SRM model
     hh_Q1 = dim(input$data$hh_focal_set)[2]-1
     hh_Q2 = dim(input$data$hh_target_set)[2]-1
     hh_Q3 = dim(input$data$hh_within_set)[2]-1
     hh_Q4 = dim(input$data$hh_between_set)[3]-1

     hh_results_srm_focal = matrix(NA, nrow=(1+hh_Q1) , ncol=6)
     hh_results_srm_target = matrix(NA, nrow=(1+hh_Q2) , ncol=6)
     hh_results_srm_within = matrix(NA, nrow=(2+hh_Q3) , ncol=6)
     hh_results_srm_between = matrix(NA, nrow=(1+hh_Q4) , ncol=6)

    ######### Calculate all hh focal effects
     hh_results_srm_focal[1,] = sum_stats("household focal effects sd", samples$srm_model_samples$hh_focal_target_sd[,1], HPDI)
     if(hh_Q1>0){
     coeff_names = colnames(input$data$hh_focal_set)[-1]
        for(i in 1:hh_Q1){
     hh_results_srm_focal[1+i,] = sum_stats(paste0("household focal effects coeffs (out-degree), ", coeff_names[i] ), samples$srm_model_samples$hh_focal_coeffs[,i], HPDI)
        }
      }

      results_list[[5]] = hh_results_srm_focal

    ######### Calculate all hh target effects
     hh_results_srm_target[1,] = sum_stats("household target effects sd", samples$srm_model_samples$hh_focal_target_sd[,2], HPDI)
     if(hh_Q2>0){
     coeff_names = colnames(input$data$hh_target_set)[-1]
        for(i in 1:hh_Q2){
     hh_results_srm_target[1+i,] = sum_stats(paste0("household target effects coeffs (in-degree), ", coeff_names[i] ), samples$srm_model_samples$hh_target_coeffs[,i], HPDI)
        }
      }

      results_list[[6]] = hh_results_srm_target

    ######### Calculate all hh within effects
     hh_results_srm_within[1,] = sum_stats("household within effects sd", samples$srm_model_samples$hh_focal_target_sd[,3], HPDI)
     hh_results_srm_within[2,] = sum_stats("household within effects mu", c(samples$srm_model_samples$hh_within_mu), HPDI)

     if(hh_Q3>0){
     coeff_names = colnames(input$data$hh_within_set)[-1]
        for(i in 1:hh_Q3){
     hh_results_srm_within[2+i,] = sum_stats(paste0("household within effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$hh_within_coeffs[,i], HPDI)
        }
      }

      results_list[[7]] = hh_results_srm_within

    ######### Calculate all hh dyad effects
     hh_results_srm_between[1,] = sum_stats("household between effects sd", c(samples$srm_model_samples$hh_dyadic_sd), HPDI)
     if(hh_Q4>0){
     coeff_names = dimnames(input$data$hh_between_set)[[3]][-1]
        for(i in 1:hh_Q4){
     hh_results_srm_between[1+i,] = sum_stats(paste0("household between effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$hh_between_coeffs[,i], HPDI)
        }
      }
     results_list[[8]] = hh_results_srm_between

    ######## Other hh effects

    samples$srm_model_samples$hh_focal_target_Corr = samples$srm_model_samples$hh_focal_target_L

    for(i in 1: dim(samples$srm_model_samples$hh_focal_target_Corr)[1]){
     samples$srm_model_samples$hh_focal_target_Corr[i,,] = samples$srm_model_samples$hh_focal_target_L[i,,] %*% t(samples$srm_model_samples$hh_focal_target_L[i,,])  
    }

     hh_results_srm_base = matrix(NA, nrow=4, ncol=6)
     hh_results_srm_base[1,] = sum_stats("household focal-target effects rho (generalized reciprocity)", samples$srm_model_samples$hh_focal_target_Corr[,2,1], HPDI)
     hh_results_srm_base[2,] = sum_stats("household focal-within effects rho", samples$srm_model_samples$hh_focal_target_Corr[,3,1], HPDI)
     hh_results_srm_base[3,] = sum_stats("household target-within effects rho", samples$srm_model_samples$hh_focal_target_Corr[,3,2], HPDI)
     hh_results_srm_base[4,] = sum_stats("household dyadic effects rho (dyadic reciprocity)", samples$srm_model_samples$hh_dyadic_L[,2,1], HPDI)

                                                                    
     results_list[[9]] = hh_results_srm_base

   ############# Finally, merge all effects into a list
     for(i in 1:9)
     colnames(results_list[[i]]) = c("Variable", "Median", "HPDI:L","HPDI:H","Mean","SD") 

     names(results_list) = c(c( "Focal effects: Out-degree", "Target effects: In-degree", "Dyadic effects", "Other estimates"),
                             c( "Household: Focal effects/Out-degree", "Household: Target effects/In-degree", "Household: Within effects","Household: Between effects", "Household: Other estimates")
                             )
          
   results_out = rbind( results_srm_focal, results_srm_target, results_srm_dyadic, results_srm_base,
                        hh_results_srm_focal, hh_results_srm_target, hh_results_srm_within, hh_results_srm_between, hh_results_srm_base)
   
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


