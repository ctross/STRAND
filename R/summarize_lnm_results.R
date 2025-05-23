#' Organize Stan output and provide summaries of model parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters.
#'
#' @param input A STRAND model object, obtained by fitting a latent network model.
#' @param include_samples An indicator for the user to specify where raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A STRAND results object including summary table, a summary list, and samples.
#' @export
#' @examples
#' \dontrun{
#' res = summarize_lnm_results(input = fit)
#' }
#'

summarize_lnm_results = function(input, include_samples=TRUE, HPDI=0.9){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_lnm_results() requires a fitted object of class: STRAND Model Object. Please use fit_latent_network_model() to run your model.")
    }

    if(attributes(input)$fit_type != "mcmc"){
      if(attributes(input)$fit_type == "vb"){
         warning("Final, publication-ready model fits for STRAND models should always be produced using MCMC! Variational inference via Pathfinder can be used in Stan
              during experimental model runs, but final inferences should be based on MCMC sampling. In our tests, Pathfinder results are decently similar to MCMC results, 
              but often failed to recover strong true effects. ")  
         } else{
         stop("Fitted results can only be reorganized for STRAND model objects fit using MCMC. Variational inference or optimization can be used in Stan
              during experimental model runs, but final inferences should be based on MCMC sampling.")  
      }    
    }

    ###################################################### Create samples 
    stanfit = posterior::as_draws_rvars(input$fit$draws())

    ################### Measurement model
    false_positive_rate = posterior::draws_of(stanfit$"false_positive_rate")
    recall_of_true_ties = posterior::draws_of(stanfit$"recall_of_true_ties")
    theta_mean = posterior::draws_of(stanfit$"theta_mean")

    fpr_sigma = posterior::draws_of(stanfit$"fpr_sigma")
    rtt_sigma = posterior::draws_of(stanfit$"rtt_sigma")
    theta_sigma = posterior::draws_of(stanfit$"theta_sigma")

    fpr_raw = posterior::draws_of(stanfit$"fpr_raw")
    rtt_raw = posterior::draws_of(stanfit$"rtt_raw")
    theta_raw = posterior::draws_of(stanfit$"theta_raw")
     

    if(dim(input$data$fpr_set)[2]>1)
    fpr_effects = posterior::draws_of(stanfit$"fpr_effects")

    if(dim(input$data$rtt_set)[2]>1)
    rtt_effects = posterior::draws_of(stanfit$"rtt_effects")

    if(dim(input$data$theta_set)[2]>1)
    theta_effects = posterior::draws_of(stanfit$"theta_effects") 

    if(dim(input$data$block_set)[2]>0)
    block_effects = posterior::draws_of(stanfit$"block_effects") 

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

    ################### Network model
    sr_sigma = posterior::draws_of(stanfit$"sr_sigma")
    sr_L = posterior::draws_of(stanfit$"sr_L") 
    sr_raw = posterior::draws_of(stanfit$"sr_raw")

    dr_L = posterior::draws_of(stanfit$"dr_L")
    dr_raw = posterior::draws_of(stanfit$"dr_raw") 
    dr_sigma = posterior::draws_of(stanfit$"dr_sigma")
    
    if(dim(input$data$focal_set)[2]>1)
    focal_effects = posterior::draws_of(stanfit$"focal_effects")   
    if(dim(input$data$target_set)[2]>1)
    target_effects = posterior::draws_of(stanfit$"target_effects")  
    if(dim(input$data$dyad_set)[3]>1)
    dyad_effects = posterior::draws_of(stanfit$"dyad_effects") 

    measurement_samples = list(
            false_positive_rate_intercept=false_positive_rate,
            false_positive_rate_sd=fpr_sigma,
            false_positive_rate_random_effects=fpr_raw,
  
            recall_of_true_ties_intercept=recall_of_true_ties,
            recall_of_true_ties_sd=rtt_sigma,
            recall_of_true_ties_random_effects=rtt_raw,

            theta_intercept=theta_mean,
            theta_sd=theta_sigma,
            theta_random_effects=theta_raw
        )


    if(dim(input$data$fpr_set)[2]>1)
    measurement_samples$false_positive_rate_coeffs = fpr_effects

    if(dim(input$data$rtt_set)[2]>1)
    measurement_samples$recall_of_true_ties_coeffs = rtt_effects

    if(dim(input$data$theta_set)[2]>1)
    measurement_samples$theta_coeffs = theta_effects

    srm_samples = list(
            block_parameters=B,

            focal_target_sd=sr_sigma,
            focal_target_L=sr_L,
            focal_target_random_effects=sr_raw,

            dyadic_sd = dr_sigma,
            dyadic_L = dr_L,
            dyadic_random_effects=dr_raw
        )

    if(dim(input$data$focal_set)[2]>1)
    srm_samples$focal_coeffs = focal_effects

    if(dim(input$data$target_set)[2]>1)
    srm_samples$target_coeffs = target_effects

    if(dim(input$data$dyad_set)[3]>1)
    srm_samples$dyadic_coeffs = dyad_effects

    samples = list(measurement_model_samples=measurement_samples, srm_model_samples=srm_samples)

    if(input$return_predicted_network == TRUE){
        samples$predicted_network_sample = posterior::draws_of(stanfit$"p_tie_out") 
        }

    ###################################################### Create summary stats      
     results_list = list()
    ################### FPR model
     Q = dim(input$data$fpr_set)[2]-1
     results_fpr = matrix(NA, nrow=(4+(Q*2)), ncol=7)

     results_fpr[1,] = sum_stats("false positive rate intercept, layer 1", samples$measurement_model_samples$false_positive_rate_intercept[,1], HPDI)
     results_fpr[2,] = sum_stats("false positive rate intercept, layer 2", samples$measurement_model_samples$false_positive_rate_intercept[,2], HPDI)
     results_fpr[3,] = sum_stats("false positive rate sd, layer 1", samples$measurement_model_samples$false_positive_rate_sd[,1], HPDI)
     results_fpr[4,] = sum_stats("false positive rate sd, layer 2", samples$measurement_model_samples$false_positive_rate_sd[,2], HPDI)

     if(Q>0){
     coeff_names = colnames(input$data$fpr_set)[-1]
        for(i in 1:Q){
     results_fpr[4+i,] = sum_stats(paste0("false positive rate coeffs, layer 1, ", coeff_names[i] ), samples$measurement_model_samples$false_positive_rate_coeffs[,1,i], HPDI)
     results_fpr[4+i+Q,] = sum_stats(paste0("false positive rate coeffs, layer 2, ", coeff_names[i] ), samples$measurement_model_samples$false_positive_rate_coeffs[,2,i], HPDI)
        }
     }

     results_list[[1]] = results_fpr

    ################### RTT model
     Q = dim(input$data$rtt_set)[2]-1
     results_rtt = matrix(NA, nrow=(4+(Q*2)), ncol=7)

     results_rtt[1,] = sum_stats("recall rate of true ties intercept, layer 1", samples$measurement_model_samples$recall_of_true_ties_intercept[,1], HPDI)
     results_rtt[2,] = sum_stats("recall rate of true ties intercept, layer 2", samples$measurement_model_samples$recall_of_true_ties_intercept[,2], HPDI)
     results_rtt[3,] = sum_stats("recall rate of true ties sd, layer 1", samples$measurement_model_samples$recall_of_true_ties_sd[,1], HPDI)
     results_rtt[4,] = sum_stats("recall rate of true ties sd, layer 2", samples$measurement_model_samples$recall_of_true_ties_sd[,2], HPDI)

     if(Q>0){
     coeff_names = colnames(input$data$rtt_set)[-1]
        for(i in 1:Q){
     results_rtt[4+i,] = sum_stats(paste0("recall rate of true ties coeffs, layer 1, ", coeff_names[i] ), samples$measurement_model_samples$recall_of_true_ties_coeffs[,1,i], HPDI)
     results_rtt[4+i+Q,] = sum_stats(paste0("recall rate of true ties coeffs, layer 2, ", coeff_names[i] ), samples$measurement_model_samples$recall_of_true_ties_coeffs[,2,i], HPDI)
        }
     }

     results_list[[2]] = results_rtt

    ################### Theta model
     Q = dim(input$data$theta_set)[2]-1
     results_theta = matrix(NA, nrow=(2+(Q)), ncol=7)

     results_theta[1,] = sum_stats("theta intercept, layer 1 to 2", c(samples$measurement_model_samples$theta_intercept), HPDI)
     results_theta[2,] = sum_stats("theta sd, layer 1 to 2", c(samples$measurement_model_samples$theta_sd), HPDI)

     if(Q>0){
     coeff_names = colnames(input$data$theta_set)[-1]
        for(i in 1:Q){
     results_theta[2+i,] = sum_stats(paste0("theta coeffs, layer 1 to 2, ", coeff_names[i] ), samples$measurement_model_samples$theta_coeffs[,i], HPDI)
        }
     }
   
    results_list[[3]] = results_theta

    measurement_results = rbind(results_fpr, results_rtt, results_theta)
 
    ################### SRM model
    Q1 = dim(input$data$focal_set)[2]-1
    Q2 = dim(input$data$target_set)[2]-1
    Q3 = dim(input$data$dyad_set)[3]-1

     results_srm_focal = matrix(NA, nrow=(1+Q1) , ncol=7)
     results_srm_target = matrix(NA, nrow=(1+Q2) , ncol=7)
     results_srm_dyadic = matrix(NA, nrow=(1+Q3) , ncol=7)

     ######### Calculate all focal effects
     results_srm_focal[1,] = sum_stats("focal effects sd", samples$srm_model_samples$focal_target_sd[,1], HPDI)
     if(Q1>0){
     coeff_names = colnames(input$data$focal_set)[-1]
        for(i in 1:Q1){
     results_srm_focal[1+i,] = sum_stats(paste0("focal effects coeffs (out-degree), ", coeff_names[i] ), samples$srm_model_samples$focal_coeffs[,i], HPDI)
        }
      }

      results_list[[4]] = results_srm_focal

     ######### Calculate all target effects
     results_srm_target[1,] = sum_stats("target effects sd", samples$srm_model_samples$focal_target_sd[,2], HPDI)
     if(Q2>0){
     coeff_names = colnames(input$data$target_set)[-1]
        for(i in 1:Q2){
     results_srm_target[1+i,] = sum_stats(paste0("target effects coeffs (in-degree), ", coeff_names[i] ), samples$srm_model_samples$target_coeffs[,i], HPDI)
        }
      }

      results_list[[5]] = results_srm_target
      
     ######### Calculate all dyad effects
     results_srm_dyadic[1,] = sum_stats("dyadic effects sd", c(samples$srm_model_samples$dyadic_sd), HPDI)
     if(Q3>0){
     coeff_names = dimnames(input$data$dyad_set)[[3]][-1]
        for(i in 1:Q3){
     results_srm_dyadic[1+i,] = sum_stats(paste0("dyadic effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$dyadic_coeffs[,i], HPDI)
        }
      }
     results_list[[6]] = results_srm_dyadic

     ######### Calculate all other effects
     results_srm_base = matrix(NA, nrow=2 + dim(block_effects)[2], ncol=7)
     results_srm_base[1,] = sum_stats("focal-target effects rho (generalized recipocity)", samples$srm_model_samples$focal_target_L[,2,1], HPDI)
     results_srm_base[2,] = sum_stats("dyadic effects rho (dyadic recipocity)", samples$srm_model_samples$dyadic_L[,2,1], HPDI)

     group_ids_character_df = cbind(rep("Any",input$data$N_id),attr(input$data, "group_ids_character"))

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
     
     results_list[[7]] = results_srm_base

  ############# Finally, merge all effects into a list

     for(i in 1:7)
     colnames(results_list[[i]]) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 

     names(results_list) = c("False positive rate", "Recall of true ties", "Theta: question-order effects", "Focal effects: Out-degree", "Target effects: In-degree", "Dyadic effects", "Other estimates")
          
   results_out = rbind(measurement_results, results_srm_focal, results_srm_target,results_srm_dyadic, results_srm_base)
   
   df = data.frame(results_out)
   colnames(df) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 


   res_final = list(summary=df, summary_list=results_list)

  if(include_samples==TRUE){
    res_final$samples=samples
   }

   print(results_list)

    attr(res_final, "class") = "STRAND Results Object"
    return(res_final)
}
