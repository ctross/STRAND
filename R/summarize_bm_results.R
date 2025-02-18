#' Organize Stan output and provide summaries of model parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters.
#'
#' @param input A STRAND model object, obtained by fitting a stochastic block model.
#' @param include_samples An indicator for the user to specify whether raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A STRAND results object including summary table, a summary list, and samples.
#' @export
#' @examples
#' \dontrun{
#' res = summarize_bm_results(input = fit)
#' }
#'

summarize_bm_results = function(input, include_samples=TRUE, HPDI=0.9){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_bm_results() requires a fitted object of class: STRAND Model Object. Please use fit_block_model() to run your model.")
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

    outcome_mode = input$data$outcome_mode 

    ################### Extract effects from sample set
    if(dim(input$data$block_set)[2]>0)
    block_effects = posterior::draws_of(stanfit$"block_effects")

    if(dim(input$data$focal_set)[2]>1)
    focal_effects = posterior::draws_of(stanfit$"focal_effects")

    if(dim(input$data$target_set)[2]>1)
    target_effects = posterior::draws_of(stanfit$"target_effects") 

    if(dim(input$data$dyad_set)[3]>1)
    dyad_effects = posterior::draws_of(stanfit$"dyad_effects")

    error_sigma = posterior::draws_of(stanfit$"error_sigma")

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
     
    ################## Now build up full list of samples
    srm_samples = list(
       block_parameters = B,
       error_sd = error_sigma
        )

    if(dim(input$data$focal_set)[2]>1)
    srm_samples$focal_coeffs = focal_effects

    if(dim(input$data$target_set)[2]>1)
    srm_samples$target_coeffs = target_effects

    if(dim(input$data$dyad_set)[3]>1)
    srm_samples$dyadic_coeffs = dyad_effects

    samples = list(srm_model_samples=srm_samples)

    if(input$return_predicted_network == TRUE){
        samples$predicted_network_sample = posterior::draws_of(stanfit$"p")
        }

    ###################################################### Create summary stats      
     ######### Prepare results arrays
     results_list = list()

     Q1 = dim(input$data$focal_set)[2]-1
     Q2 = dim(input$data$target_set)[2]-1
     Q3 = dim(input$data$dyad_set)[3]-1

     results_srm_focal = matrix(NA, nrow=(1+Q1) , ncol=7)
     results_srm_target = matrix(NA, nrow=(1+Q2) , ncol=7)
     results_srm_dyadic = matrix(NA, nrow=(1+Q3) , ncol=7)

     ######### Calculate all focal effects
     results_srm_focal[1,1] = "focal effects sd"
     if(Q1>0){
     coeff_names = colnames(input$data$focal_set)[-1]
        for(i in 1:Q1){
     results_srm_focal[1+i,] = sum_stats(paste0("focal effects coeffs (out-degree), ", coeff_names[i] ), samples$srm_model_samples$focal_coeffs[,i], HPDI)
        }
      }

      results_list[[1]] = results_srm_focal

     ######### Calculate all target effects
     results_srm_target[1,1] = "target effects sd"
     if(Q2>0){
     coeff_names = colnames(input$data$target_set)[-1]
        for(i in 1:Q2){
     results_srm_target[1+i,] = sum_stats(paste0("target effects coeffs (in-degree), ", coeff_names[i] ), samples$srm_model_samples$target_coeffs[,i], HPDI)
        }
      }

      results_list[[2]] = results_srm_target

     ######### Calculate all dyadic effects
     results_srm_dyadic[1,1] = "dyadic effects sd"
     if(Q3>0){
     coeff_names = dimnames(input$data$dyad_set)[[3]][-1]
        for(i in 1:Q3){
     results_srm_dyadic[1+i,] = sum_stats(paste0("dyadic effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$dyadic_coeffs[,i], HPDI)
        }
      }

     results_list[[3]] = results_srm_dyadic

     ######### Calculate all block-model effects
     results_srm_base = matrix(NA, nrow = (1 + dim(block_effects)[2]), ncol=7)

        if(outcome_mode == 4){
        results_srm_base[1,] = sum_stats("error sd", c(samples$srm_model_samples$error_sd), HPDI)
        } else{
        results_srm_base[1,] = c("error sd - ", rep(NA,6))     
        }

     
 
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
      results_srm_base[ 1 + ticker,] = sum_stats(paste0("offset, ", group_ids_character[b1], " to ", group_ids_character[b2]), 
                                                                         samples$srm_model_samples$block_parameters[[q]][,b1,b2], HPDI)
     }}

     }

     results_list[[4]] = results_srm_base
    
    ############# Finally, merge all effects into a list
     for(i in 1:4)
     colnames(results_list[[i]]) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 

     names(results_list) = c( "Focal effects: Out-degree", "Target effects: In-degree", "Dyadic effects", "Other estimates")
          
     results_out = rbind( results_srm_focal, results_srm_target,results_srm_dyadic, results_srm_base)
   
     df = data.frame(results_out)
     colnames(df) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 

     df = df[complete.cases(df),]

     res_final = list(summary=df, summary_list=results_list)

     if(include_samples==TRUE){
       res_final$samples=samples
       }

    print(results_list)

    attr(res_final, "class") = "STRAND Results Object"
    return(res_final)
}
