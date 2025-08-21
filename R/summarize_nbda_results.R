#' Organize Stan output and provide summaries of model parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters.
#'
#' @param input A STRAND model object, obtained by fitting an NBDA model.
#' @param include_samples An indicator for the user to specify where raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A STRAND results object including summary table, a summary list, and samples.
#' @export
#' @examples
#' \dontrun{
#' res = summarize_nbda_results(input = fit)
#' }
#'

summarize_nbda_results = function(input, include_samples=TRUE, HPDI=0.9){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_nbda_results() requires a fitted object of class: STRAND Model Object. Please use fit_NBDA_model() to run your model.")
    }

    if(attributes(input)$model_type != "NBDA"){
        stop("summarize_nbda_results() requires a fitted object of class: NBDA. Please use fit_NBDA_model() to run your model.")
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

    ################### Network model
    alpha = posterior::draws_of(stanfit$"alpha")
    sigma = posterior::draws_of(stanfit$"sigma") 
    eta = posterior::draws_of(stanfit$"eta")
    lambda = posterior::draws_of(stanfit$"lambda")

    if(dim(input$data$ind_focal_set)[2]>1)
    ind_focal_effects = posterior::draws_of(stanfit$"ind_focal_effects") 
    
    if(dim(input$data$block_set)[2]>0)
    block_effects = posterior::draws_of(stanfit$"block_effects")

    if(dim(input$data$focal_set)[2]>1)
    focal_effects = posterior::draws_of(stanfit$"focal_effects") 

    if(dim(input$data$target_set)[2]>1)
    target_effects = posterior::draws_of(stanfit$"target_effects")  

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
            block_parameters = B,
            lambda = lambda,

            alpha = alpha,
            sigma = sigma,
            eta = eta
        )

    if(dim(input$data$ind_focal_set)[2]>1)
    srm_samples$ind_focal_coeffs = ind_focal_effects

    if(dim(input$data$focal_set)[2]>1)
    srm_samples$focal_coeffs = focal_effects

    if(dim(input$data$target_set)[2]>1)
    srm_samples$target_coeffs = target_effects

    if(dim(input$data$dyad_set)[3]>1)
    srm_samples$dyadic_coeffs = dyad_effects

    samples = list(srm_model_samples=srm_samples)


    ###################################################### Create summary stats   
     results_list = list()

    ################### SRM model
    Q0 = dim(input$data$ind_focal_set)[2]-1
    Q1 = dim(input$data$focal_set)[2]-1
    Q2 = dim(input$data$target_set)[2]-1
    Q3 = dim(input$data$dyad_set)[3]-1

     results_srm_ind_focal = matrix(NA, nrow=((Q0)+1) , ncol=7)
     results_srm_focal = matrix(NA, nrow=(Q1) , ncol=7)
     results_srm_target = matrix(NA, nrow=(Q2) , ncol=7)
     results_srm_dyadic = matrix(NA, nrow=(Q3) , ncol=7)

     results_srm_ind_focal[1,] = sum_stats("Individual learning propensity (focal), intercept", samples$srm_model_samples$lambda, HPDI)

     if(Q0>0){
     coeff_names = colnames(input$data$ind_focal_set)[-1]
        for(i in 1:Q0){
     results_srm_ind_focal[i+1,] = sum_stats(paste0("Individual learning propensity (focal), ", coeff_names[i] ), samples$srm_model_samples$ind_focal_coeffs[,i], HPDI)
        }
      }

      results_list[[1]] = results_srm_ind_focal


     if(Q1>0){
     coeff_names = colnames(input$data$focal_set)[-1]
        for(i in 1:Q1){
     results_srm_focal[i,] = sum_stats(paste0("Social learning propensity (focal), ", coeff_names[i] ), samples$srm_model_samples$focal_coeffs[,i], HPDI)
        }
      }

      results_list[[2]] = results_srm_focal


     if(Q2>0){
     coeff_names = colnames(input$data$target_set)[-1]
        for(i in 1:Q2){
     results_srm_target[i,] = sum_stats(paste0("Social learning propensity (target), ", coeff_names[i] ), samples$srm_model_samples$target_coeffs[,i], HPDI)
        }
      }

      results_list[[3]] = results_srm_target


     if(Q3>0){
     coeff_names = dimnames(input$data$dyad_set)[[3]][-1]
        for(i in 1:Q3){
     results_srm_dyadic[i,] = sum_stats(paste0("Social learning propensity (dyadic), ", coeff_names[i] ), samples$srm_model_samples$dyadic_coeffs[,i], HPDI)
        }
      }

     results_list[[4]] = results_srm_dyadic

      ######### Calculate all block effects
     results_srm_base = matrix(NA, nrow=dim(block_effects)[2], ncol=7)
 
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
      results_srm_base[ticker,] = sum_stats(paste0("Social learning propensity (block), ", group_ids_character[b1], " to ", group_ids_character[b2]), 
                                                                         samples$srm_model_samples$block_parameters[[q]][,b1,b2], HPDI)
     }}

     }
     
     results_list[[5]] = results_srm_base


     results_srm_base_extra = matrix(NA, nrow=3, ncol=7)
     
     if(input$data$ces_settings == 1){
     results_srm_base_extra[1,] = sum_stats("Out-to-In share (CES)", samples$srm_model_samples$alpha, HPDI)
     results_srm_base_extra[2,] = sum_stats("Elasticity of substitution (CES)", samples$srm_model_samples$sigma, HPDI)
     results_srm_base_extra[3,] = sum_stats("Returns to scale (CES)", samples$srm_model_samples$eta, HPDI)
     }

    if(input$data$ces_settings == 2){
     results_srm_base_extra[1,] = sum_stats("Out-to-In share (CES)", samples$srm_model_samples$alpha, HPDI)
     results_srm_base_extra[2,] = sum_stats("Elasticity of substitution (CES)", samples$srm_model_samples$sigma, HPDI)
     results_srm_base_extra[3,1] = "Returns to scale (CES)"
     }

    if(input$data$ces_settings == 3){
     results_srm_base_extra[1,] = sum_stats("Out-to-In share (CES)", samples$srm_model_samples$alpha, HPDI)
     results_srm_base_extra[2,1] = "Elasticity of substitution (CES)"
     results_srm_base_extra[3,1] = "Returns to scale (CES)"
     }
     
     results_list[[6]] = results_srm_base_extra


     for(i in 1:6)
     colnames(results_list[[i]]) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 

     names(results_list) = c( "Individual learning propensity", "Social learning propensity: Focal", "Social learning propensity: Target", 
                              "Social learning propensity: Dyadic", "Social learning propensity: Block", "CES function parameters")
          
   results_out = rbind( results_srm_focal, results_srm_target,results_srm_dyadic, results_srm_base)
   
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
