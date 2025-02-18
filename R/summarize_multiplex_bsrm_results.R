#' Organize Stan output and provide summaries of model parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters
#'
#' @param input A STRAND model object, obtained by fitting a multiplex combined stochastic block and social relations model.
#' @param include_samples An indicator for the user to specify where raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A STRAND results object including summary table, a summary list, and samples.
#' @export
#' @examples
#' \dontrun{
#' res = summarize_multiplex_bsrm_results(input=fit)
#' }
#'

summarize_multiplex_bsrm_results = function(input, include_samples=TRUE, HPDI=0.9){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_multiplex_bsrm_results() requires a fitted object of class: STRAND Model Object.")
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
    N_responses = input$data$N_responses
    layer_names = attr(input$data,"layer_names")

    outcome_mode = input$data$outcome_mode 

    stanfit = posterior::as_draws_rvars(input$fit$draws())

    ################### Network model parameters
    sr_sigma = posterior::draws_of(stanfit$"sr_sigma")
    sr_L = posterior::draws_of(stanfit$"sr_L") 
    sr_raw = posterior::draws_of(stanfit$"sr_raw")

    dr_L = posterior::draws_of(stanfit$"dr_L")
    dr_raw = posterior::draws_of(stanfit$"dr_raw") 
    dr_sigma = posterior::draws_of(stanfit$"dr_sigma")
    error_sigma = posterior::draws_of(stanfit$"error_sigma")

    G_corr = posterior::draws_of(stanfit$"G_corr")
    D_corr = posterior::draws_of(stanfit$"D_corr")
    
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
    B_multi = list()
    for(l in 1:N_responses){
    B = list()
    for(q in 1:input$data$N_group_vars){
      B[[q]] = array(NA, c(dim(block_effects)[1], input$data$N_groups_per_var[q], input$data$N_groups_per_var[q]  ))

      for(s in 1:dim(block_effects)[1]){
       B[[q]][s,,] =  array(block_effects[s,l,(block_indexes[q]+1):(block_indexes[q+1])], c(input$data$N_groups_per_var[q], input$data$N_groups_per_var[q]))
      }
    } 
    B_multi[[l]] = B  
    }

    names(B_multi) = layer_names

    srm_samples = list(
            block_parameters=B_multi,

            focal_target_sd=sr_sigma,
            focal_target_L=sr_L,
            focal_target_random_effects=sr_raw,
            
            error_sd = error_sigma,
            dyadic_sd = dr_sigma,
            dyadic_L = dr_L,
            dyadic_random_effects=dr_raw,

            G_corr=G_corr,
            D_corr=D_corr
        )

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
     Q1 = dim(input$data$focal_set)[2]-1
     Q2 = dim(input$data$target_set)[2]-1
     Q3 = dim(input$data$dyad_set)[3]-1

     results_srm_focal = matrix(NA, nrow=(N_responses + Q1*N_responses), ncol=7)
     results_srm_target = matrix(NA, nrow=(N_responses + Q2*N_responses), ncol=7)
     results_srm_dyadic = matrix(NA, nrow=(N_responses + Q3*N_responses), ncol=7)

    ######### Calculate all focal effects
     for(q in 1:N_responses){
      results_srm_focal[q,] = sum_stats(paste0("focal effects sd - ", layer_names[q]), samples$srm_model_samples$focal_target_sd[,q], HPDI)
     }

     if(Q1>0){
     coeff_names = colnames(input$data$focal_set)[-1]
        for(q in 1:N_responses){
        for(i in 1:Q1){
     results_srm_focal[N_responses + i + Q1*(q-1),] = sum_stats(paste0("focal effects coeffs (out-degree), ", layer_names[q], " - ", coeff_names[i] ), samples$srm_model_samples$focal_coeffs[,q,i], HPDI)
        }}
      }

      results_list[[1]] = results_srm_focal

    ######### Calculate all target effects
     for(q in 1:N_responses){
      results_srm_target[q,] = sum_stats(paste0("target effects sd - ", layer_names[q]), samples$srm_model_samples$focal_target_sd[,N_responses+q], HPDI)
     }

     if(Q2>0){
     coeff_names = colnames(input$data$target_set)[-1]
        for(q in 1:N_responses){
        for(i in 1:Q2){
     results_srm_target[N_responses + i + Q2*(q-1),] = sum_stats(paste0("target effects coeffs (in-degree), ", layer_names[q], " - ", coeff_names[i] ), samples$srm_model_samples$target_coeffs[,q,i], HPDI)
        }}
      }

      results_list[[2]] = results_srm_target

    ######### Calculate all dyad effects
     for(q in 1:N_responses){
      results_srm_dyadic[q,] = sum_stats(paste0("dyadic effects sd - ", layer_names[q]), samples$srm_model_samples$dyadic_sd[,q], HPDI)
     }

     if(Q3>0){
     coeff_names = dimnames(input$data$dyad_set)[[3]][-1]
        for(q in 1:N_responses){
        for(i in 1:Q3){
     results_srm_dyadic[N_responses + i + Q3*(q-1),] = sum_stats(paste0("dyadic effects coeffs, ", layer_names[q], " - ", coeff_names[i] ), samples$srm_model_samples$dyadic_coeffs[,q,i], HPDI)
        }}
      }


     results_list[[3]] = results_srm_dyadic

    ######### Calculate all corr and block effects
    layer_names_long = c(paste0(layer_names, " (out)"), paste0(layer_names, " (in)"))
    layer_names_long2 = c(paste0(layer_names, " (i to j)"), paste0(layer_names, " (j to i)"))
     results_srm_base = matrix(NA, nrow = (1+ (2*N_responses*(2*N_responses-1)) + N_responses*dim(block_effects)[3]), ncol=7)
     ticker = 1

     if(outcome_mode == 4){
        results_srm_base[ticker,] = sum_stats(paste0("error sd - ", layer_names[q]), samples$srm_model_samples$error_sd[,q], HPDI)
        } else{
        results_srm_base[ticker,] = c(paste0("error sd - ", layer_names[q]), rep(NA,6))     
        }

     for(m in 1:((N_responses*2)-1)){
     for(n in (m+1):(N_responses*2)){
        ticker = ticker + 1 
        results_srm_base[ticker,] = sum_stats(paste0("Generalized reciprocity - ", layer_names_long[m], " - ", layer_names_long[n]), samples$srm_model_samples$G_corr[,m,n], HPDI)
        }}

     for(m in 1:((N_responses*2)-1)){
     for(n in (m+1):(N_responses*2)){
        ticker = ticker + 1 
        results_srm_base[ticker,] = sum_stats(paste0("Dyadic reciprocity - ", layer_names_long2[m], " - ", layer_names_long2[n]), samples$srm_model_samples$D_corr[,m,n], HPDI)
        }}


     group_ids_character_df = cbind(rep("Any",input$data$N_id),attr(input$data, "group_ids_character"))

    if(is.null(colnames(group_ids_character_df))){
        colnames(group_ids_character_df) = paste0("(NoBlockingVars)", 1:ncol(group_ids_character_df))
     }

     if(is.null(dim(attr(input$data, "group_ids_character")))){
        colnames(group_ids_character_df) = c("(Intercept)","Scrap")
     }
     
     colnames(group_ids_character_df)[1] = "(Intercept)"
     in_IDs = colnames(input$data$block_set)
     all_IDs = colnames(group_ids_character_df)
     group_ids_character_df = group_ids_character_df[,match(in_IDs, all_IDs), drop = FALSE]

     group_id_levels = append("Any", attr(input$data, "group_ids_levels"), 1)
     names(group_id_levels)[1]= "(Intercept)"
     
     for(l in 1:N_responses){
       B_set = samples$srm_model_samples$block_parameters[[l]]
     for(q in 1:input$data$N_group_vars){
      group_ids_character = levels(factor(group_ids_character_df[,q]))
      test_sorting = group_id_levels[[which(names(group_id_levels) == colnames(group_ids_character_df)[q])]]
      if(all(group_ids_character==test_sorting)==FALSE){
        stop("Factors not sorted correctly.")
      }

      
      for(b1 in 1:input$data$N_groups_per_var[q]){
      for(b2 in 1:input$data$N_groups_per_var[q]){
       ticker = ticker + 1  
      results_srm_base[ticker,] = sum_stats(paste0("offset, ", layer_names[l], " - ", group_ids_character[b1], " to ", group_ids_character[b2]), 
                                                                         B_set[[q]][,b1,b2], HPDI)
     }}

     }}
     
     results_list[[4]] = results_srm_base

   ############# Finally, merge all effects into a list
     for(i in 1:4)
     colnames(results_list[[i]]) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 

     names(results_list) = c( "Focal effects: Out-degree", "Target effects: In-degree", "Dyadic effects", "Other estimates")
          
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
