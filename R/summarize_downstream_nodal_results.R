#' Organize Stan output and provide summaries of model parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters.
#'
#' @param input A STRAND model object, obtained by fitting a downstream nodal model.
#' @param include_samples An indicator for the user to specify where raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A STRAND results object including summary table, a summary list, and samples.
#' @export
#' @examples
#' \dontrun{
#' res = summarize_downstream_nodal_results(input = fit)
#' }
#'

summarize_downstream_nodal_results = function(input, include_samples=TRUE, HPDI=0.9){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_downstream_nodal_results() requires a fitted object of class: STRAND Model Object. Please use fit_downstream_nodal_model() to run your model.")
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

    ################### Model parameters
    sr_raw = posterior::draws_of(stanfit$"sr_raw")
    error_sigma = posterior::draws_of(stanfit$"error_sigma")

    if(dim(input$data$focal_set)[2]>1)
    focal_effects = posterior::draws_of(stanfit$"focal_effects") 

    alpha = posterior::draws_of(stanfit$"alpha")
    kappa = posterior::draws_of(stanfit$"kappa") 

    srm_samples = list(
            sender_receiver_effects=sr_raw,
            error_sd = error_sigma,
            intercept=alpha,
            sender_slope=kappa[,1],
            receiver_slope=kappa[,2]            
        )

    if(dim(input$data$focal_set)[2]>1)
    srm_samples$focal_coeffs = focal_effects

    samples = list(srm_model_samples=srm_samples)

    ###################################################### Create summary stats 
     results_list = list()

    ################### Regression model
     Q1 = dim(input$data$focal_set)[2]-1

     results_srm_focal = matrix(NA, nrow=(4+Q1) , ncol=7)

    ######### Calculate all focal effects
     results_srm_focal[1,] = sum_stats("intercept", samples$srm_model_samples$intercept, HPDI)

     if(input$data$outcome_mode %in% c(4,5,6,7)){
     results_srm_focal[2,] = sum_stats("sigma", samples$srm_model_samples$error_sd, HPDI)
      } else{
        results_srm_focal[2,1] = "sigma"
      }

     if(input$data$Z[1]==1){
     results_srm_focal[3,] = sum_stats("out-strength", samples$srm_model_samples$sender_slope, HPDI)
      } else{
        results_srm_focal[3,1] = "out-strength"
      }

     if(input$data$Z[2]==1){
     results_srm_focal[4,] = sum_stats("in-strength", samples$srm_model_samples$receiver_slope, HPDI)
      } else{
        results_srm_focal[4,1] = "in-strength"
      }

     if(Q1>0){
     coeff_names = colnames(input$data$focal_set)[-1]
        for(i in 1:Q1){
     results_srm_focal[4+i,] = sum_stats(coeff_names[i], samples$srm_model_samples$focal_coeffs[,i], HPDI)
        }
      }

      results_list[[1]] = results_srm_focal


   ############# Finally, merge all effects into a list
     for(i in 1:1)
     colnames(results_list[[i]]) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 

     names(results_list) = c("Downstream nodal model")
          
   results_out = rbind(results_srm_focal)
   
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
