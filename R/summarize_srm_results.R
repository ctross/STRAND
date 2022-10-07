#' Organize Stan output and provide summaries of model parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters
#'
#' @param 
#' input A STRAND model object, obtained by fitting a social relations model.
#' @param 
#' include_samples An indicator for the user to specify where raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param 
#' HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A STRAND results object including summary table, a summary list, and samples.
#' @export
#' @examples
#' \dontrun{
#' res = summarize_srm_results(input=fit)
#' }
#'

# Should change to allow users to specify HPDI intervals

summarize_srm_results = function(input, include_samples=TRUE, HPDI=0.9){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_srm_results() requires a fitted object of class: STRAND Model Object. Please use fit_social_relations_model() to run your model.")
    }

    if(attributes(input)$fit_type != "mcmc"){
        stop("Fitted results can only be reorganized for STRAND model objects fit using MCMC. Variational inference or optimization can be used in Stan
              during experimental model runs, but final inferences should be based on MCMC sampling.")   
    }

    ###################################################### Create samples 
    fit = input$fit
    stanfit = rstan::read_stan_csv(fit$output_files())

    ################### Network model
    B = rstan::extract(stanfit, pars="B")$B  

    sr_sigma = rstan::extract(stanfit, pars="sr_sigma")$sr_sigma  
    sr_L = rstan::extract(stanfit, pars="sr_L")$sr_L  
    sr_raw = rstan::extract(stanfit, pars="sr_raw")$sr_raw  

    dr_sigma = rstan::extract(stanfit, pars="dr_sigma")$dr_sigma  
    dr_L = rstan::extract(stanfit, pars="dr_L")$dr_L  
    dr_raw = rstan::extract(stanfit, pars="dr_raw")$dr_raw  
    
    if(dim(input$data$focal_set)[2]>1)
    focal_effects = rstan::extract(stanfit, pars="focal_effects")$focal_effects  
    if(dim(input$data$target_set)[2]>1)
    target_effects = rstan::extract(stanfit, pars="target_effects")$target_effects  
    if(dim(input$data$dyad_set)[3]>1)
    dyad_effects = rstan::extract(stanfit, pars="dyad_effects")$dyad_effects  

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

    samples = list(srm_model_samples=srm_samples)

   if(input$return_predicted_network == TRUE){
        samples$predicted_network_sample = rstan::extract(stanfit, pars="p")$p  
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

     results_srm_focal = matrix(NA, nrow=(1+Q1) , ncol=6)
     results_srm_target = matrix(NA, nrow=(1+Q2) , ncol=6)
     results_srm_dyadic = matrix(NA, nrow=(1+Q3) , ncol=6)

     results_srm_focal[1,] = sum_stats("focal effects sd", samples$srm_model_samples$focal_target_sd[,1], HPDI)
     if(Q1>0){
     coeff_names = colnames(input$data$focal_set)[-1]
        for(i in 1:Q1){
     results_srm_focal[1+i,] = sum_stats(paste0("focal effects coeffs (out-degree), ", coeff_names[i] ), samples$srm_model_samples$focal_coeffs[,i], HPDI)
        }
      }

      results_list[[1]] = results_srm_focal

     results_srm_target[1,] = sum_stats("target effects sd", samples$srm_model_samples$focal_target_sd[,2], HPDI)
     if(Q2>0){
     coeff_names = colnames(input$data$target_set)[-1]
        for(i in 1:Q2){
     results_srm_target[1+i,] = sum_stats(paste0("target effects coeffs (in-degree), ", coeff_names[i] ), samples$srm_model_samples$target_coeffs[,i], HPDI)
        }
      }

      results_list[[2]] = results_srm_target

     results_srm_dyadic[1,] = sum_stats("dyadic effects sd", c(samples$srm_model_samples$dyadic_sd), HPDI)
     if(Q3>0){
     coeff_names = dimnames(input$data$dyad_set)[[3]][-1]
        for(i in 1:Q3){
     results_srm_dyadic[1+i,] = sum_stats(paste0("dyadic effects coeffs, ", coeff_names[i] ), samples$srm_model_samples$dyadic_coeffs[,i], HPDI)
        }
      }
     results_list[[3]] = results_srm_dyadic

     results_srm_base = matrix(NA, nrow=3, ncol=6)
     results_srm_base[1,] = sum_stats("focal-target effects rho (generalized recipocity)", samples$srm_model_samples$focal_target_L[,2,1], HPDI)
     results_srm_base[2,] = sum_stats("dyadic effects rho (dyadic recipocity)", samples$srm_model_samples$dyadic_L[,2,1], HPDI)
     results_srm_base[3,] = sum_stats("intercept, any to any", samples$srm_model_samples$block_parameters[,1,1], HPDI)
     
     results_list[[4]] = results_srm_base

     for(i in 1:4)
     colnames(results_list[[i]]) = c("Variable", "Median", "HPDI:0.05","HPDI:0.95","Mean","SD") 

     names(results_list) = c( "Focal efffects: Out-degree", "Target effects: In-degree", "Dyadic effects", "Other estimates")
          
   results_out = rbind( results_srm_focal, results_srm_target,results_srm_dyadic, results_srm_base)
   
   df = data.frame(results_out)
   colnames(df) = c("Variable", "Median", "HPDI:0.05","HPDI:0.95","Mean","SD") 

   res_final = list(summary=df, summary_list=results_list)

  if(include_samples==TRUE){
    res_final$samples=samples
   }

   print(results_list)

    attr(res_final, "class") = "STRAND Results Object"
    return(res_final)
}



