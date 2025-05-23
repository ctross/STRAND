#' Organize Stan output and provide summaries of model parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters.
#'
#' @param input A STRAND model object, obtained by fitting a model using a STRAND function.
#' @param include_samples An indicator for the user to specify where raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A STRAND results object including summary table, a summary list, and samples.
#' @export
#' @examples
#' \dontrun{
#' res = summarize_strand_results(input = fit)
#' }
#'

summarize_strand_results = function(input, include_samples=TRUE, HPDI=0.9){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("summarize_strand_results() requires a fitted object of class: STRAND Model Object.")
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

    if(attributes(input)$model_type == "SRM"){
       res = summarize_srm_results(input=input, include_samples=include_samples, HPDI=HPDI)
      }

    if(attributes(input)$model_type == "SBM"){
       res = summarize_bm_results(input=input, include_samples=include_samples, HPDI=HPDI)
      }

    if(attributes(input)$model_type == "SRM+SBM"){
       res = summarize_bsrm_results(input=input, include_samples=include_samples, HPDI=HPDI)
      }

    if(attributes(input)$model_type == "SRM+SBM+ME"){
       res = summarize_bsrm_results_with_measurement_bias(input=input, include_samples=include_samples, HPDI=HPDI)
      }

    if(attributes(input)$model_type == "LNM"){
       res = summarize_lnm_results(input=input, include_samples=include_samples, HPDI=HPDI)
      }

    if(attributes(input)$model_type == "LNM+Flows"){
       res = summarize_lnmf_results(input=input, include_samples=include_samples, HPDI=HPDI)
      }

    if(attributes(input)$model_type == "Multiplex"){
       res = summarize_multiplex_bsrm_results(input=input, include_samples=include_samples, HPDI=HPDI)
      }

    if(attributes(input)$model_type == "Longitudinal"){
       res = summarize_longitudinal_bsrm_results(input=input, include_samples=include_samples, HPDI=HPDI)
      }

    return(res)
}



