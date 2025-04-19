#' A function to create priors for STRAND models
#' 
#' This function allows users to change the numerical values of priors (but not their distributions). Users should only change these values if they have read and understood the source code. 
#'
#' @param priors_to_change A named vector of priors to change. To get parameter names, run: make_priors(include_rownames=TRUE). Values must be 2-vectors. If the prior has only one parameter, inlcude 0 as the second element.
#' @param include_rownames Should rownames be printed?
#' @return A STRAND priors matrix.
#' @export
#' @examples
#' \dontrun{
#' make_priors(priors_to_change=list(sr_L=5),include_rownames=TRUE)
#' make_priors()
#' }
#' 

make_priors = function(priors_to_change=NULL, include_rownames=FALSE){

 parameter_names = c("false_positive_rate", "recall_of_true_ties", "theta_mean",  
                     "fpr_sigma", "rtt_sigma", "theta_sigma",
                     "fpr_effects", "rtt_effects", "theta_effects",
                     "B_ingroup", "B_outgroup",
                     "focal_effects", "target_effects", "dyad_effects",
                     "sr_sigma", "dr_sigma", 
                     "sr_L", "dr_L",
                     "penalty",
                     "effect_max",
                     "effect_decay",
                     "flow_rate",
                     "gaussian_error_priors"
                     )
   
 priors = matrix(0, nrow=23, ncol=2)
  priors[1,1] = -3
  priors[1,2] = 1.5

  priors[2,1] = 3
  priors[2,2] = 1.5

  priors[3,1] = -1.5
  priors[3,2] = 1

  priors[4,1] = 1
  priors[5,1] = 1
  priors[6,1] = 1

  priors[7,1] = 0
  priors[7,2] = 2.5

  priors[8,1] = 0 
  priors[8,2] = 2.5

  priors[9,1] = 0
  priors[9,2] = 2.5

  priors[10,1] = 0.1
  priors[10,2] = 2.5

  priors[11,1] = 0.01
  priors[11,2] = 2.5

  priors[12,1] = 0
  priors[12,2] = 2.5

  priors[13,1] = 0 
  priors[13,2] = 2.5 

  priors[14,1] = 0
  priors[14,2] = 2.5

  priors[15,1] = 0
  priors[15,2] = 2.5

  priors[16,1] = 0
  priors[16,2] = 2.5

  priors[17,1] = 2.5
  priors[18,1] = 2.5

  priors[19,1] = 1.5

  priors[20,1] = 3
  priors[20,2] = 1

  priors[21,1] = 2

  priors[22,1] = 3
  priors[22,2] = 12

  priors[23,1] = 0
  priors[23,2] = 2.5

  if(!is.null(priors_to_change)){
  for(i in 1:length(priors_to_change)){
    if(length(priors_to_change[[i]])==1){
   priors[which(parameter_names==names(priors_to_change)[i] ),1] = priors_to_change[[i]]
    }
   if(length(priors_to_change[[i]])==2){
   priors[which(parameter_names==names(priors_to_change)[i] ),1:2] = priors_to_change[[i]]
    }
  }}

  if(include_rownames == TRUE){
   rownames(priors) = parameter_names
  }

  return(priors)
}
