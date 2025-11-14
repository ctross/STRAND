#' A function to compute relevant data on missings for Stan imputation
#'
#' Simpler helper function.
#'
#' @param data A data object from make_strand_data after its proccessed by a model function call.
#' @return process_missings(data)
#' @export

process_missings = function(data){
      ######### Find missings in focal, target, and dyadic vars
      if(!is.null(data$focal_set)){
        locations_missing_focal_set = which(is.na(data$focal_set),arr.ind=TRUE)
        N_missing_focal_set = nrow(locations_missing_focal_set)
        focal_lims = matrix(NA, nrow=2, ncol=dim(data$focal_set)[2])
        focal_lims[1,] = colMins(data$focal_set, margin=2)
        focal_lims[2,] = colMaxs(data$focal_set, margin=2)
        data$focal_set[is.na(data$focal_set)] = -9999999
        data$N_missing_focal_set = N_missing_focal_set  
        data$locations_missing_focal_set = locations_missing_focal_set  
        data$focal_lims = focal_lims[, -1, drop=FALSE]

            if(N_missing_focal_set>0){
            warning("Missing data detected in the focal predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
            }
      }

      if(!is.null(data$target_set)){
        locations_missing_target_set = which(is.na(data$target_set),arr.ind=TRUE)
        N_missing_target_set = nrow(locations_missing_target_set)
        target_lims = matrix(NA, nrow=2, ncol=dim(data$target_set)[2])
        target_lims[1,] = colMins(data$target_set, margin=2)
        target_lims[2,] = colMaxs(data$target_set, margin=2)
        data$target_set[is.na(data$target_set)] = -9999999
        data$N_missing_target_set = N_missing_target_set  
        data$locations_missing_target_set = locations_missing_target_set  
        data$target_lims = target_lims[, -1, drop=FALSE]

            if(N_missing_target_set>0){
            warning("Missing data detected in the target predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
            }
        }

      if(!is.null(data$dyad_set)){
        locations_missing_dyad_set = which(is.na(data$dyad_set),arr.ind=TRUE)
        N_missing_dyad_set = nrow(locations_missing_dyad_set)
        dyad_lims = matrix(NA, nrow=2, ncol=dim(data$dyad_set)[3])
        dyad_lims[1,] = colMins(data$dyad_set, margin=3)
        dyad_lims[2,] = colMaxs(data$dyad_set, margin=3) 
        data$dyad_set[is.na(data$dyad_set)] = -9999999 
        data$N_missing_dyad_set = N_missing_dyad_set  
        data$locations_missing_dyad_set = locations_missing_dyad_set  
        data$dyad_lims = dyad_lims[, -1, drop=FALSE]

            if(N_missing_dyad_set>0){
            warnings("Missing data detected in the dyadic predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
            }
        }

      if(!is.null(data$sampling_set)){
        locations_missing_sampling_set = which(is.na(data$sampling_set),arr.ind=TRUE)
        N_missing_sampling_set = nrow(locations_missing_sampling_set)
        sampling_lims = matrix(NA, nrow=2, ncol=dim(data$sampling_set)[2])
        sampling_lims[1,] = colMins(data$sampling_set, margin=2)
        sampling_lims[2,] = colMaxs(data$sampling_set, margin=2) 
        data$sampling_set[is.na(data$sampling_set)] = -9999999
        data$N_missing_sampling_set = N_missing_sampling_set
        data$locations_missing_sampling_set = locations_missing_sampling_set 
        data$sampling_lims = sampling_lims[, -1, drop=FALSE]

            if(N_missing_sampling_set>0){
            warning("Missing data detected in the sampling predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
            }
        }            

      if(!is.null(data$censoring_set)){
        locations_missing_censoring_set = which(is.na(data$censoring_set),arr.ind=TRUE)
        N_missing_censoring_set = nrow(locations_missing_censoring_set)
        censoring_lims = matrix(NA, nrow=2, ncol=dim(data$censoring_set)[2])
        censoring_lims[1,] = colMins(data$censoring_set, margin=2)
        censoring_lims[2,] = colMaxs(data$censoring_set, margin=2)
        data$censoring_set[is.na(data$censoring_set)] = -9999999 
        data$N_missing_censoring_set = N_missing_censoring_set
        data$locations_missing_censoring_set = locations_missing_censoring_set    
        data$censoring_lims = censoring_lims[, -1, drop=FALSE]

            if(N_missing_censoring_set>0){
            warning("Missing data detected in the censoring predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
            }
        }

      if(!is.null(data$fpr_set)){
        locations_missing_fpr_set = which(is.na(data$fpr_set),arr.ind=TRUE)
        N_missing_fpr_set = nrow(locations_missing_fpr_set)
        fpr_lims = matrix(NA, nrow=2, ncol=dim(data$fpr_set)[2])
        fpr_lims[1,] = colMins(data$fpr_set, margin=2)
        fpr_lims[2,] = colMaxs(data$fpr_set, margin=2)
        data$fpr_set[is.na(data$fpr_set)] = -9999999
        data$N_missing_fpr_set = N_missing_fpr_set  
        data$locations_missing_fpr_set = locations_missing_fpr_set  
        data$fpr_lims = fpr_lims[, -1, drop=FALSE]

            if(N_missing_fpr_set>0){
            warning("Missing data detected in the fpr predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
            }
        }

      if(!is.null(data$rtt_set)){
        locations_missing_rtt_set = which(is.na(data$rtt_set),arr.ind=TRUE)
        N_missing_rtt_set = nrow(locations_missing_rtt_set)
        rtt_lims = matrix(NA, nrow=2, ncol=dim(data$rtt_set)[2])
        rtt_lims[1,] = colMins(data$rtt_set, margin=2)
        rtt_lims[2,] = colMaxs(data$rtt_set, margin=2)
        data$rtt_set[is.na(data$rtt_set)] = -9999999
        data$N_missing_rtt_set = N_missing_rtt_set  
        data$locations_missing_rtt_set = locations_missing_rtt_set  
        data$rtt_lims = rtt_lims[, -1, drop=FALSE]

            if(N_missing_rtt_set>0){
            warning("Missing data detected in the rtt predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
            }
        }

      if(!is.null(data$theta_set)){
        locations_missing_theta_set = which(is.na(data$theta_set),arr.ind=TRUE)
        N_missing_theta_set = nrow(locations_missing_theta_set)
        theta_lims = matrix(NA, nrow=2, ncol=dim(data$theta_set)[2])
        theta_lims[1,] = colMins(data$theta_set, margin=2)
        theta_lims[2,] = colMaxs(data$theta_set, margin=2)
        data$theta_set[is.na(data$theta_set)] = -9999999
        data$N_missing_theta_set = N_missing_theta_set  
        data$locations_missing_theta_set = locations_missing_theta_set  
        data$theta_lims = theta_lims[, -1, drop=FALSE]

            if(N_missing_theta_set>0){
            warning("Missing data detected in the theta predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
            }
        }

      ######## Use mask to blank out any missings in outcomes, exposure, or mask
        locations_missing_outcomes = which(is.na(data$outcomes))
        locations_missing_exposure = which(is.na(data$exposure))
        locations_missing_mask = which(is.na(data$mask))

        if(length(locations_missing_outcomes)>0){
            warning("Missing data detected in the outcome layer. These outcomes will be sliced out of the model using a masking layer.") 
            data$outcomes[locations_missing_outcomes] = -99 
            data$outcomes_real[locations_missing_outcomes] = -99   
        }

        if(length(locations_missing_exposure)>0){
            warning("Missing data detected in the exposure layer. The corresponding outcomes will be sliced out of the model using a masking layer.") 
            data$exposure[locations_missing_exposure] = -99 
        }

        if(length(locations_missing_mask)>0){
            warning("Missing data detected in the mask layer. The corresponding outcomes will be sliced out of the model using a masking layer. NAs in the mask become 1s.") 
            data$mask[locations_missing_mask] = -99     
        }

        data$mask[locations_missing_outcomes] = 1    
        data$mask[locations_missing_exposure] = 1  
        data$mask[locations_missing_mask] = 1  

        if("LNM" %in% attributes(data)$"supported_models"){
         mirror_mask = ifelse(data$mask[,,1]==1 | t(data$mask[,,2])==1, 1, 0)   
         data$mask[,,1] = mirror_mask
         data$mask[,,2] = mirror_mask
        }

        ######## Use mask to blank out any missings in outcomes, exposure, or mask for ME models
        locations_missing_sampled = which(is.na(data$sampled))
        locations_missing_sampled_exposure = which(is.na(data$sampled_exposure))
        locations_missing_sampled_mask = which(is.na(data$sampled_mask))

        if(length(locations_missing_sampled)>0){
            warning("Missing data detected in the ME sampled layer. These outcomes will be sliced out of the model using a masking layer.") 
        }

        if(length(locations_missing_sampled_exposure)>0){
            warning("Missing data detected in the ME sampled exposure layer. The corresponding outcomes will be sliced out of the model using a masking layer.") 
        }

        if(length(locations_missing_sampled_mask)>0){
            warning("Missing data detected in the ME sampled mask layer. The corresponding outcomes will be sliced out of the model using a masking layer. NAs in the mask become 1s.") 
        }

        data$sampled[locations_missing_sampled] = -99  
        data$sampled_exposure[locations_missing_sampled_exposure] = -99        
        data$sampled_mask[locations_missing_sampled_mask] = -99     

        data$sampled_mask[locations_missing_sampled] = 1    
        data$sampled_mask[locations_missing_sampled_exposure] = 1  
        data$sampled_mask[locations_missing_sampled_mask] = 1  



        locations_missing_detected = which(is.na(data$detected))
        locations_missing_detected_exposure = which(is.na(data$detected_exposure))
        locations_missing_detected_mask = which(is.na(data$detected_mask))

        if(length(locations_missing_detected)>0){
            warning("Missing data detected in the ME detected layer. These outcomes will be sliced out of the model using a masking layer.") 
        }

        if(length(locations_missing_detected_exposure)>0){
            warning("Missing data detected in the ME detected exposure layer. The corresponding outcomes will be sliced out of the model using a masking layer.") 
        }

        if(length(locations_missing_detected_mask)>0){
            warning("Missing data detected in the ME detected mask layer. The corresponding outcomes will be sliced out of the model using a masking layer. NAs in the mask become 1s.") 
        }

        data$detected[locations_missing_detected] = -99  
        data$detected_exposure[locations_missing_detected_exposure] = -99        
        data$detected_mask[locations_missing_sampled_mask] = -99     

        data$detected_mask[locations_missing_detected] = 1    
        data$detected_mask[locations_missing_detected_exposure] = 1  
        data$detected_mask[locations_missing_detected_mask] = 1  


        ######## For block_set we need to do a single resampling and then print stong warning
        block_miss = sum(is.na(data$block_set))

        if(block_miss>0){
          warning("Missing data detected in the block predictors. We do not currently support Bayesian imputation of discrete missings. Missing data will be not imputed using a Bayesian model.
            Missing data will be imputed a single time by randomly sampling from the non-missing data in the same column. This is NOT a rigorous treatment of missingness. 
            However, if missings are sparse, and randomly distributed, this model should be OK. If many data are missing, or missing is non-random, do not use this hack.") 
            
            for(k in 1:ncol(data$block_set)){
                if(sum(is.na(data$block_set[,k]))>0){
              data$block_set[is.na(data$block_set[,k]),k] = sample(na.omit(data$block_set[,k]), size = sum(is.na(data$block_set[,k])), replace = TRUE)
            }}
        }

        return(data)
    } 


