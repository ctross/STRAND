#' A function to compute relevant data on missings for Stan imputation
#'
#' Simpler helper function.
#'
#' @param data A data object from make_strand_data after its proccessed by a model function call.
#' @return process_missings(data)
#' @export

process_missings = function(data){
      ######### Find missings in focal, target, and dyadic vars
        locations_missing_focal_set = which(is.na(data$focal_set),arr.ind=TRUE)
        locations_missing_target_set = which(is.na(data$target_set),arr.ind=TRUE)
        locations_missing_dyad_set = which(is.na(data$dyad_set),arr.ind=TRUE)

        N_missing_focal_set = nrow(locations_missing_focal_set)
        N_missing_target_set = nrow(locations_missing_target_set)
        N_missing_dyad_set = nrow(locations_missing_dyad_set)

        focal_lims = matrix(NA, nrow=2, ncol=dim(data$focal_set)[2])
        target_lims = matrix(NA, nrow=2, ncol=dim(data$target_set)[2])
        dyad_lims = matrix(NA, nrow=2, ncol=dim(data$dyad_set)[3])

        focal_lims[1,] = colMins(data$focal_set, margin=2)
        focal_lims[2,] = colMaxs(data$focal_set, margin=2)

        target_lims[1,] = colMins(data$target_set, margin=2)
        target_lims[2,] = colMaxs(data$target_set, margin=2)

        dyad_lims[1,] = colMins(data$dyad_set, margin=3)
        dyad_lims[2,] = colMaxs(data$dyad_set, margin=3)
        
        if(N_missing_focal_set>0){
            warning("Missing data detected in the focal predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
        }

        if(N_missing_target_set>0){
            warning("Missing data detected in the target predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
        }
        
        if(N_missing_dyad_set>0){
            warnings("Missing data detected in the dyadic predictors. If missings are sparse, and randomly distributed, this model should constitute a rigorous treatment of the missing data.
            Missing data will be imputed using a Bayesian model. Missing data parameters will have flat, continous prior support over the range of observed data (min to max).
            If any of these assumptions are problematic, please consider writing your own bespoke model with a generative model of the observation process.") 
        }

      ######### Code -99 into missing spots now (these will get replaced with parameters in Stan, and will NOT be used as data)
        data$focal_set[is.na(data$focal_set)] = -9999999
        data$target_set[is.na(data$target_set)] = -9999999
        data$dyad_set[is.na(data$dyad_set)] = -9999999 
       
      ######## Add to data stack
        data$N_missing_focal_set = N_missing_focal_set  
        data$N_missing_target_set = N_missing_target_set  
        data$N_missing_dyad_set = N_missing_dyad_set  

        data$locations_missing_focal_set = locations_missing_focal_set  
        data$locations_missing_target_set = locations_missing_target_set  
        data$locations_missing_dyad_set = locations_missing_dyad_set  

        data$focal_lims = focal_lims[,-1]
        data$target_lims = target_lims[,-1]
        data$dyad_lims = dyad_lims[,-1]

      ######## Use mask to blank out any missings in outcomes, exposure, or mask
        locations_missing_outcomes = which(is.na(data$outcomes))
        locations_missing_exposure = which(is.na(data$exposure))
        locations_missing_mask = which(is.na(data$mask))

        if(length(locations_missing_outcomes)>0){
            warning("Missing data detected in the outcome layer. These outcomes will be sliced out of the model using a masking layer.") 
        }

        if(length(locations_missing_exposure)>0){
            warning("Missing data detected in the exposure layer. The corresponding outcomes will be sliced out of the model using a masking layer.") 
        }

        if(length(locations_missing_mask)>0){
            warning("Missing data detected in the mask layer. The corresponding outcomes will be sliced out of the model using a masking layer. NAs in the mask become 1s.") 
        }

        data$outcomes[locations_missing_outcomes] = -99  
        data$exposure[locations_missing_exposure] = -99        
        data$mask[locations_missing_mask] = -99     

        data$mask[locations_missing_outcomes] = 1    
        data$mask[locations_missing_exposure] = 1  
        data$mask[locations_missing_mask] = 1  


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


