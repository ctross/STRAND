#' An internal function to check that a list of STRAND data objects can be merged into a longitudinal object
#' 
#' @param data A list of data objects of class STRAND prepared using the make_strand_data() function. The data object must include all covariates used in the formulas of the fit_longitudinal_model() function.
#' @param pass 1 or 2. First pass checks raw data, second pass checks some computed values also.
#' @return Returns an error if objects do not match.
#' @export
#' 

run_longitudinal_checks = function(data, pass){
   if(pass==1){
    check_set = matrix(NA, nrow=length(data),ncol=12)
    checks = c("N_networktypes", "N_id", "N_responses", "N_periods", "N_individual_predictors", "N_dyadic_predictors", 
             "N_block_predictors", "individual_predictors", "block_predictors", "dyadic_predictors",  
             "outcome_mode", "N_groups_per_block_type")
     }
   if(pass==2){ 
    check_set = matrix(NA, nrow=length(data),ncol=18)
    checks = c("N_networktypes", "N_id", "N_responses", "N_periods", "N_individual_predictors", "N_dyadic_predictors", 
             "N_block_predictors", "individual_predictors", "block_predictors", "dyadic_predictors",  
             "outcome_mode", "N_groups_per_block_type", "N_params", "N_group_vars","N_groups_per_var", 
             "max_N_groups", "group_ids_levels","priors" )
     }

  for(i in 1:length(data)){
   if(!"Longitudinal" %in% attr(data[[i]], "supported_models")){stop("For longitudinal models, set 'longitudinal = TRUE' in call to make_strand_data().")}

   check_set[i,1] = ifelse(data[[1]]$N_networktypes == data[[i]]$N_networktypes,1,0)
   check_set[i,2] = ifelse(data[[1]]$N_id == data[[i]]$N_id,1,0)
   check_set[i,3] = ifelse(data[[1]]$N_responses == data[[i]]$N_responses,1,0)
   check_set[i,4] = ifelse(data[[1]]$N_periods == data[[i]]$N_periods,1,0)
   check_set[i,5] = ifelse(data[[1]]$N_individual_predictors == data[[i]]$N_individual_predictors,1,0)
   check_set[i,6] = ifelse(data[[1]]$N_dyadic_predictors == data[[i]]$N_dyadic_predictors,1,0)
   check_set[i,7] = ifelse(data[[1]]$N_block_predictors == data[[i]]$N_block_predictors,1,0)
   check_set[i,8] = ifelse(all(colnames(data[[1]]$individual_predictors)==colnames(data[[i]]$individual_predictors)),1,0)
   check_set[i,9] = ifelse(all(colnames(data[[1]]$block_predictors)==colnames(data[[i]]$block_predictors)),1,0)
   check_set[i,10] = ifelse(all(names(data[[1]]$dyadic_predictors)==names(data[[i]]$dyadic_predictors)),1,0)
   check_set[i,11] = ifelse(data[[1]]$outcome_mode == data[[i]]$outcome_mode,1,0)
   check_set[i,12] = ifelse(all(data[[1]]$N_groups_per_block_type == data[[i]]$N_groups_per_block_type),1,0)
   }
   
   if(pass==2){ 
    for(i in 1:length(data)){ 
   check_set[i,13] = ifelse(all(data[[1]]$N_params == data[[i]]$N_params),1,0)
   check_set[i,14] = ifelse(all(data[[1]]$N_group_vars == data[[i]]$N_group_vars),1,0)
   check_set[i,15] = ifelse(all(data[[1]]$N_groups_per_var == data[[i]]$N_groups_per_var),1,0)
   check_set[i,16] = ifelse(all(data[[1]]$max_N_groups == data[[i]]$max_N_groups),1,0)
  
    temp1 = c()
    for(j in 1:length(attr(data[[1]], "group_ids_levels"))){ 
     temp1[j] = all(attr(data[[1]], "group_ids_levels")[[j]]==attr(data[[i]], "group_ids_levels")[[j]])
     }
    check_set[i,17] = all(temp1)
    check_set[i,18] = ifelse(all(data[[1]]$priors == data[[i]]$priors),1,0)
     }
    }

  cschecks = ifelse(colSums(check_set) == length(data),1,0)

  for(i in 1:length(cschecks)){
    if(cschecks[i]==0){
        stop(paste0("Error, match not possible for ", checks[i]))
    }
  }
}
