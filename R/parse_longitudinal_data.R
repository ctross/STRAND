#' An internal function to parse a list of STRAND data objects into a longitudinal object and then apply regression models.
#' 
#' @param data A data objects of class STRAND, prepared using the make_strand_data() function. The data object must include all covariates used in the formulas of the fit_longitudinal_model function.
#' @param block_regression A formula for the block-level predictors. This should be specified as in lm(), e.g.: ~ Ethnicity + Sex. Dont use interactions, however.
#' @param focal_regression A formula for the predictors of out-degree (i.e., focal effects, or the effects of individual covariates on outgoing ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param target_regression A formula for the predictors of in-degree (i.e., target effects, or the effects of individual covariates on incoming ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param dyad_regression A formula for the predictors of dyadic relationships. This should be specified as in lm(), e.g.: ~ Kinship + Friendship
#' @return A STRAND data object.
#' @export
#' 

parse_longitudinal_data = function(data,
                                   block_regression,
                                   focal_regression,
                                   target_regression,
                                   dyad_regression
                                    ){
    ############################################################################# Prepare data and parse formulas
     ind_names = colnames(data$individual_predictors)
     dyad_names = names(data$dyadic_predictors)

    if(data$N_individual_predictors==0 & focal_regression != ~ 1){
        stop("No individual covariate data has been provided. focal_regression must equal ~ 1 ")
    }

    if(data$N_individual_predictors==0 & target_regression != ~ 1){
        stop("No individual covariate data has been provided. target_regression must equal ~ 1 ")
    }

    if(data$N_dyadic_predictors==0 & dyad_regression != ~ 1){
        stop("No individual covariate data has been provided. dyad_regression must equal ~ 1 ")
    }

    if(data$N_block_predictors==0 & block_regression != ~ 1){
        stop("No block covariate data has been provided. block_regression must equal ~ 1 ")
    }

    ################################################################ Dyad model matrix
     if(data$N_dyadic_predictors>0){
     dyad_dims = c(data$N_id, data$N_id, length(dyad_names))

     dyad_dat = list()
     for(i in 1:dyad_dims[3]){
      dyad_dat[[i]] = c(data$dyadic_predictors[[i]])  
     }

     #dyad_dat = do.call(rbind.data.frame, dyad_dat)
     #dyad_dat = as.data.frame(do.call(cbind, dyad_dat))
     dyad_dat = do.call(data.frame, dyad_dat)
     colnames(dyad_dat) = dyad_names
     dyad_model_matrix = model.matrix( dyad_regression , dyad_dat )

     dyad_dat_out = array(NA, c(dyad_dims[1], dyad_dims[2], ncol(dyad_model_matrix)))
     for(i in 1:ncol(dyad_model_matrix)){
      dyad_dat_out[,,i] = matrix(dyad_model_matrix[,i], nrow=dyad_dims[1], ncol=dyad_dims[2])  
     }

     dimnames(dyad_dat_out)[[3]] = colnames(dyad_model_matrix)
     data$dyad_set = dyad_dat_out
     } else{
      data$dyad_set = array(1, c(data$N_id, data$N_id, 1))
     }

     ################################################################ Individual model matrix
     if(data$N_individual_predictors>0){
      data$focal_set = model.matrix( focal_regression , data$individual_predictors )
      data$target_set = model.matrix( target_regression , data$individual_predictors )
     } else{
      data$focal_set = matrix(1,nrow=data$N_id, ncol=1)
      data$target_set = matrix(1,nrow=data$N_id, ncol=1)
     }
    
    data$N_params = c(ncol(data$focal_set), ncol(data$target_set), dim(data$dyad_set)[3])

    ################################################################ Block model matrix
     if(data$N_block_predictors>0){
      data$block_set = model.matrix( block_regression , data$block_predictors )
     } else{
      data$block_set = as.array(matrix(1, nrow=data$N_id, ncol=1))
      colnames(data$block_set) = "(Intercept)"
     }

     data$N_group_vars = ncol(data$block_set) 
     data$N_groups_per_var = rep(NA, data$N_group_vars)

     for(i in 1:data$N_group_vars){
      data$N_groups_per_var[i] = length(unique(data$block_set[,i]))  
     }

     data$N_groups_per_var = as.array(data$N_groups_per_var)
  
     data$max_N_groups = max(data$N_groups_per_var)

    return(data)
}
