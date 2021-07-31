#' A function to create a STRAND data object
#'
#' This function will organize network data and covariates into a form that can be used by STRAND for model fitting. All 
#' STRAND model fitting functions require their data to be supplied in the format exported here.
#'
#' @param 
#' self_report A list of primary network data (e.g., self reports). Each entry in the list must be an adjacency matrix. This will be a list of length 1 for single-sampled networks
#' and a list of length 2 for double-sampled networks. Data is presumed to be organized such that self_report[[1]][i,j] represents i's reports of transfers from i to j, and self_report[[2]][i,j]
#' represents i's reports of transfers from j to i.
#' @param 
#' ground_truth A list of secondary network data about equivalent latent relationships (i.e., from focal observations). Each entry in the list must be an adjacency matrix. 
#' Data is presumed to be organized such that ground_truth[[t]][i,j] represents observed transfers from i to j at time-point t.
#' @param 
#' group_ids A vector of group IDs (e.g., ethnicity, class, religion, etc.) corresponding to the individuals in the 'self_report' network(s). This should be provided as a factor.
#' @param 
#' individual_covariates An N_id by N_parameters dataframe of all individual-level covariates that are to be included in the model.
#' @param 
#' dyadic_covariates A list of N_id by N_id by N_dyadic_parameters matrices.
#' @return A list of data formatted for use by STRAND models.
#' @export
#' @examples
#' \dontrun{
#' model_dat = make_strand_data(self_report=LoanData)
#' }
#'

make_strand_data = function(self_report, ground_truth=NULL, group_ids=NULL, individual_covariates=NULL, dyadic_covariates=NULL){

         ############################################################################# Check inputs
         # Check self_report data
         if(is.null(self_report)) stop("self_report must be a list of matrices.")
         if(!is.list(self_report)) stop("self_report must be a list of matrices.")
         if(!length(self_report) %in% c(1,2)) stop("self_report must be a list length 1 or 2.")
         for(i in 1:length(self_report))
         if(!all(self_report[[i]] %in% c(0,1))) stop("self_report must be binary 0 or 1")

         # Check ground_truth data
         if(!is.null(ground_truth)){ 
         if(!is.list(ground_truth)) stop("ground_truth must be a list of matrices.")
         }

         # Check group_ids data
         if(!is.null(group_ids)){ 
         if(!is.factor(group_ids)) stop("group_ids must be a factor.")
         }

         # Check individual_covariates data
         if(!is.null(individual_covariates)){ 
         if(!is.data.frame(individual_covariates)) stop("individual_covariates must be a data frame.")
         }

         # Check self_report data
         if(!is.null(dyadic_covariates)){ 
         if(!is.list(dyadic_covariates)) stop("dyadic_covariates must be a list of matrices.")
         } 
          
          # need to add checks that rownames and colnames of all data types match

         ############################################################################# Process data
         N_id =  dim(self_report[[1]])[1]
         N_responses = length(self_report)

         outcomes = array(NA, c(N_id, N_id, N_responses))
         
         for(i in 1:length(self_report))
         outcomes[,,i] = self_report[[i]]

        if(is.null(group_ids)){
         group_ids = rep(1, N_id)
         N_groups = 1
         group_ids_character = rep("Any", N_id)
          } else{
         group_ids_character = as.character(group_ids) 
         group_ids = as.numeric(group_ids)
         N_groups = max(as.numeric(group_ids)) 
         }
 
        if(is.null(ground_truth)){
         N_networktypes = N_responses
         N_periods = 0
         flows=0
         } else{
         N_networktypes = N_responses + 1
         N_periods = length(ground_truth)

         flows = array(NA, c(N_id, N_id, N_periods))
         
         for(i in 1:length(ground_truth))
         flows[,,i] = ground_truth[[i]]

         }

        if(is.null(individual_covariates)){
          N_individual_predictors = 0
         } else{
          N_individual_predictors = dim(individual_covariates)[2]  
          individual_predictors = individual_covariates
         }

        if(is.null(dyadic_covariates)){
          N_dyadic_predictors = 0
         } else{
          N_dyadic_predictors = length(dyadic_covariates)
          dyadic_predictors = dyadic_covariates
         }

         ############################################################################# Determine legal models
         if(N_responses==1){
          if(N_groups>1){
            supported_models = c("SRM", "SBM", "SRM+SBM")
            } else{
            supported_models = c("SRM")
            }
         } 

         if(N_responses==2 & N_networktypes==2){
          supported_models = c("LNM")
         }

         if(N_responses==2 & N_networktypes==3){
          supported_models = c("LNM","LNM+Flows")
         }


   model_dat = list(
     N_networktypes = N_networktypes,                                               
     N_id = N_id,                                                         
     N_groups = N_groups,                                                  
     N_responses = N_responses,  
     N_periods=N_periods,                                              
     N_individual_predictors = N_individual_predictors,                                       
     N_dyadic_predictors =  N_dyadic_predictors,    
     group_ids = group_ids,                                        
     outcomes = outcomes,  
     flows=flows,                         
     individual_predictors = individual_predictors,      
     dyadic_predictors = dyadic_predictors
     )

   attr(model_dat, "class") = "STRAND Data Object"
   attr(model_dat, "supported_models") = supported_models
   attr(model_dat, "group_ids_character") = unique(group_ids_character)
   
  return(model_dat)
}
