#' A function to create a STRAND data object
#'
#' This function organizes network data and covariates into a form that can be used by STRAND for model fitting. All 
#' STRAND model fitting functions require their data to be supplied in the format exported here.
#'
#' @param 
#' self_report A list of primary network data (e.g., self reports). Each entry in the list must be an adjacency matrix. This will be a list of length 1 for single-sampled networks
#' and a list of length 2 for double-sampled networks. Data is presumed to be organized such that self_report[[1]][i,j] represents i's reports of transfers from i to j, and self_report[[2]][i,j]
#' represents i's reports of transfers from j to i. Data should be binary, 0 or 1, unless an alternative outcome_mode is provided. If outcome_mode="poisson", then data can be integer values.
#' If an exposure variable is provided, self_report can take integer values and outcome_mode="binomial" can be set.
#' @param
#' outcome_mode Can be either "bernoulli", "binomial", or "poisson", based on the kind of network data being modeled.
#' @param 
#' ground_truth A list of secondary network data about equivalent latent relationships (i.e., from focal observations). Each entry in the list must be an adjacency matrix. 
#' Data is presumed to be organized such that ground_truth[[t]][i,j] represents observed transfers from i to j at time-point t.
#' @param 
#' block_covariates A vector of group IDs (e.g., ethnicity, class, religion, etc.) corresponding to the individuals in the 'self_report' network(s). This should be provided as a factor.
#' @param 
#' individual_covariates An N_id by N_parameters dataframe of all individual-level covariates that are to be included in the model.
#' @param 
#' dyadic_covariates A list of N_id by N_id by N_dyadic_parameters matrices.
#' @param 
#' exposure A list of matrices matched to the self_report matrices. If self_report is a count data set with binomial outcomes, then this variable holds the sample size information.
#' @return A list of data formatted for use by STRAND models.
#' @export
#' @examples
#' \dontrun{
#' model_dat = make_strand_data(self_report=LoanData)
#' }
#'

make_strand_data = function(outcome=NULL, self_report=NULL, outcome_mode="bernoulli", ground_truth=NULL, block_covariates=NULL, individual_covariates=NULL, dyadic_covariates=NULL, exposure=NULL){

         ############################################################################# Check inputs
         ###################### Outcome mode
         outcome_mode_numeric = NULL

         if(outcome_mode=="bernoulli"){
          outcome_mode_numeric = 1
         }

         if(outcome_mode=="binomial"){
          outcome_mode_numeric = 2
         }

         if(outcome_mode=="poisson"){
          outcome_mode_numeric = 3
         }

         if(is.null(outcome_mode_numeric)) stop("outcome_mode not supported")

         # Renames self-report if needed
         if(is.null(self_report)){
          self_report = outcome
         }

         # Check self_report data
         if(is.null(exposure)){
         if(is.null(self_report)) stop("self_report must be a list of matrices.")
         if(!is.list(self_report)) stop("self_report must be a list of matrices.")
         if(!length(self_report) %in% c(1,2)) stop("self_report must be a list length 1 or 2.")
         for(i in 1:length(self_report)){
          if(outcome_mode=="bernoulli"){
          if(!all(self_report[[i]] %in% c(0,1))) stop("self_report must be binary 0 or 1")
          }
         }
         }

         if(!is.null(exposure)){
           if(length(self_report) != length(exposure)) stop("self_report and exposure must be lists of matrices equal in length.")
         }

         # Check ground_truth data
         if(!is.null(ground_truth)){ 
         if(!is.list(ground_truth)) stop("ground_truth must be a list of matrices.")
         }

         # Check block_covariates data
         if(!is.null(block_covariates)){ 
         if(!is.data.frame(block_covariates)) stop("block_covariates must be a data frame.")
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

         if(is.null(exposure)){
          exposure_on = 0
          exposure_risk = array(0, c(N_id, N_id, N_responses))
         } else{
          exposure_on = 1
          exposure_risk = array(NA, c(N_id, N_id, N_responses))
          for(i in 1:length(exposure)){
          exposure_risk[,,i] = exposure[[i]]
          }
         }
         
         for(i in 1:length(self_report)){
          outcomes[,,i] = self_report[[i]]
         }

        if(is.null(block_covariates)){
         block_covariates = rep(1, N_id)
         N_groups_per_type = 1
         N_block_types = 0
         group_ids_character = rep("Any", N_id)
         group_ids = block_covariates
         group_ids_levels = "No Blocks"
          } else{
          
         N_block_types = length(block_covariates[1,]) 
         N_groups_per_type = rep(NA, N_block_types)
         group_ids_character = array(NA, c(N_id, N_block_types))
         group_ids = array(NA, c(N_id, N_block_types))
         group_ids_levels = vector("list", N_block_types)

         for(i in 1:N_block_types){
          N_groups_per_type[i] = max(as.numeric(block_covariates[,i]))
          group_ids_character[,i] = as.character(block_covariates[,i]) 
          group_ids[,i] = as.numeric(block_covariates[,i])
          group_ids_levels[[i]] = levels(block_covariates[,i])
         }
          group_ids = data.frame(group_ids)
          colnames(group_ids) = colnames(block_covariates)
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
          individual_predictors = 0
         } else{
          N_individual_predictors = dim(individual_covariates)[2]  
          individual_predictors = individual_covariates
         }

        if(is.null(dyadic_covariates)){
          N_dyadic_predictors = 0
          dyadic_predictors = 0
         } else{
          N_dyadic_predictors = length(dyadic_covariates)
          dyadic_predictors = dyadic_covariates
         }

         ############################################################################# Determine legal models
         if(N_responses==1){
          if(max(N_groups_per_type)>1){
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
     N_responses = N_responses,  
     N_periods=N_periods,                                              
     N_individual_predictors = N_individual_predictors,                                       
     N_dyadic_predictors =  N_dyadic_predictors,                                           
     outcomes = outcomes,  
     flows=flows,                         
     individual_predictors = individual_predictors,      
     dyadic_predictors = dyadic_predictors,
     N_block_predictors = N_block_types,
     N_groups_per_block_type = N_groups_per_type,
     block_predictors = group_ids,
     outcome_mode=outcome_mode_numeric,
     exposure=exposure_risk
     )

   attr(model_dat, "class") = "STRAND Data Object"
   attr(model_dat, "supported_models") = supported_models
   attr(model_dat, "group_ids_character") = group_ids_character
   attr(model_dat, "group_ids_levels") = group_ids_levels
   colnames(attr(model_dat, "group_ids_character"))=colnames(model_dat$block_predictors)
   
  return(model_dat)
}
