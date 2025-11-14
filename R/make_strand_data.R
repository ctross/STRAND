#' A function to create a STRAND data object
#'
#' This function organizes network data and covariates into a form that can be used by STRAND for model fitting. All 
#' STRAND model fitting functions require their data to be supplied in the format exported here.
#'
#' @param outcome A named list of primary network data (e.g., self reports). Each entry in the list must be an adjacency matrix. This will be a list of length 1 for single-sampled networks
#' and a list of length 2 for double-sampled networks. For douple-sampled networks, data is presumed to be organized such that outcome[[1]][i,j] represents i's reports of transfers from i to j, and outcome[[2]][i,j]
#' represents i's reports of transfers from j to i. Data should be binary, 0 or 1, unless an alternative outcome_mode is provided. If outcome_mode="poisson", then data can be integer values.
#' If an exposure variable is provided, outcome can take integer values and outcome_mode="binomial" can be set. For multiplex models, the list can be longer than 2, but all layers must be single sampled, and you must set multiplex = TRUE.
#' @param self_report A named list of primary network data (e.g., self reports). This is a deprecated alias for outcome above.
#' @param outcome_mode Can be either "bernoulli", "binomial", "poisson", or "gaussian" based on the kind of network data being modeled.
#' @param link_mode Can be either "logit", "probit", "log", or "identity"; "log" must be used with "poisson" outcomes; "identity" must be used with "gaussian" outcomes; "logit" is default for Bernoulli and Binomial outcomes; "probit" is basically the same as "logit", but the model fits slower.
#' @param ground_truth A list of secondary network data about equivalent latent relationships (i.e., from focal observations). Each entry in the list must be an adjacency matrix. 
#' Data is presumed to be organized such that ground_truth[[t]][i,j] represents observed transfers from i to j at time-point t.
#' @param block_covariates A data.frame of group IDs (e.g., ethnicity, class, religion, etc.) corresponding to the individuals in the 'self_report' network(s). Each variable should be provided as a factor.
#' @param individual_covariates An N_id by N_parameters dataframe of all individual-level covariates that are to be included in the model.
#' @param dyadic_covariates A named list of N_id by N_id by N_dyadic_parameters matrices.
#' @param exposure A named list of matrices matched to the self_report matrices. If self_report is a count data set with binomial outcomes, then this variable holds the sample size information.
#' @param m_e_data A list of integer vectors: list(sampled=sampled, sampled_exposure=sampled_exposure, sampled_mask=sampled_mask, detected=detected, detected_exposure=detected_exposure, detected_mask=detected_mask), to be used in measurement error models.
#' @param mask A list of matrices matched to the self_report matrices. If mask[i,j,m]==0, then ties between i and j in layer m are detectable. If mask[i,j,m]==1, then i to j ties in layer m are censored (e.g., if i and j were monkeys kept in different enclosures).
#' @param diffusion_outcome An N-vector of outcome data for a trait diffusing over a network.
#' @param diffusion_exposure An N-vector matched with the diffusion_outcome matrix. If diffusion_outcome is a count data set with binomial outcomes, then this variable holds the sample size information.
#' @param diffusion_mask An N-vector of indicators for diffusion outcomes that were masked.
#' @param directed If TRUE, then STRAND will treat the outcomes as directed. If set to FALSE, STRAND will treat the outcomes as undirected; this leads to some addition checks on model definition.
#' @param multiplex If TRUE, then all layers in outcome are modeled jointly.
#' @param imputation If TRUE, then checks for NAs in data are omitted, and supported models will impute the missings.
#' @param longitudinal If TRUE, then checks for longitudinal data structure are performed.
#' @param check_data_organization If TRUE, then checks that all colnames and rownames match. This will catch missorted data.
#' @param check_standardization If TRUE, then checks that all individual and dyadic variables are standardized to have SDs in the range of standardization_threshold. Standardization is important in STRAND so that priors have equal strength across all predictors.
#' @param standardization_threshold If check_standardization is TRUE, then individual and dyadic predictors must have SDs in this range.
#' @return A list of data formatted for use by STRAND models.
#' @export
#' @examples
#' \dontrun{
#' model_dat = make_strand_data(self_report=LoanData)
#' }
#'

make_strand_data = function(outcome=NULL, self_report=NULL, outcome_mode=NULL, link_mode=NULL, ground_truth=NULL, block_covariates=NULL, individual_covariates=NULL, dyadic_covariates=NULL, 
                            exposure=NULL, m_e_data = NULL, mask=NULL, diffusion_outcome = NULL, diffusion_exposure = NULL, diffusion_mask = NULL, multiplex = FALSE, longitudinal = FALSE, 
                            directed = TRUE, imputation = FALSE, check_data_organization = TRUE, check_standardization = TRUE, standardization_threshold = c(0.5, 2)){

         ############################################################################# Check inputs
         # Renames self-report if needed
         if(is.null(self_report)){
          self_report = outcome
         }

         if(longitudinal==TRUE & length(self_report)>1){
           stop("Longitudinal models require single-layer networks.")
         }

         ######################################################## Check types
         if(is.null(outcome_mode)){
          stop("outcome_mode must be either bernoulli, binomial, poisson, or gaussian.")
         }

         if(is.null(link_mode)){
          stop("link_mode must be either logit, probit, log, or identity.")
         }

         if(!is.null(outcome) & !is.list(outcome)){
          stop("outcome must be a list of adjacency matrices, even if length is one.")
         }

         if(!is.null(self_report) & !is.list(self_report)){
          stop("self_report must be a list of adjacency matrices, even if length is one.")
         }

         if(!is.null(outcome_mode) & (!outcome_mode %in% c("bernoulli", "binomial", "poisson", "gaussian"))){
          stop("outcome_mode must be either bernoulli, binomial, poisson, or gaussian.")
         }

         if(!is.null(link_mode) & (!link_mode %in% c("logit", "probit", "log", "identity"))){
          stop("link_mode must be either logit, probit, log, or identity.")
         }

         if(!is.null(ground_truth) & !is.list(ground_truth)){
          stop("ground_truth must be a list of adjacency matrices, even if length is one.")
         }

         if(!is.null(block_covariates) & !is.data.frame(block_covariates)){
          stop("block_covariates must be a data.frame.")
         }

         if(!is.null(individual_covariates) & !is.data.frame(individual_covariates)){
          stop("individual_covariates must be a data.frame.")
         }

         if(!is.null(dyadic_covariates) & !is.list(dyadic_covariates)){
          stop("dyadic_covariates must be a list of adjacency matrices, even if length is one.")
         }

        if(!is.null(exposure) & !is.list(exposure)){
          stop("exposure must be a list of adjacency matrices, even if length is one.")
         }

        if(!is.null(m_e_data) & !is.list(m_e_data)){
          stop("m_e_data must be a list.")
         }

        if(!is.null(mask) & !is.list(mask)){
          stop("mask must be a list of adjacency matrices, even if length is one.")
         }

        if(!is.null(diffusion_outcome) & !is.vector(diffusion_outcome)){
          stop("diffusion_outcome must be a vector.")
         }

        if(!is.null(diffusion_exposure) & !is.vector(diffusion_exposure)){
          stop("diffusion_exposure must be a vector.")
         }

        if(!is.null(diffusion_mask) & !is.vector(diffusion_mask)){
          stop("diffusion_mask must be a vector.")
         }

         if(!multiplex %in% c(TRUE, FALSE)){
          stop("multiplex must be either TRUE or FALSE.")
         }

         if(!longitudinal %in% c(TRUE, FALSE)){
          stop("longitudinal must be either TRUE or FALSE.")
         }

         if(!check_data_organization %in% c(TRUE, FALSE)){
          stop("check_data_organization must be either TRUE or FALSE.")
         }
      
      ################################## Check names
         if(!is.null(outcome) & is.list(outcome) & is.null(names(outcome))){
          stop("outcome must be a *named* list of adjacency matrices, even if length is one. Please provide a name/label for each of the outcome matrices.  ")
         }

         if(!is.null(self_report) & is.list(self_report) & is.null(names(self_report))){
          stop("self_report must be a *named* list of adjacency matrices, even if length is one. Please provide a name/label for each of the outcome matrices. ")
         }

         if(!is.null(ground_truth) & is.list(ground_truth) & is.null(names(ground_truth))){
          stop("ground_truth must be a *named* list of adjacency matrices, even if length is one. Please provide a name/label for each of the ground_truth matrices. ")
         }

         if(!is.null(dyadic_covariates) & is.list(dyadic_covariates) & is.null(names(dyadic_covariates))){
          stop("dyadic_covariates must be a *named* list of adjacency matrices, even if length is one.Please provide a name/label for each of the dyadic_covariates matrices. ")
         }

         if(!is.null(exposure) & is.list(exposure) & is.null(names(exposure))){
          stop("exposure must be a *named* list of adjacency matrices, even if length is one. Please provide a name/label (matched to those of 'outcome') for each of the exposure matrices.")
         }

         if(!is.null(mask) & is.list(mask) & is.null(names(mask))){
          stop("mask must be a *named* list of adjacency matrices, even if length is one. Please provide a name/label (matched to those of 'outcome') for each of the mask matrices.")
         }


      ###################### Outcome mode
         if(is.null(outcome_mode)) stop(" 'outcome_mode' must be declared.")
         if(is.null(link_mode)) stop(" 'link_mode' must be declared.")

         outcome_mode_numeric = NULL
         link_mode_numeric = NULL

         if(outcome_mode=="bernoulli"){
          outcome_mode_numeric = 1

          if(!link_mode %in% c("logit", "probit")){stop("If outcome_mode is 'bernoulli', you must set link_mode to 'logit' or 'probit'.")}

          if(link_mode == "logit"){
            link_mode_numeric = 1
          } 

          if(link_mode == "probit"){
            link_mode_numeric = 2
          } 

         }

         if(outcome_mode=="binomial"){
          if(is.null(exposure)){stop("If outcome is binomial, an exposure variable must be provided.")}
          outcome_mode_numeric = 2

          if(!link_mode %in% c("logit", "probit")){stop("If outcome_mode is 'binomial', you must set link_mode to 'logit' or 'probit'.")}

          if(link_mode == "logit"){
            link_mode_numeric = 1
          } 

          if(link_mode == "probit"){
            link_mode_numeric = 2
          } 

         }

         if(outcome_mode=="poisson"){
          if(link_mode != "log"){stop("If outcome_mode is 'poisson', you must set link_mode to 'log'.")}
          outcome_mode_numeric = 3
          link_mode_numeric = 3
         }

         if(outcome_mode=="gaussian"){
          if(link_mode != "identity"){stop("If outcome_mode is 'gaussian', you must set link_mode to 'identity'.")}
          outcome_mode_numeric = 4
          link_mode_numeric = 4
         }

         if(is.null(outcome_mode_numeric)) stop("outcome_mode not supported")

         N_id = dim(self_report[[1]])[1]

         if(!is.null(m_e_data)){ 
          if(!all(names(m_e_data)==c("sampled", "sampled_exposure", "sampled_mask", "detected", "detected_exposure", "detected_mask")))
          stop("m_e_data must be named list: sampled, sampled_exposure, sampled_mask, detected, detected_exposure, detected_mask, in this order.")
         }else{
          m_e_data = list(sampled=rep(0,N_id), sampled_exposure=rep(0,N_id), sampled_mask=rep(0,N_id), detected=rep(0,N_id), detected_exposure=rep(0,N_id), detected_mask=rep(0,N_id))
         }

         layer_names = names(self_report)

         # Check self_report data
         if(is.null(exposure)){
         if(is.null(self_report)) stop("outcome must be a list of matrices.")
         if(!is.list(self_report)) stop("outcome must be a list of matrices.")
         for(i in 1:length(self_report)){
          if(outcome_mode=="bernoulli"){
          if(!all(self_report[[i]] %in% c(0,1,NA))) stop("self_report must be binary 0 or 1")
          }
         }
         }

         if(!is.null(exposure)){
           if(length(self_report) != length(exposure)) stop("outcome and exposure must be lists of matrices equal in length.")
           if(sum(names(exposure)==names(outcome)) != length(names(exposure))) stop("Names of exposure and outcome must match. Order matters.")
         }

        if(!is.null(mask)){
           if(length(self_report) != length(mask)) stop("outcome and mask must be lists of matrices equal in length.")
           if(sum(names(mask)==names(outcome)) != length(names(mask))) stop("Names of mask and outcome must match. Order matters.")
         }

         # Check ground_truth data
         if(!is.null(ground_truth)){ 
         if(!is.list(ground_truth)) stop("ground_truth must be a list of matrices.")
         }

         # Check block_covariates data
         if(!is.null(block_covariates)){ 
         if(!is.data.frame(block_covariates)) stop("block_covariates must be a data frame.")
             for(i in 1:dim(block_covariates)[2]){
                 if(!is.factor(block_covariates[,i])) stop("block_covariates must be factor variables.")
             }     
         }

         # Check individual_covariates data
         if(!is.null(individual_covariates)){ 
         if(!is.data.frame(individual_covariates)) stop("individual_covariates must be a data frame.")
         }

         # Check self_report data
         if(!is.null(dyadic_covariates)){ 
         if(!is.list(dyadic_covariates)) stop("dyadic_covariates must be a list of matrices.")
         } 
          
         ###################### Check rownames and colnames
         if(check_data_organization==TRUE){
          node_names = rownames(self_report[[1]])

         for(i in 1:length(self_report)){
          if(is.null(rownames(self_report[[i]]))) stop("All outcome layers must have rownames.")
          if(is.null(colnames(self_report[[i]]))) stop("All outcome layers must have colnames.")
          if(all(rownames(self_report[[i]]) == node_names)==FALSE) stop("All outcome matrices must have the same colnames.")
          if(all(colnames(self_report[[i]]) == node_names)==FALSE) stop("All outcome matrices must have the same rownames.")
         }
        
         if(!is.null(dyadic_covariates)){ 
         for(i in 1:length(dyadic_covariates)){
          if(is.null(rownames(dyadic_covariates[[i]]))) stop("All dyadic layers must have rownames.")
          if(is.null(colnames(dyadic_covariates[[i]]))) stop("All dyadic layers must have colnames.")
          if(all(rownames(dyadic_covariates[[i]]) == node_names)==FALSE) stop("All dyadic matrices must have the same colnames and match outcomes.")
          if(all(colnames(dyadic_covariates[[i]]) == node_names)==FALSE) stop("All dyadic matrices must have the same rownames and match outcomes.")
         }}

         if(!is.null(mask)){ 
         for(i in 1:length(mask)){
          if(is.null(rownames(mask[[i]]))) stop("All mask layers must have rownames.")
          if(is.null(colnames(mask[[i]]))) stop("All mask layers must have colnames.")
          if(all(rownames(mask[[i]]) == node_names)==FALSE) stop("All mask matrices must have the same colnames and match outcomes.")
          if(all(colnames(mask[[i]]) == node_names)==FALSE) stop("All mask matrices must have the same rownames and match outcomes.")
         }}

         if(!is.null(exposure)){ 
         for(i in 1:length(exposure)){
          if(is.null(rownames(exposure[[i]]))) stop("All exposure layers must have rownames.")
          if(is.null(colnames(exposure[[i]]))) stop("All exposure layers must have colnames.")
          if(all(rownames(exposure[[i]]) == node_names)==FALSE) stop("All exposure matrices must have the same colnames and match outcomes.")
          if(all(colnames(exposure[[i]]) == node_names)==FALSE) stop("All exposure matrices must have the same rownames and match outcomes.")
         }}

         if(!is.null(individual_covariates)){ 
          if(is.null(rownames(individual_covariates))) stop("individual_covariates must have rownames.")
          if(all(rownames(individual_covariates) == node_names)==FALSE) stop("individual_covariates must have same rownames as outcome matrix.")
         }

         if(!is.null(block_covariates)){ 
          if(is.null(rownames(block_covariates))) stop("block_covariates must have rownames.")
          if(all(rownames(block_covariates) == node_names)==FALSE) stop("block_covariates must have same rownames as outcome matrix.")
         }
         }

         ############################################################################# Process data
         N_id =  dim(self_report[[1]])[1]
         N_responses = length(self_report)

         outcomes = array(NA, c(N_id, N_id, N_responses))
         mask_mat = array(NA, c(N_id, N_id, N_responses))

         if(is.null(mask)){
            for(i in 1:length(self_report)){
              mask_mat[,,i] = matrix(0, nrow=N_id, ncol=N_id)
             }} else{
            for(i in 1:length(self_report)){
              mask_mat[,,i] = mask[[i]]
             }
         }


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
         group_ids_character = array("Any", c(N_id, 1))
         group_ids = block_covariates
         group_ids_levels = "No Blocks"
          } else{
          
         N_block_types = length(block_covariates[1,]) 
         N_groups_per_type = rep(NA, N_block_types)
         group_ids_character = array(NA, c(N_id, N_block_types))
         group_ids = array(NA, c(N_id, N_block_types))
         group_ids_levels = vector("list", N_block_types)

         for(i in 1:N_block_types){
          N_groups_per_type[i] = max(as.numeric(block_covariates[,i]),na.rm=TRUE)
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

         if(check_standardization == TRUE){
            for(i in 1:dim(individual_covariates)[2]){
              if(is_numeric_not_binary(individual_covariates[,i])==TRUE){
                sd_test = sd(individual_covariates[,i],na.rm=TRUE)
                if(sd_test > standardization_threshold[1] & sd_test < standardization_threshold[2]){
                    # All good
                  } else{
                    stop(paste0("Column ", i, " in individual_covariates has not been standardized. Use standardize_strand(). \n Setting: 'check_standardization = FALSE' turns off this error check. \n Predictors in STRAND models should generally be standardized, igonore at your own risk."))
                  }
              }
            }
          }


         }

        if(is.null(dyadic_covariates)){
          N_dyadic_predictors = 0
          dyadic_predictors = 0
         } else{
          N_dyadic_predictors = length(dyadic_covariates)
          dyadic_predictors = dyadic_covariates

          if(check_standardization == TRUE){
            for(i in 1:length(dyadic_covariates)){
              if(is_numeric_not_binary(dyadic_covariates[[i]])==TRUE){
                sd_test = sd(dyadic_covariates[[i]],na.rm=TRUE)
                if(sd_test > standardization_threshold[1] & sd_test < standardization_threshold[2]){
                    # All good
                  } else{
                    stop(paste0("Element ", i, " in dyadic_covariates has not been standardized. Use standardize_strand(). \n Setting: 'check_standardization = FALSE' turns off this error check. \n Predictors in STRAND models should generally be standardized, igonore at your own risk."))
                  }
              }
            }
          }


         }

         ############################################################################# Determine legal models
         if(multiplex==FALSE){
          if(longitudinal==TRUE){
            supported_models = c("Longitudinal")
            } else{
         if(N_responses==1){
          if(is.null(m_e_data)){
            supported_models = c("SRM", "SBM", "SRM+SBM")
            } else{
            supported_models = c("SRM", "SBM", "SRM+SBM","SRM+SBM+ME")
            }
         } 

         if(N_responses==2 & N_networktypes==2){
          supported_models = c("LNM")
         }

         if(N_responses==2 & N_networktypes==3){
          supported_models = c("LNM","LNM+Flows")
         }
         }} else{
          if(N_responses > 1){
           supported_models = c("Multiplex")
          }
         }

   ############################################ Diffusion models
   if(is.null(diffusion_outcome)){
    diffusion_outcome = rep(0, N_id)
         } 

   if(is.null(diffusion_exposure)){
    diffusion_exposure = rep(1, N_id)
         } 

   if(is.null(diffusion_mask)){
    diffusion_mask = rep(0, N_id)
         } 

   ############################################ Check if directed
    directed_flags = rep(NA, N_responses)
    for(q in 1:N_responses){
     directed_flags[q] = ifelse(all(outcomes[,,q] == t(outcomes[,,q])) & sum(outcomes[,,q])>0, "undirected", "directed")
         }
     directed_flag = ifelse(all(directed_flags == "directed"), "directed", "undirected") 

      if(directed_flag=="directed" & directed == FALSE){
        stop("All outcome networks are directed, but 'directed' argument set to FALSE.")
      }

      if(directed_flag=="undirected" & directed == TRUE){
       warning("At least one outcome network is undirected, but 'directed' argument set to TRUE. 
       Consider setting 'directed=FALSE' so that STRAND performs relevant checks.
       If outcome networks are undirected, then dyadic covariates should be symmetric, and focal and target regressions equivalent.") 
      }
      
      if(N_dyadic_predictors>0){
          dyadic_directed_flags = rep(NA, N_dyadic_predictors)
        for(q in 1:N_dyadic_predictors){
          dyadic_directed_flags[q] = ifelse(all(dyadic_predictors[[q]] == t(dyadic_predictors[[q]])), "undirected", "directed") 
        }
        dyadic_directed_flag = ifelse(any(dyadic_directed_flags == "directed"), "directed", "undirected") 

      if(dyadic_directed_flag=="directed" & directed == FALSE){
        warning("At least one dyadic covariate layer is directed, but 'directed' argument set to FALSE. 
        If outcome networks are undirected, then dyadic covariates should typically be symmetric too. Rethink your model.") 
      }

      }
      

   ###################################### Merge all
   model_dat = list(
     N_networktypes = N_networktypes,                                               
     N_id = N_id,                                                                                                          
     N_responses = N_responses,  
     N_periods=N_periods,                                              
     N_individual_predictors = N_individual_predictors,                                       
     N_dyadic_predictors =  N_dyadic_predictors,                                           
     outcomes = outcomes,  
     outcomes_real = outcomes,  
     flows=flows,                         
     individual_predictors = individual_predictors,      
     dyadic_predictors = dyadic_predictors,
     N_block_predictors = N_block_types,
     N_groups_per_block_type = N_groups_per_type,
     block_predictors = group_ids,
     outcome_mode=outcome_mode_numeric,
     link_mode=link_mode_numeric,
     exposure=exposure_risk,
     mask=mask_mat,
     sampled = m_e_data$sampled, 
     detected = m_e_data$detected,
     sampled_exposure = m_e_data$sampled_exposure, 
     detected_exposure = m_e_data$detected_exposure,
     sampled_mask = m_e_data$sampled_mask, 
     detected_mask = m_e_data$detected_mask,
     diffusion_outcomes = diffusion_outcome,
     diffusion_exposure = diffusion_exposure,
     diffusion_mask = diffusion_mask,
     imputation = ifelse(imputation==TRUE, 1, 0)
     )
    
    if(imputation == FALSE){
    for(q in 1:length(model_dat)){
      if(sum(is.na(model_dat[[q]]))>0){
       stop(paste0("Variable: ", names(model_dat)[q], " contains missing values. Set: imputation=TRUE, or rebuild your data objects."))
      }
     }
    }

   attr(model_dat, "class") = "STRAND Data Object"
   attr(model_dat, "supported_models") = supported_models
   attr(model_dat, "layer_names") = layer_names
   attr(model_dat, "group_ids_character") = group_ids_character
   attr(model_dat, "group_ids_levels") = group_ids_levels
   attr(model_dat, "directed") = ifelse(directed == TRUE, "directed", "undirected")
   colnames(attr(model_dat, "group_ids_character")) = colnames(model_dat$block_predictors)
   names(attr(model_dat, "group_ids_levels")) = colnames(attr(model_dat, "group_ids_character"))
   
  return(model_dat)
}
