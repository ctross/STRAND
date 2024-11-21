#' An internal function to make a longitudinal data object using the STRAND framework
#' 
#' @param long_data A list of data objects of class STRAND prepared using the make_strand_data() function. The data objects must include all covariates used in the formulas listed below.
#' @param block_regression A formula for the block-level predictors. This should be specified as in lm(), e.g.: ~ Ethnicity + Sex. Dont use interactions, however.
#' @param focal_regression A formula for the predictors of out-degree (i.e., focal effects, or the effects of individual covariates on outgoing ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param target_regression A formula for the predictors of in-degree (i.e., target effects, or the effects of individual covariates on incoming ties). This should be specified as in lm(), e.g.: ~ Age * Education
#' @param dyad_regression A formula for the predictors of dyadic relationships. This should be specified as in lm(), e.g.: ~ Kinship + Friendship
#' @return A STRAND data object.
#' @export
#' @examples
#' \dontrun{
#' fit = make_longitudinal_data(long_data,
#'                              block_regression = ~ Ethnicity,
#'                              focal_regression = ~ Age * NoFood,
#'                              target_regression = ~ Age * NoFood,
#'                              dyad_regression = ~ Relatedness + Friends * SameSex
#'                               )
#' }
#' 

make_longitudinal_data = function(long_data,
                                  block_regression,
                                  focal_regression,
                                  target_regression,
                                  dyad_regression
                                  ){
    ############################################################################# Run checks that each STRAND data object matches one another
     if(is.null(names(long_data))){
      stop("Please ensure that long_data is a named list of STRAND data objects.")
     }

     run_longitudinal_checks(long_data, pass=1)

    ############################################################################# Resolve data list into STRAND data object
    thin_dat = NULL
     for(i in 1:length(long_data)){
      thin_dat[[i]] = parse_longitudinal_data(
      data=long_data[[i]],
      block_regression = block_regression,
      focal_regression = focal_regression,
      target_regression = target_regression,
      dyad_regression = dyad_regression
     )
     }

     run_longitudinal_checks(thin_dat, pass=2)

    dyad_set = array(NA, c(dim(thin_dat[[1]]$dyad_set),length(thin_dat)))
    focal_set = array(NA, c(dim(thin_dat[[1]]$focal_set),length(thin_dat)))
    target_set = array(NA, c(dim(thin_dat[[1]]$target_set),length(thin_dat)))
    block_set  = array(NA, c(dim(thin_dat[[1]]$block_set),length(thin_dat)))

    outcomes  = array(NA, c(dim(thin_dat[[1]]$outcomes)[c(1:2)],length(thin_dat)))
    exposure  = array(NA, c(dim(thin_dat[[1]]$exposure)[c(1:2)],length(thin_dat)))
    mask  = array(NA, c(dim(thin_dat[[1]]$mask)[c(1:2)],length(thin_dat)))

    for(i in 1:length(thin_dat)){
     dyad_set[ , , , i] = thin_dat[[i]]$dyad_set
     focal_set[ , , i] = thin_dat[[i]]$focal_set
     target_set[ , , i] = thin_dat[[i]]$target_set
     block_set[ , , i]  = thin_dat[[i]]$block_set

     outcomes[ , , i]  = thin_dat[[i]]$outcomes
     exposure[ , , i]  = thin_dat[[i]]$exposure
     mask[ , , i]  = thin_dat[[i]]$mask
    }

    dimnames(dyad_set) = dimnames(thin_dat[[1]]$dyad_set)
    dimnames(focal_set) = dimnames(thin_dat[[1]]$focal_set)
    dimnames(target_set) = dimnames(thin_dat[[1]]$target_set)
    dimnames(block_set) = dimnames(thin_dat[[1]]$block_set)

    dimnames(outcomes) = dimnames(thin_dat[[1]]$outcomes)
    dimnames(exposure) = dimnames(thin_dat[[1]]$exposure)
    dimnames(mask) = dimnames(thin_dat[[1]]$mask)

    dat_out = thin_dat[[1]]

    dat_out$dyad_set = dyad_set
    dat_out$focal_set = focal_set
    dat_out$target_set = target_set
    dat_out$block_set = block_set
    dat_out$outcomes = outcomes
    dat_out$exposure = exposure
    dat_out$mask = mask
    
    attr(dat_out, "class") = "STRAND Data Object"
    attr(dat_out, "supported_models") = c("Longitudinal")
    attr(dat_out, "layer_names") = names(long_data)

    data = dat_out
    data$N_responses = length(long_data)

    return(data)
}
