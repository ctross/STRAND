#' Calculate the contrast between two block estimates for each layer
#' 
#' This is a helper function to compute contrast between samples
#'
#' @param input A STRAND model object 
#' @param focal The outcome of interest, will be of the form "X to Y", where X and Y are factor levels from the block model.
#' @param base The base case, will be of the form "X to Y", where X and Y are factor levels from the block model.
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @return A matrix of output.
#' @export
#'

process_block_parameters = function(input, focal, base, HPDI=0.9){
     if(attributes(input)$class != "STRAND Model Object"){
        stop("process_block_parameters() requires a fitted object of class: STRAND Model Object.")
    }

    if(attributes(input)$fit_type != "mcmc"){
        stop("Fitted results can only be processed for STRAND model objects fit using MCMC. Variational inference or optimization can be used in Stan
              during experimental model runs, but final inferences should be based on MCMC sampling.")   
    }

    ###################################################### Create samples 
    N_responses = input$data$N_responses
    layer_names = attr(input$data,"layer_names")

    stanfit = posterior::as_draws_rvars(input$fit$draws())

    ################### Block model parameters
    if(dim(input$data$block_set)[2]>0)
    block_effects_raw = posterior::draws_of(stanfit$"block_effects")

    if(attr(input, "model_type") %in% c("LNM","LNM+Flows")){
       N_responses = 1
       layer_names = layer_names[1]
    }

    if(attr(input, "model_type") %in% c("SBM", "SRM+SBM","SRM+SBM+ME","LNM","LNM+Flows")){
       dims = dim(block_effects_raw)
       dims = c(dims[1],1,dims[2])
       block_effects = array(NA, dims)
       block_effects[,1,] = block_effects_raw
    }

    if(attr(input, "model_type") %in% c("Longitudinal","Multiplex")){
      # Nothing to do  
      block_effects = block_effects_raw
    }

    if(!attr(input, "model_type") %in% c("Longitudinal","Multiplex","SBM", "SRM+SBM","SRM+SBM+ME","LNM","LNM+Flows")){
       stop("process_block_parameters() requires a fitted object of class: STRAND Model Object.")
    }

    ################### Get index data for block-model samples
    block_indexes = c()
    block_indexes[1] = 0
    for(q in 1:input$data$N_group_vars){ 
    block_indexes[1+q] = input$data$N_groups_per_var[q]*input$data$N_groups_per_var[q] + block_indexes[q]
    }

    ################### Convert the block-model effects into an array form
    B_multi = list()
    for(l in 1:N_responses){
    B = list()
    for(q in 1:input$data$N_group_vars){
      B[[q]] = array(NA, c(dim(block_effects)[1], input$data$N_groups_per_var[q], input$data$N_groups_per_var[q]  ))

      for(s in 1:dim(block_effects)[1]){
       B[[q]][s,,] =  array(block_effects[s,l,(block_indexes[q]+1):(block_indexes[q+1])], c(input$data$N_groups_per_var[q], input$data$N_groups_per_var[q]))
      }
    } 
    B_multi[[l]] = B  
    }

    names(B_multi) = layer_names

    srm_samples = list(
            block_parameters=B_multi
        )

    samples = list(srm_model_samples=srm_samples)

    ###################################################### Create summary stats 
     sum_stats_2 = function(y1, y2, y3, y4, y5, x, z){
      bob = rep(NA, 11)
       dig = 3
      bob[1] = y1
      bob[2] = y2
      bob[3] = y3
      bob[4] = y4
      bob[5] = round(median(x),dig)
      bob[6] = round(HPDI(x, z)[1],dig)
      bob[7] = round(HPDI(x, z)[2],dig)
      bob[8] = round(mean(x),dig)
      bob[9] = round(sd(x),dig)
      bob[10] = y5
      bob[11] = round(bayesian_p(x),dig)

      return(bob)
      }
     
     results_list = list()

    ######### Calculate all block effects
    results_srm_base = matrix(NA, nrow = dim(block_effects)[1], ncol = N_responses*dim(block_effects)[3])
    names_of_block_vars = rep(NA, N_responses*dim(block_effects)[3])
    names_of_layers = rep(NA, N_responses*dim(block_effects)[3])
     ticker = 0

     group_ids_character_df = cbind(rep("Any",input$data$N_id),attr(input$data, "group_ids_character"))

    if(is.null(colnames(group_ids_character_df))){
        colnames(group_ids_character_df) = paste0("(NoBlockingVars)", 1:ncol(group_ids_character_df))
     }

     if(is.null(dim(attr(input$data, "group_ids_character")))){
        colnames(group_ids_character_df) = c("(Intercept)","Scrap")
     }
     
     colnames(group_ids_character_df)[1] = "(Intercept)"
     in_IDs = colnames(input$data$block_set)
     all_IDs = colnames(group_ids_character_df)
     group_ids_character_df = group_ids_character_df[,match(in_IDs, all_IDs), drop = FALSE]

     group_id_levels = append("Any", attr(input$data, "group_ids_levels"), 1)
     names(group_id_levels)[1]= "(Intercept)"
     
     for(l in 1:N_responses){
       B_set = samples$srm_model_samples$block_parameters[[l]]
     for(q in 1:input$data$N_group_vars){
      group_ids_character = levels(factor(group_ids_character_df[,q]))
      test_sorting = group_id_levels[[which(names(group_id_levels) == colnames(group_ids_character_df)[q])]]
      if(all(group_ids_character==test_sorting)==FALSE){
        stop("Factors not sorted correctly.")
      }

      for(b1 in 1:input$data$N_groups_per_var[q]){
      for(b2 in 1:input$data$N_groups_per_var[q]){
       ticker = ticker + 1  

      names_of_block_vars[ticker] = paste0(group_ids_character[b1], " to ", group_ids_character[b2])
      names_of_layers[ticker] = layer_names[l]
      results_srm_base[,ticker] = B_set[[q]][,b1,b2]
                                                                         
     }}

     }}
  ################ Process
  results_srm_out = matrix(NA, nrow=N_responses, ncol=11)

     for(l in 1:N_responses){
       base_wave = results_srm_base[,which(names_of_layers == layer_names[l])]
       nbv = names_of_block_vars[which(names_of_layers == layer_names[l])]

       contrast_wave = base_wave[,which(nbv == focal)] - base_wave[,which(nbv == base)] 

      results_srm_out[l,] = sum_stats_2(paste0("Contrast: ", focal, ", versus ", base), layer_names[l], focal, base, l, contrast_wave, HPDI)

     }   
    
    colnames(results_srm_out) = c("Variable", "Layer", "Category" , "Base" , "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","LayerNumeric","P") 


   ############# Finally, return
    return(results_srm_out)
}


