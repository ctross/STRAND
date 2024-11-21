#' Organize Stan output and provide summaries of variance parameters
#' 
#' This is a function to organize Stan output and provide summaries of key model parameters
#'
#' @param input A STRAND model object with social relations model parameters.
#' @param n_partitions Should variance on latent scale be partioned into 3 factors (focal, target, dyadic+error) as in amen? Or 4 factors (focal, target, dyadic, error)?
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @param include_samples An indicator for the user to specify if raw samples, or only the summary statistics should be returned. Samples can take up a lot of space.
#' @param include_reciprocity Should reciprocity estimates be returned?.
#' @param mode Should the dyadic correlation "cor" be computed, or the dyadic covariance "cov", or the adjusted dyadic+error correlation "adj".
#' @param plot_pairs A Boolean to plot pairs of the variance terms.
#' @param plot_type Should the raw variances be plotted: "raw", or should the variance partition coefficients: "vpcs"?
#' @param plot_layer The name of the layer to plot.
#' @param verbose If TRUE prints the results.
#' @return A STRAND results object.
#' @export
#' @examples
#' \dontrun{
#' res = strand_VPCs(input=fit)
#' }
#'

strand_VPCs = function(input, n_partitions = 4, HPDI=0.9, include_reciprocity=FALSE, mode="cor", include_samples=FALSE, plot_pairs=FALSE, plot_type = "vpcs", plot_layer = NULL, verbose=FALSE){
    if(attributes(input)$class != "STRAND Model Object"){
        stop("strand_VPCs() requires a fitted object of class: STRAND Model Object.")
    }

    if(attributes(input)$fit_type != "mcmc"){
      if(attributes(input)$fit_type == "vb"){
         warning("Final, publication-ready model fits for STRAND models should always be produced using MCMC! Variational inference via Pathfinder can be used in Stan
              during experimental model runs, but final inferences should be based on MCMC sampling. In our tests, Pathfinder results are decently similar to MCMC results, 
              but often failed to recover strong true effects. ")  
         } else{
         stop("Fitted results can only be reorganized for STRAND model objects fit using MCMC. Variational inference or optimization can be used in Stan
              during experimental model runs, but final inferences should be based on MCMC sampling.")  
      }    
    }

    if(!n_partitions %in% c(3,4)){stop("n_partitions must be 3 or 4")}

    ###################################################### Create samples 
    N_responses = input$data$N_responses
    layer_names = attr(input$data,"layer_names")

    if(attr(input, "model_type") %in% c("LNM","LNM+Flows")){
       N_responses = 1
       layer_names = layer_names[1]
    }

    stanfit = posterior::as_draws_rvars(input$fit$draws())

    if(input$data$link_mode==1){
       base_sd = sqrt(0.33333 * (3.14159^2))
    }

    if(input$data$link_mode==2){
       base_sd = 1
    }

    if(input$data$link_mode==3){
       stop("Automatic VPCs not yet deployable for Poisson models.")
    }


    ################### Network model parameters
    sr_sigma = posterior::draws_of(stanfit$"sr_sigma")
    dr_sigma = posterior::draws_of(stanfit$"dr_sigma")

    G_corr = posterior::draws_of(stanfit$"G_corr")
    D_corr = posterior::draws_of(stanfit$"D_corr")
    lims = c(-1,1)

    #################################### Prep for alternate
     if(mode %in% c("cov", "adj")){
     new = D_corr
     dr_sigma = posterior::draws_of(stanfit$"dr_sigma")

    if(input$data$link_mode==1){
       base_sd = sqrt(0.33333 * (3.14159^2))
    }

    if(input$data$link_mode==2){
       base_sd = 1
    }
    
    if(mode == "cov"){
      for(q in 1:dim(D_corr)[1]){
       new[q,,] = diag(rep(dr_sigma[q,],2)) %*% D_corr[q,,] %*% diag(rep(dr_sigma[q,],2))
      }

     med_cov = apply(new, 2:3, median) 
     lims = lims*max(abs(c(min(med_cov), max(med_cov))))
    } 

    if(mode == "adj"){
      if(input$data$link_mode==3){
       stop("Automatic adjustment not yet deployable for Poisson models.")
       }

      for(q in 1:dim(D_corr)[1]){
        sigma_scrap = rep(dr_sigma[q,],2)
        sigma_scrap = sigma_scrap / sqrt(sigma_scrap^2 + rep(base_sd, length(sigma_scrap))^2) # Check this is right
        new[q,,] = diag(sigma_scrap) %*% D_corr[q,,] %*% diag(sigma_scrap)
      }
    }  
     D_corr = new
    }

    #############

    srm_samples = list(
            focal_target_sd=sr_sigma,
            dyadic_sd = dr_sigma,
            G_corr=G_corr,
            D_corr=D_corr
        )

    samples = list(srm_model_samples=srm_samples)

    ######################################################### Get VPCs
    results_list_vcs3 = list()
    results_list_vcs4 = list()

    results_list_vpcs3 = list()
    results_list_vpcs4 = list()

    ######### Calculate all variance componants
    error_sd = rep(base_sd, length(samples$srm_model_samples$dyadic_sd[,1]))

     for(q in 1:N_responses){
      ################################### Three-way  
      vcs = cbind(samples$srm_model_samples$focal_target_sd[,q]^2, samples$srm_model_samples$focal_target_sd[,N_responses+q]^2, samples$srm_model_samples$dyadic_sd[,q]^2 + error_sd^2)
      vcs_total = rowSums(vcs)
      vcst = cbind(vcs_total,vcs_total,vcs_total)

      vcs_normalized = vcs/vcst

      results_list_vcs3[[q]] = vcs
      results_list_vpcs3[[q]] = vcs_normalized

      ################################### Four-way  
      vcs = cbind(samples$srm_model_samples$focal_target_sd[,q]^2, samples$srm_model_samples$focal_target_sd[,N_responses+q]^2, samples$srm_model_samples$dyadic_sd[,q]^2, error_sd^2)
      vcs_total = rowSums(vcs)
      vcst = cbind(vcs_total,vcs_total,vcs_total,vcs_total)

      vcs_normalized = vcs/vcst

      results_list_vcs4[[q]] = vcs
      results_list_vpcs4[[q]] = vcs_normalized
     }


    ###################################################### Create summary stats  
     results_list_3_vcs = list()
     results_list_4_vcs = list()

     results_list_3_vpcs = list()
     results_list_4_vpcs = list()

    ################### SRM model
     results_vcs3 = matrix(NA, nrow=3, ncol=7)
     results_vpcs3 = matrix(NA, nrow=3, ncol=7)

     results_vcs4 = matrix(NA, nrow=4, ncol=7)
     results_vpcs4 = matrix(NA, nrow=4, ncol=7)

    ######### Calculate all corr and block effects
    layer_names_long = c(paste0(layer_names, " (out)"), paste0(layer_names, " (in)"))
    layer_names_long2 = c(paste0(layer_names, " (i to j)"), paste0(layer_names, " (j to i)"))
     results_srm_base = matrix(NA, nrow = (2*N_responses*(2*N_responses-1)), ncol=7)
     ticker = 0

     for(m in 1:((N_responses*2)-1)){
     for(n in (m+1):(N_responses*2)){
        ticker = ticker + 1 
        results_srm_base[ticker,] = sum_stats(paste0("Generalized reciprocity - ", layer_names_long[m], " - ", layer_names_long[n]), samples$srm_model_samples$G_corr[,m,n], HPDI)
        }}

     for(m in 1:((N_responses*2)-1)){
     for(n in (m+1):(N_responses*2)){
        ticker = ticker + 1 
        results_srm_base[ticker,] = sum_stats(paste0("Dyadic reciprocity - ", layer_names_long2[m], " - ", layer_names_long2[n]), samples$srm_model_samples$D_corr[,m,n], HPDI)
        }}

    ######### Calculate all focal effects
     for(q in 1:N_responses){
      results_vcs3[1,] = sum_stats(paste0("focal var - ", layer_names[q]), results_list_vcs3[[q]][,1], HPDI)
      results_vcs3[2,] = sum_stats(paste0("target var - ", layer_names[q]), results_list_vcs3[[q]][,2], HPDI)
      results_vcs3[3,] = sum_stats(paste0("dyadic+error var - ", layer_names[q]), results_list_vcs3[[q]][,3], HPDI)

      results_vpcs3[1,] = sum_stats(paste0("focal var - ", layer_names[q]), results_list_vpcs3[[q]][,1], HPDI)
      results_vpcs3[2,] = sum_stats(paste0("target var - ", layer_names[q]), results_list_vpcs3[[q]][,2], HPDI)
      results_vpcs3[3,] = sum_stats(paste0("dyadic+error var - ", layer_names[q]), results_list_vpcs3[[q]][,3], HPDI)

      results_vcs4[1,] = sum_stats(paste0("focal var - ", layer_names[q]), results_list_vcs4[[q]][,1], HPDI)
      results_vcs4[2,] = sum_stats(paste0("target var - ", layer_names[q]), results_list_vcs4[[q]][,2], HPDI)
      results_vcs4[3,] = sum_stats(paste0("dyadic var - ", layer_names[q]), results_list_vcs4[[q]][,3], HPDI)
      results_vcs4[4,] = sum_stats(paste0("error var - ", layer_names[q]), results_list_vcs4[[q]][,4], HPDI)

      results_vpcs4[1,] = sum_stats(paste0("focal var - ", layer_names[q]), results_list_vpcs4[[q]][,1], HPDI)
      results_vpcs4[2,] = sum_stats(paste0("target var - ", layer_names[q]), results_list_vpcs4[[q]][,2], HPDI)
      results_vpcs4[3,] = sum_stats(paste0("dyadic var - ", layer_names[q]), results_list_vpcs4[[q]][,3], HPDI)
      results_vpcs4[4,] = sum_stats(paste0("error var - ", layer_names[q]), results_list_vpcs4[[q]][,4], HPDI)

      colnames(results_vcs3) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 
      colnames(results_vpcs3) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 
      colnames(results_vcs4) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 
      colnames(results_vpcs4) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 
      colnames(results_srm_base) = c("Variable", "Median", paste("HPDI", (1-HPDI)/2, sep=":"), paste("HPDI", (1+HPDI)/2, sep=":"), "Mean","SD","P") 

     results_list_3_vcs[[q]] = results_vcs3
     results_list_4_vcs[[q]]  = results_vcs4

     results_list_3_vpcs[[q]]  = results_vpcs3
     results_list_4_vpcs[[q]]  = results_vpcs4
     }

   ############# Finally, merge all effects into a list
     names(results_list_3_vcs) = layer_names
     names(results_list_4_vcs) = layer_names

     names(results_list_3_vpcs) = layer_names
     names(results_list_4_vpcs) = layer_names

     if(include_reciprocity==FALSE){
      results_srm_base = NULL
     }
          
   results_out_3 = list( "Raw Variance Componants" = results_list_3_vcs, "Variance Partition Coefficients" = results_list_3_vpcs, "Reciprocity" = results_srm_base)
   results_out_4 = list( "Raw Variance Componants" = results_list_4_vcs, "Variance Partition Coefficients" = results_list_4_vpcs, "Reciprocity" = results_srm_base)
   
   if(n_partitions==3){
    if(verbose==TRUE){
     print(results_out_3)
     }

      if(include_samples==TRUE){
       results_out_3$RawSamples = results_list_vcs3
       results_out_3$VPCsSamples = results_list_vpcs3
      }

    res_final = results_out_3
   }

   if(n_partitions==4){
    if(verbose==TRUE){
     print(results_out_4)
     }

     if(include_samples==TRUE){
      results_out_4$RawSamples = results_list_vcs4
      results_out_4$VPCsSamples = results_list_vpcs4
     }

    res_final = results_out_4
   }

   if(plot_pairs==TRUE){
    if(include_samples==FALSE) {stop("To plot pairs, please set include_samples to TRUE.")}
    if(plot_type == "vpcs"){
      pairs(res_final$VPCsSamples[[which(layer_names==plot_layer)]])
    }

    if(plot_type == "raw"){
      pairs(res_final$RawSamples[[which(layer_names==plot_layer)]])
    }
    } 
 
    return(res_final)
     
 }

