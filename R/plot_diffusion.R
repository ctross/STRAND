#' A function to simulate a network-based diffusion proccess using the STRAND framework
#' 
#' This function allows users to simulate data using a NBDA model. The user must supply a list of STRAND data objects, a set of parameters,
#' and a series of formulas following standard lm() style syntax. 
#'
#' @param long_data A list of data objects of class STRAND prepared using the make_strand_data() function. The data objects must include all covariates and trait diffusion data used in the formulas listed below.
#' @param palette A 2-vector of plot colors. First element is for the low value of scale, second is for the high value.
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' plot_diffusion(long_data)
#' }
#' 

plot_diffusion = function(long_data, palette = c("#13667B", "#35ABC0")){
  X = Y = Z = NULL
  data = make_longitudinal_data(long_data = long_data,
                               block_regression = ~ 1,
                               focal_regression = ~ 1,
                               target_regression = ~ 1,
                               dyad_regression = ~ 1
                                  )

 diffusion_outcomes = matrix(NA, nrow = data$N_id, ncol = data$N_responses)

  for(t in 1:data$N_responses){
   diffusion_outcomes[,t] = long_data[[t]]$diffusion_outcomes 
    }

  # Heatmap 
   cols_diffusion_outcomes = rows_diffusion_outcomes = matrix(NA, nrow=data$N_id, ncol=data$N_responses) 

   for(i in 1:data$N_id){
    rows_diffusion_outcomes[i,] = c(1:data$N_responses)
   }

   for(i in 1:data$N_responses){
    cols_diffusion_outcomes[,i] = c(1:data$N_id)
   }

  plot_dat = data.frame(X = c(rows_diffusion_outcomes), Y=c(cols_diffusion_outcomes), Z=c(diffusion_outcomes))

  # Viz
  ggplot2::ggplot(plot_dat, ggplot2::aes(X, Y, fill = Z)) + 
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient(low=palette[1], high=palette[2]) +
    ggplot2::theme(legend.position="none") + ggplot2::xlab("Time-step") + ggplot2::ylab("Individual")

}
