#' Plot longitudinal random effects - dyadic reciprocity
#'
#' @param input A STRAND model object, obtained by fitting a longitudinal combined stochastic block and social relations model.
#' @param results A results data.frame or matrix with colnames: c("Variable", "Layer", "Target" , "Base" , "Median", "L", "H", "Mean","SD","LayerNumeric") 
#' @param plot Should a plot be displayed?
#' @param save_plot Should a plot be exported to working directory? If so, set save_plot="desired_filename.pdf".
#' @param height Height of exported plot.
#' @param width Width of exported plot.
#' @param palette Override the default palette with a 3-vector of color codes.
#' @return A figure or tabluar data to make a figure.
#' @export
#'

longitudinal_plot_b = function(input, results, plot = TRUE, save_plot = NULL, height=6, width=6, palette=NULL){
  H = L = LayerNumeric = Median = Target = NULL
 
 if(is.null(palette)){
  palette = terrain.colors(nrow(results)) 
    } 

  colnames(results) = c("Variable", "Layer", "Target" , "Base" , "Median", "L", "H", "Mean","SD","LayerNumeric")   
   
p = ggplot2::ggplot(results, ggplot2::aes(x=LayerNumeric, y=as.numeric(Median), ymin=as.numeric(L), ymax=as.numeric(H), group=Target, color=Target))+ 
     ggplot2::geom_linerange(size=1, position = ggplot2::position_dodge(width = 0.3)) + 
     ggplot2::geom_point(size=2, position = ggplot2::position_dodge(width = 0.3))+
     ggplot2::geom_hline(ggplot2::aes(yintercept=0),color="black",linetype="dashed")+
     ggplot2::labs(y="Effect size", x="Time step") + 
     ggplot2::theme(strip.text.x = ggplot2::element_text(size=12,face="bold"), 
      strip.text.y = ggplot2::element_text(size=12,face="bold"),
      axis.text = ggplot2::element_text(size=12),
      axis.title = ggplot2::element_text(size=14, face="bold"))+
     ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
    # ggplot2::coord_flip() + 
     ggplot2::theme(panel.spacing = grid::unit(1, "lines")) + ggplot2::scale_color_manual(values = palette) + 
     ggplot2::theme(legend.position="bottom") + ggplot2::theme(legend.title = ggplot2::element_blank())

 if(!is.null(save_plot)){
  ggplot2::ggsave(save_plot, p, height=height, width=width)
 }

 if(plot == TRUE){
  plot(p)
  return(p)  
 }

}



