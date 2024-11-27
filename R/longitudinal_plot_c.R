#' Plot longitudinal random effects - dyadic reciprocity
#'
#' @param input A STRAND model object, obtained by fitting a longitudinal combined stochastic block and social relations model.
#' @param parameters Name of the parameter as included in the function call. E.g., ~ Age + Sex... Then: parameters=c("Age", "Age", "Sex", "Sex")
#' @param type For each parameter, indicate its type: "focal", "target", or "dyadic". And, from the above: type = =c("focal", "target", "focal", "target")
#' @param normalized Should effects be normalized?
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @param plot Should a plot be displayed?
#' @param export_as_table Should the tabular data rather than a ggplot object be returned?
#' @param save_plot Should a plot be exported to working directory? If so, set save_plot="desired_filename.pdf".
#' @param height Height of exported plot.
#' @param width Width of exported plot.
#' @param palette Override the default palette with a 3-vector of color codes.
#' @return A figure or tabluar data to make a figure.
#' @export
#' @examples
#' \dontrun{
#' res = longitudinal_plot_c(input = fit, parameters, type)
#' }
#'

longitudinal_plot_c = function(input, parameters, type, normalized = FALSE, HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = NULL, height=6, width=6, palette=NULL){
  layer_names = attr(input$data,"layer_names")
  res = summarize_longitudinal_bsrm_results(input, quiet = TRUE)

  #### Make names
  long_names = c()
  short_names = c() 
  time_point = c() 
  time_point_int = c() 
  extra_short_names = c() 
  Median = l = h = c() 
  type_set = c() 
  name_types = c("focal effects coeffs (out-degree), ", "target effects coeffs (in-degree), ", "dyadic effects coeffs, ")
  ticker = 0
 
  for(i in 1:length(parameters)){
  for(q in 1:length(layer_names)){
   ticker = 1 + ticker
  
    if(type[i]=="focal"){
     long_names[ticker] = paste0(name_types[1], layer_names[q], " - ", parameters[i] )
     short_names[ticker] = paste0("Focal - ", layer_names[q], " - ", parameters[i] )
     time_point[ticker] = layer_names[q] 
     time_point_int[ticker] = q
     extra_short_names[ticker] = paste0("Focal - ", parameters[i] )
     type_set[ticker] = "Focal"
    }

    if(type[i]=="target"){
     long_names[ticker] = paste0(name_types[2], layer_names[q], " - ", parameters[i] )
     short_names[ticker] = paste0("Target - ", layer_names[q], " - ", parameters[i] )
     time_point[ticker] = layer_names[q] 
     time_point_int[ticker] = q
     extra_short_names[ticker] = paste0("Target - ", parameters[i] )
     type_set[ticker] = "Target"
    }

    if(type[i]=="dyadic"){
     long_names[ticker] = paste0(name_types[3], layer_names[q], " - ", parameters[i] )
     short_names[ticker] = paste0("Dyadic - ", layer_names[q], " - ", parameters[i] )
     time_point[ticker] = layer_names[q] 
     time_point_int[ticker] = q
     extra_short_names[ticker] = paste0("Dyadic - ", parameters[i] )
     type_set[ticker] = "Dyadic"
    }
  }}

 if(is.null(palette)){
  palette = terrain.colors(length(parameters)) 
    } 
   
   res_short = res$summary

   in_set = data.frame(Variable=long_names, short_names, time_point, time_point_int, extra_short_names, type_set) 

   df = merge(res_short, in_set, by="Variable")

   df2 = df 

   colnames(df2) = c("Variable", "Median", "l", "h", "Mean", "SD", "P", "short_names", "time_point", "time_point_int", "extra_short_names", "type")


   df2$Diff = as.numeric(df2$h)-as.numeric(df2$l)   

   if(normalized==TRUE){
     df2$Median = as.numeric(df2$Median)/df2$Diff
     df2$l = as.numeric(df2$l)/df2$Diff
     df2$h =  as.numeric(df2$h)/df2$Diff
     }


p = ggplot2::ggplot(df2, ggplot2::aes(x=time_point_int, y=as.numeric(Median), ymin=as.numeric(l), ymax=as.numeric(h), group=extra_short_names, color=extra_short_names))+ 
     ggplot2::geom_linerange(size=1, position = ggplot2::position_dodge(width = 0.3)) + ggplot2::facet_grid(cols = ggplot2::vars(type)) +
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

 if(export_as_table == TRUE){
  return(df)   
 }

}



