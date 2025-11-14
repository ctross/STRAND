#' Plot multiplex random effects - generalized reciprocity
#'
#' @param input A STRAND model object, obtained by fitting a multiplex combined stochastic block and social relations model.
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
#' res = multiplex_plot_g(input = fit)
#' }
#'

multiplex_plot_g = function(input, HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = NULL, height=6, width=6, palette=NULL){
 if(is.null(palette)){
  palette = c("#7D370D", "#FBFEF9", "#114B47")  
    }

  if(attr(input, "fit_type")=="numpyro"){
         samps = convert_posterior(input$fit$get_samples())
         corr = samps$G_corr 
    }else{
         stanfit = posterior::as_draws_rvars(input$fit$draws())
         corr = posterior::draws_of(stanfit$"G_corr")
      }

 N_responses = input$data$N_responses
 layer_names = attr(input$data,"layer_names")

 rs_m = apply(corr, 2:3, median)
 rs_l = apply(corr, 2:3, HPDI, prob=HPDI)[1,,]
 rs_h = apply(corr, 2:3, HPDI, prob=HPDI)[2,,]

 rs_m[lower.tri(rs_m)] = NA
 diag(rs_m) = NA

 rs_l[lower.tri(rs_l)] = NA
 diag(rs_l) = NA

 rs_h[lower.tri(rs_h)] = NA
 diag(rs_h) = NA

 rs_type = rs_m
 rs_type[which(!is.na(rs_type))] = "Cross"

 for(m in 1:(N_responses-1)){
 for(n in (m+1):N_responses){
  rs_type[m,n] = "Sender"
  rs_type[N_responses+m,N_responses+n] = "Receiver"
 }}

 rs_m = c(rs_m)
 rs_l = c(rs_l)
 rs_h = c(rs_h)

 rs_type = c(rs_type)  

 # Prep for CI figure
  substrRight = function(x, n){
  x = as.character(x)
  substr(x, nchar(x)-n+1, nchar(x))
  }

  substrLeft = function(x, n){
  x = as.character(x)
  substr(x, 1, nchar(x)-n+1)
  }

 names_outcomes = c(paste0(layer_names, "\n(sender)"), paste0(layer_names, "\n(receiver)"))

 measure1 = factor(rep(names_outcomes, each=(N_responses*2)))
 measure2 = factor(rep(names_outcomes, (N_responses*2)))

 measure1 = factor(measure1,names_outcomes)
 measure2 = factor(measure2,rev(names_outcomes))

 r_if_sig = ifelse(rs_l > 0 | rs_h < 0, round(rs_m, 2), NA)

 df = data.frame(rs_m=rs_m, l=rs_l, h=rs_h, measure1=measure1, measure2=measure2, r_if_sig=r_if_sig, rs_type=rs_type)


 df$Group1 = ifelse(substrRight(df$measure1,4)=="der)" & substrRight(df$measure2,4)=="der)", "Sender",
                   ifelse(substrRight(df$measure1,4)=="ver)" & substrRight(df$measure2,4)=="ver)", "Receiver",
                   "Cross"
                   ))

 p1 = ggplot2::ggplot(df, ggplot2::aes(measure1, measure2, fill=rs_m, label=round(r_if_sig,2),color=rs_type)) +
 ggplot2::geom_tile(ggplot2::aes(width=0.96, height=0.96), size=1.69) +
 ggplot2::labs(x = NULL, y = NULL, fill = "Correlation", title="Generalized reciprocity estimates", subtitle="Only reliable correlation coefficients shown") + 
 ggplot2::scale_fill_gradient2(mid=palette[2],low=palette[1],high=palette[3], limits=c(-1,1)) +
 ggplot2::geom_text(color="black",size=12*0.36) +
 ggplot2::theme_classic() +
 ggplot2::scale_x_discrete(expand=c(0,0)) +
 ggplot2::scale_y_discrete(expand=c(0,0)) + ggplot2::guides(color = "none") + ggplot2::theme(plot.title = ggplot2::element_text(size = 14))  +
 ggplot2::scale_color_manual(values=c("grey75","grey15","grey15","grey75")) + ggplot2::theme(axis.text = ggplot2::element_text(size = 14))

 if(!is.null(save_plot)){
  ggplot2::ggsave(save_plot, p1, height=height, width=width)
 }

 if(plot == TRUE){
  plot(p1)
  return(p1)  
 }

 if(export_as_table == TRUE){
  return(df)   
 }

}
