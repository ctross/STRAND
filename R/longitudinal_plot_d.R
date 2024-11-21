#' Plot longitudinal random effects - dyadic reciprocity
#'
#' @param input A STRAND model object, obtained by fitting a longitudinal combined stochastic block and social relations model.
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @param mode Should the dyadic correlation "cor" be plotted, or the dyadic covariance "cov", or the adjusted dyadic+error correlation "adj".
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
#' res = longitudinal_plot_d(input = fit)
#' }
#'

longitudinal_plot_d = function(input, HPDI=0.9, mode="cor", plot = TRUE, export_as_table = FALSE, save_plot = NULL, height=6, width=6, palette=NULL){
 if(is.null(palette)){
  palette = c("#975e3d", "black", "#406e6b")  
    } 

 stanfit = posterior::as_draws_rvars(input$fit$draws())
 corr = posterior::draws_of(stanfit$"D_corr")
 lims = c(-1,1)

  if(mode %in% c("cov", "adj")){
  new = corr
  dr_sigma = posterior::draws_of(stanfit$"dr_sigma")

    if(input$data$link_mode==1){
       base_sd = sqrt(0.33333 * (3.14159^2))
    }

    if(input$data$link_mode==2){
       base_sd = 1
    }
    
    if(mode == "cov"){
      for(q in 1:dim(corr)[1]){
       new[q,,] = diag(rep(dr_sigma[q,],2)) %*% corr[q,,] %*% diag(rep(dr_sigma[q,],2))
      }

     med_cov = apply(new, 2:3, median) 
     lims = lims*max(abs(c(min(med_cov), max(med_cov))))
    } 

    if(mode == "adj"){
      if(input$data$link_mode==3){
       stop("Automatic adjustment not yet deployable for Poisson models.")
       }

      for(q in 1:dim(corr)[1]){
        sigma_scrap = rep(dr_sigma[q,],2)
        sigma_scrap = sigma_scrap / sqrt(sigma_scrap^2 + rep(base_sd, length(sigma_scrap))^2) # Check this is right
        new[q,,] = diag(sigma_scrap) %*% corr[q,,] %*% diag(sigma_scrap)
      }
    }  
   corr = new
    }

 N_responses = input$data$N_responses
 layer_names = attr(input$data,"layer_names")

 names_outcomes = c(paste0(layer_names, "\n(i to j)"), paste0(layer_names, "\n(j to i)"))

 rs_m = apply(corr, 2:3, median)
 rs_l = apply(corr, 2:3, HPDI, prob=HPDI)[1,,]
 rs_h = apply(corr, 2:3, HPDI, prob=HPDI)[2,,]

 rs_m[lower.tri(rs_m)] = NA
 diag(rs_m) = NA

 rs_l[lower.tri(rs_l)] = NA
 diag(rs_l) = NA

 rs_h[lower.tri(rs_h)] = NA
 diag(rs_h) = NA

  m1 = m2 = rs_type = rs_m
 for(i in 1:(2*N_responses)){
   m1[i,] = names_outcomes[i]
   m2[,i] = names_outcomes[i]
 }

 rs_type[which(!is.na(rs_type))] = "Between-node"

 for(m in 1:(N_responses-1)){
 for(n in (m+1):N_responses){
  rs_type[m,n] = "Within-node"
  rs_type[n,m+N_responses] = NA
  rs_type[N_responses+m,N_responses+n] = NA

  rs_m[n,m+N_responses] = NA
  rs_m[N_responses+m,N_responses+n] = NA

  rs_l[n,m+N_responses] = NA
  rs_l[N_responses+m,N_responses+n] = NA

  rs_h[n,m+N_responses] = NA
  rs_h[N_responses+m,N_responses+n] = NA
 }}

 for(m in 1:N_responses){
  rs_type[m,m+N_responses] = "Dyadic Reciprocity"
 }

 rs_m = c(rs_m)
 rs_l = c(rs_l)
 rs_h = c(rs_h)

 m1 = c(m1)
 m2 = c(m2)

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

 r_if_sig = ifelse(rs_l > 0 | rs_h < 0, round(rs_m, 2), NA)
 Time2 = Time = l = h = c()

 df = data.frame(rs_m=rs_m, l=rs_l, h=rs_h, m1=m1, m2=m2, r_if_sig=r_if_sig, rs_type=rs_type)

 df_sub = df[which(df$m1==names_outcomes[1]),]

 df_sub$Time = c(layer_names,layer_names)
 df_sub$Time2 = c(1:length(layer_names),1:length(layer_names))-1

 df_sub = df_sub[which(!is.na(df_sub$rs_m)),]

p = ggplot2::ggplot(df_sub, ggplot2::aes(x=Time2, y=as.numeric(rs_m), ymin=as.numeric(l), ymax=as.numeric(h), group=rs_type, color=rs_type))+ 
     ggplot2::geom_linerange(size=1, position = ggplot2::position_dodge(width = 0.3))+
     ggplot2::geom_point(size=2, position = ggplot2::position_dodge(width = 0.3))+
     ggplot2::geom_hline(ggplot2::aes(yintercept=0),color="black",linetype="dashed")+
     ggplot2::labs(y="Dyadic reciprocity correlations", x="Time-step lags") + 
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
  return(df_sub)   
 }

}



