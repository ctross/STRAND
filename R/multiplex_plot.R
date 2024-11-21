#' Plot multiplex random effects in matrix form
#'
#' @param input A STRAND model object, obtained by fitting a multiplex combined stochastic block and social relations model.
#' @param type Plot dyadic reciprocity ("dyadic") or generalized reciprocity ("generalized")?.
#' @param mode For dyadic plots, should the dyadic correlation "cor" be plotted, or the dyadic covariance "cov", or the adjusted dyadic+error correlation "adj".
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
#' res = multiplex_plot(input = fit)
#' }
#'

multiplex_plot = function(input, type="dyadic", mode="cor", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = NULL, height=6, width=6, palette=NULL){
 if(!type %in% c("dyadic", "generalized", "dyad", "general")) stop("type nust be 'dyadic' or 'generalized'. Or, a short form: 'dyad', or 'general'.")

 if(!attr(input, "model_type") %in% c("Longitudinal", "Multiplex")) stop("multiplex_plot() can only be fit to models fit using  fit_multiplex_model() or fit_longitudinal_model().")

 if(type=="dyadic" | type=="dyad"){
  multiplex_plot_d(input=input, HPDI=HPDI, mode=mode, plot = plot, export_as_table = export_as_table, save_plot = save_plot, height=height, width=width, palette=palette)
 }

 if(type=="generalized" | type=="general"){
 multiplex_plot_g(input=input, HPDI=HPDI, plot = plot, export_as_table = export_as_table, save_plot = save_plot, height=height, width=width, palette=palette)
 }

}
