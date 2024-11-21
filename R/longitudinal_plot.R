#' Plot longitudinal random effects in time-lag form
#'
#' @param input A STRAND model object, obtained by fitting a longitudinal combined stochastic block and social relations model.
#' @param type Plot dyadic reciprocity ("dyadic") or generalized reciprocity ("generalized")?. Can use "custom", if the user passes in an appropriate data frame.
#' @param mode For dyadic plots, should the dyadic correlation "cor" be plotted, or the dyadic covariance "cov", or the adjusted dyadic+error correlation "adj".
#' @param parameter_set A labeled list of paramters to plot. E.g.: parameter_set = list(focal="Age", target="Age", focal="Sex", target="Sex",  dyadic = "Relatedness"). The 
#' names must be focal, target, or dyadic, and the quoted text must be paramters given in the regression equations. E.g., focal_regression = ~ Age + Sex.
#' @param results A results data.frame or matrix with colnames: c("Variable", "Layer", "Target" , "Base" , "Median", "L", "H", "Mean","SD","LayerNumeric") 
#' @param normalized Should effects be normalized? Only valid for type="coefficient".
#' @param HPDI Highest Posterior Density Interval. Ranges in (0,1).
#' @param plot Should a plot be displayed?
#' @param export_as_table Should the tabular data rather than a ggplot object be returned?
#' @param save_plot Should a plot be exported to working directory? If so, set save_plot="desired_filename.pdf".
#' @param height Height of exported plot.
#' @param width Width of exported plot.
#' @param palette Override the default palette with a vector of color codes. "dyadic" need 3 colors, "generalized" needs 5 colors, and "coefficient" needs as many colors as supplied parameters.
#' @return A figure or tabluar data to make a figure.
#' @export
#' @examples
#' \dontrun{
#' res = longitudinal_plot(input = fit)
#' }
#'

longitudinal_plot = function(input, type="dyadic", mode="cor", parameter_set, results, normalized = FALSE,  HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = NULL, height=6, width=6, palette=NULL){
 if(!type %in% c("dyadic", "generalized", "coefficient", "dyad", "general", "coeff", "custom")) stop("type must be: 'coefficient', 'dyadic' or 'generalized'. Or, a short form: 'dyad', 'general', 'coeff'. Last option is: 'custom'. ")

 if(attr(input, "model_type") != "Longitudinal") stop("longitudinal_plot() can only be fit to models fit using fit_longitudinal_model().")

 if(type=="dyadic" | type=="dyad"){
  if(attr(input, "random_effects_mode") != "fixed") stop("Dyadic longitudinal_plot() can only be fit to models fit using fit_longitudinal_model() with: random_effects_mode='fixed'.")
   res = longitudinal_plot_d(input=input, HPDI=HPDI, mode=mode, plot = plot, export_as_table = export_as_table, save_plot = save_plot, height=height, width=width, palette=palette)
 }

 if(type=="generalized" | type=="general"){
  if(attr(input, "random_effects_mode") != "fixed") stop("Generalized longitudinal_plot() can only be fit to models fit using fit_longitudinal_model() with: random_effects_mode='fixed'.")
   res = longitudinal_plot_g(input=input, HPDI=HPDI, plot = plot, export_as_table = export_as_table, save_plot = save_plot, height=height, width=width, palette=palette)
 }

 if(type=="coefficient" | type=="coeff"){
  if(attr(input, "coefficient_mode") != "varying") stop("Coefficient longitudinal_plot() can only be fit to models fit using fit_longitudinal_model() with: coefficient_mode='varying'.")
   res = longitudinal_plot_c(input=input, parameters=as.vector(unlist(parameter_set)), type=names(parameter_set), normalized=normalized, HPDI=HPDI, plot = plot, export_as_table = export_as_table, save_plot = save_plot, height=height, width=width, palette=palette)
 }

 if(type=="custom"){
  if(attr(input, "coefficient_mode") != "varying") stop("Coefficient longitudinal_plot() can only be fit to models fit using fit_longitudinal_model() with: coefficient_mode='varying'.")
   res = longitudinal_plot_b(input=input, results=results, plot = plot, save_plot = save_plot, height=height, width=width, palette=palette)
 }

 return(res)

}
