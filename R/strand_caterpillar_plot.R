#' A function to create a STRAND data object
#'
#' This function organizes network data and covariates into a form that can be used by STRAND for model fitting. All 
#' STRAND model fitting functions require their data to be supplied in the format exported here.
#'
#' @param 
#' results The results object from summarize_strand_results.
#' @param
#' submodels Which submodels to plot? Supported: "False positive rate", "Recall of true ties", "Theta: question-order effects", "Focal effects: Out-degree", "Target effects: In-degree", "Dyadic effects", "Other estimates"
#' @param 
#' normalized Should effects be normalized? Do not use for correlations and variance terms.
#' @param 
#' only_slopes If TRUE then technical parameters are omited, and only slopes and intercepts are plotted.
#' @param 
#' only_technicals If TRUE then technical parameters are plotted, and slopes and intercepts are omited
#' @param 
#' site A string to indicate group or feildsite.
#' @param 
#' export_as_table If TRUE, the data are not plotted, and a data.frame is returned. This is useful for making complex ggplot objects with multiple groups/sites.
#' @return A plot or dataframe
#' @export
#' @examples
#' \dontrun{
#' vis = strand_caterpillar_plot(res, submodel=c("Focal efffects: Out-degree",
#'                                               "Target effects: In-degree",
#'                                               "Dyadic effects","Other estimates"), 
#'                                    normalized=TRUE, 
#'                                    site="XY", 
#'                                    only_technicals=TRUE, 
#'                                    export_as_table=FALSE
#'                                     )
#' }
#'

#plotting results
strand_caterpillar_plot = function(results, submodels=NULL, normalized=FALSE, only_slopes=FALSE, only_technicals=FALSE, site="BOB", export_as_table=FALSE){
  Variable = Median = LI = HI = NULL
  dat = vector("list",length(results$summary_list))

  for(k in 1:length(results$summary_list)){
   dat[[k]] = data.frame(results$summary_list[[k]])
   dat[[k]]$SubModel = names(results$summary_list)[k]
   colnames(dat[[k]]) = c("Variable", "Median", "LI", "HI", "Mean","SD", "SubModel")
   for(j in 2:6)
   dat[[k]][,j] = as.numeric(dat[[k]][,j])
  }


df = do.call(rbind, dat)

colnames(df) = c("Variable", "Median", "LI", "HI", "Mean","SD", "SubModel")


df$Submodel = factor(df$SubModel)
df$Submodel = factor(df$SubModel, levels=c("False positive rate", "Recall of true ties","Theta: question-order effects",
                                           "Focal effects: Out-degree","Target effects: In-degree","Dyadic effects", "Other estimates" ))

if(only_slopes==TRUE){
exclude=c("false positive rate intercept, layer 1",               
"false positive rate intercept, layer 2",               
"false positive rate sd, layer 1",                     
"false positive rate sd, layer 2",                            
"recall rate of true ties intercept, layer 1",          
"recall rate of true ties intercept, layer 2",          
"recall rate of true ties sd, layer 1",                 
"recall rate of true ties sd, layer 2",                 
"theta intercept, layer 1 to 2",                        
"theta sd, layer 1 to 2",                               
"focal effects sd",                                            
"target effects sd",                                            
"dyadic effects sd",                                                         
"focal-target effects rho (generalized recipocity)",    
"dyadic effects rho (dyadic recipocity)")   

df = df[which(!df$Variable %in% exclude),]
}

if(only_technicals==TRUE){
include=c("false positive rate intercept, layer 1",               
"false positive rate intercept, layer 2",               
"false positive rate sd, layer 1",                     
"false positive rate sd, layer 2",                            
"recall rate of true ties intercept, layer 1",          
"recall rate of true ties intercept, layer 2",          
"recall rate of true ties sd, layer 1",                 
"recall rate of true ties sd, layer 2",                 
"theta intercept, layer 1 to 2",                        
"theta sd, layer 1 to 2",                               
"focal effects sd",                                            
"target effects sd",                                            
"dyadic effects sd",                                                         
"focal-target effects rho (generalized recipocity)",    
"dyadic effects rho (dyadic recipocity)")   

unit=c("false positive rate intercept, layer 1",               
"false positive rate intercept, layer 2",                                  
"recall rate of true ties intercept, layer 1",          
"recall rate of true ties intercept, layer 2",                          
"theta intercept, layer 1 to 2",                                                                                                             
"focal-target effects rho (generalized recipocity)",    
"dyadic effects rho (dyadic recipocity)")

df = df[which(df$Variable %in% include),]

df$Scaling = ifelse(df$Variable %in% unit, "Rates", "Dispersion")
df$SubModel2 = ifelse(df$SubModel =="Other estimates", "Correlation", "Dispersion")
}


if(!is.null(submodels))
df = df[which(df$SubModel %in% submodels),]

df$Diff = df$HI-df$LI   

if(normalized==TRUE) {
  df$Median = df$Median/df$Diff
  df$LI = df$LI/df$Diff
  df$HI =  df$HI/df$Diff
}
 
 df$Site = site

df$Variable[which(df$Variable=="focal-target effects rho (generalized recipocity)")] = "focal-target effects rho"
df$Variable[which(df$Variable=="dyadic effects rho (dyadic recipocity)")] = "dyadic effects rho"

df$Submodel = factor(df$SubModel)
df$Submodel = factor(df$SubModel, levels=c("False positive rate", "Recall of true ties","Theta: question-order effects",
                                           "Focal effects: Out-degree","Target effects: In-degree","Dyadic effects", "Other estimates" ))

p = ggplot2::ggplot(df, ggplot2::aes(x=Variable,y=Median,ymin=LI,ymax=HI))+ 
     ggplot2::geom_linerange(size=1)+
     ggplot2::geom_point(size=2)+
     ggplot2::facet_grid( Submodel ~ ., scales = "free", space='free')+
       #facet_wrap(vars(SubModel), ncol=1,scales = "free")+
     ggplot2::geom_hline(ggplot2::aes(yintercept=0),color="black",linetype="dashed")+
     ggplot2::labs(y="Regression parameters", x="") + 
     ggplot2::theme(strip.text.x = ggplot2::element_text(size=12,face="bold"), 
      strip.text.y = ggplot2::element_text(size=12,face="bold"),
      axis.text = ggplot2::element_text(size=12),
      axis.title.y = ggplot2::element_text(size=14, face="bold"), 
      axis.title.x = ggplot2::element_blank())+
     ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
     ggplot2::coord_flip() + 
     ggplot2::theme(panel.spacing = unit(1, "lines")) 

p2 = ggplot2::ggplot(df, ggplot2::aes(x=Variable,y=Median,ymin=LI,ymax=HI))+ 
     ggplot2::geom_linerange(size=1)+
     ggplot2::geom_point(size=2)+
     ggplot2::facet_grid( SubModel2~., scales = "free")+
       #facet_wrap(vars(SubModel), ncol=1,scales = "free")+
     ggplot2::geom_hline(ggplot2::aes(yintercept=0),color="black",linetype="dashed")+
     ggplot2::labs(y="Regression parameters", x="") + 
     ggplot2::theme(strip.text.x = ggplot2::element_text(size=12,face="bold"), 
      strip.text.y = ggplot2::element_text(size=12,face="bold"),
      axis.text=ggplot2::element_text(size=12),
      axis.title.y=ggplot2::element_text(size=14, face="bold"), 
      axis.title.x=ggplot2::element_blank())+
     ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
     ggplot2::coord_flip() + 
     ggplot2::theme(panel.spacing = unit(1, "lines")) 

if(export_as_table==FALSE){
if(only_technicals==TRUE){
 return(p2)
 } else{
    return(p)
 }
}

if(export_as_table==TRUE){
  df
  }


 }
