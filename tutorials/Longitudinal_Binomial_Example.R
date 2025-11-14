#####################################################################################################
#
#   Longitudinal Network Analyses - Binomial models  
#
#####################################################################################################

 library(stringr)
 library(ggplot2)
 library(psych)

# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors


############################################## 
# Load data
 data(Baboon_Longitudinal_Data)
 d = Baboon_Longitudinal_Data

############################################## 
# Process data into longitudinal form
dat_long = NULL
# Loop over y days of data, making a day-specific network data-set
# Some days are missing outcome data. This is dealt with via the mask layer
# if outcome data is missing, mask[i,j]=1
# currently, missings in the predictors aren't supported in STRAND, but will be eventually
for(y in 1:14){
    d$Individual$Age = standardize_strand(d$Individual$Age)

    if(all(d$Mask[[y]]==1)){
    Presenting = t(d$Presenting[[y]])               # If all entries are masked, then standardize_strand() fails 
    } else{
    Presenting = standardize_strand(t(d$Presenting[[y]]))  # If there are  some data, then standardize  
    }

 # Merge data
 dat_long[[y]] = make_strand_data(
  outcome = list("Affiliative" = d$Affiliative[[y]]),
  exposure = list("Affiliative" = d$Exposure[[y]]),
  mask = list("Affiliative" = d$Mask[[y]]),
  block_covariates = NULL, 
  individual_covariates = d$Individual, 
  dyadic_covariates = list("Presenting" = Presenting),
  longitudinal = TRUE,
  outcome_mode="binomial",
  link_mode="logit",
  directed = TRUE
  )
 }
names(dat_long) = paste("Time", c(1:14))

############################################## 
# Fit model with time-varying slopes
fit_2a = fit_longitudinal_model(
 long_data=dat_long,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ Presenting,
 coefficient_mode="varying",
 random_effects_mode="fixed",
 bandage_penalty = -1,
  mode="mcmc",
 mcmc_parameters = list(seed = 1, chains = 1, init=0,
    parallel_chains = 1, refresh = 1, iter_warmup = 1000,
    iter_sampling = 1000, max_treedepth = 12)
  )

res_2a = summarize_longitudinal_bsrm_results(fit_2a)

############################################## 
# Visualize results
pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,2))
longitudinal_plot(fit_2a, type="dyadic", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,2,9,3,4))
longitudinal_plot(fit_2a, type="generalized", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,7,5,10,8,6))
longitudinal_plot(fit_2a,type="coefficient", 
    parameter_set = list(
    focal="Age", target="Age", 
    focal="SexMale", target="SexMale",
    dyadic="Presenting"),
    palette=pal,
    normalized=TRUE)


############################################## 
# Fit model with time-invariant slopes
fit_2b = fit_longitudinal_model(
 long_data=dat_long,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ Presenting,
 coefficient_mode="fixed",
 random_effects_mode="fixed",
 bandage_penalty = -1,
  mode="mcmc",
 mcmc_parameters = list(seed = 1, chains = 1, init=0,
    parallel_chains = 1, refresh = 1, iter_warmup = 100,
    iter_sampling = 1000, max_treedepth = 12),
 priors=NULL
  )

res_2b = summarize_longitudinal_bsrm_results(fit_2b)

############################################## 
# Visualize results
pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,2))
longitudinal_plot(fit_2b, type="dyadic", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,2,9,3,4))
longitudinal_plot(fit_2b, type="generalized", palette = pal, height=6, width=6.5)


############################################## 
# Make a merged plot
bab2b_set = c("focal effects coeffs (out-degree), Time 1 - Age", 
"focal effects coeffs (out-degree), Time 1 - SexMale", 
"target effects coeffs (in-degree), Time 1 - Age", 
"target effects coeffs (in-degree), Time 1 - SexMale",
"dyadic effects coeffs, Time 1 - Presenting")

to_add = res_2b$summary[which(res_2b$summary$Variable %in% bab2b_set),]

    parameter_set = list(
    focal="Age", target="Age", 
    focal="SexMale", target="SexMale",
    dyadic="Presenting")

base_set = longitudinal_plot_c(fit_2a, parameters=as.vector(unlist(parameter_set)), type=names(parameter_set),
    plot=FALSE,
    normalized=FALSE,
    export_as_table=TRUE)

to_add$short_names = c("Focal - Time 1 - Age", "Focal - Time 1 - SexMale",
                "Target - Time 1 - Age", "Target - Time 1 - SexMale",
                 "Dyadic - Time 1 - Presenting")
to_add$time_point = rep("Time 0", 5)
to_add$time_point_int = rep(0, 5)
to_add$extra_short_names = c("Focal - Age", "Focal - SexMale",
                "Target - Age", "Target - SexMale",
                 "Dyadic - Presenting")

to_add$type_set = c("Focal", "Focal", "Target", "Target", "Dyadic")

to_add$Model = "Fixed"
base_set$Model = "Varying"

full_set = rbind(base_set,to_add)


full_set2 = full_set[,c("Variable", "time_point", "extra_short_names", "extra_short_names", "Median", "HPDI:0.05", "HPDI:0.95", "Mean", "SD", "time_point_int", "type_set")]
colnames(full_set2) = c("Variable", "Layer", "Target" , "Base" , "Median", "L", "H", "Mean","SD","LayerNumeric","Type")
      
     Diff = as.numeric(full_set2$H)-as.numeric(full_set2$L)   
     full_set2$Median = as.numeric(full_set2$Median)/Diff
     full_set2$L = as.numeric(full_set2$L)/Diff
     full_set2$H =  as.numeric(full_set2$H)/Diff

     full_set2$Model = ifelse(full_set2$"LayerNumeric"==0, "Fixed", "Varying")
     full_set2$Model2 = ifelse(full_set2$"Model"=="Fixed", "Fixed", full_set2$Type)

     full_set2$Model2 = factor(full_set2$Model2)
     full_set2$Model2 = factor(full_set2$Model2, levels=c("Fixed", "Dyadic", "Focal", "Target"))

p = ggplot(full_set2, aes(x=LayerNumeric, y=as.numeric(Median), ymin=as.numeric(L), ymax=as.numeric(H), group=Target, color=Target))+ 
     geom_linerange(size=1, position = position_dodge(width = 0.3)) + facet_grid(.~Model2, scales="free", space="free") + 
     geom_point(size=2, position = position_dodge(width = 0.3))+
     geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
     labs(y="Effect size", x="Time step") + 
     theme(strip.text.x = element_text(size=12,face="bold"), 
      strip.text.y = element_text(size=12,face="bold"),
      axis.text = element_text(size=12),
      axis.title = element_text(size=14, face="bold"))+
     theme(strip.text.y = element_text(angle = 360)) + 
    # coord_flip() + 
     theme(panel.spacing = grid::unit(1, "lines")) + scale_color_manual(values = pal) + 
     theme(legend.position="bottom") + theme(legend.title = element_blank()) + scale_x_continuous(breaks=1:14,expand = c(0, 0.95))

p

