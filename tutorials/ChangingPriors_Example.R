##############################################################################################
#
#   Poisson Analysis with Custom Priors 
#
##############################################################################################

# Load libraries
library(STRAND)
library(ggplot2)

# Load data
data(Bat_Data)

# Number of minutes of blood lickings
nets = list(Lick = round(Bat_Data$Lick/60,0))

# Dyadic variables
dyad = list(Relatedness = standardize_strand(Bat_Data$Relatedness), 
            NoOpportunity = Bat_Data$NoOpportunity
              )

# Block variables
group_ids = data.frame(Sex = as.factor(Bat_Data$Individual$Sex))
rownames(group_ids) = rownames(Bat_Data$Individual)

model_dat = make_strand_data(outcome = nets,
                              block_covariates = group_ids, 
                              individual_covariates = NULL, 
                              dyadic_covariates = dyad,
                              outcome_mode = "poisson",
                              link_mode = "log"
                              )

# Default Priors Model
fit = fit_block_plus_social_relations_model(data=model_dat,
                                            block_regression = ~ Sex ,
                                            focal_regression = ~ 1,
                                            target_regression = ~ 1,
                                            dyad_regression = ~ NoOpportunity * Relatedness,
                                            mode="mcmc",
                                            mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                                        max_treedepth = 11, adapt_delta = 0.95)
                                           )
res = summarize_strand_results(fit)

# Custom, Extra-Vague Priors Model
fit2 = fit_block_plus_social_relations_model(data=model_dat,
                                            block_regression = ~ Sex ,
                                            focal_regression = ~ 1,
                                            target_regression = ~ 1,
                                            dyad_regression = ~ NoOpportunity * Relatedness,
                                            mode="mcmc",
                                            mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                                        max_treedepth = 11, adapt_delta = 0.95),
                                            priors=make_priors(
                                             priors_to_change=list(
                                               "B_ingroup"=c(0.01, 3.5), 
                                               "B_outgroup"=c(0.01, 3.5),
                                               "focal_effects"=c(0, 10), 
                                               "target_effects"=c(0, 10), 
                                               "dyad_effects"=c(0, 10),
                                               "sr_sigma"=c(2, 15), 
                                               "dr_sigma"=c(2, 15) 
                                             )
                                           )
                                          )
res2 = summarize_strand_results(fit2)

##############################################################################################
# Prep results for plotting
tab1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), 
                                 normalized=TRUE,  only_slopes=TRUE, export_as_table=TRUE)
tab2 = strand_caterpillar_plot(res2, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), 
                                 normalized=TRUE,  only_slopes=TRUE, export_as_table=TRUE)

tab1$Priors = "Default, Weakly regularizing"
tab2$Priors = "Modified, Especially vague"

tab0 = rbind(tab1, tab2)

plot_0 = ggplot(tab0, aes(x = Variable, y = Median, 
        ymin = LI, ymax = HI, color=Priors)) + geom_linerange(linewidth = 1,  position = position_dodge(0.33)) + 
        geom_point(size = 2, position = position_dodge(0.33)) + facet_grid(Submodel ~ 
        ., scales = "free", space = "free") + geom_hline(aes(yintercept = 0), 
        color = "black", linetype = "dashed") + labs(y = "Regression parameters", 
        x = "") + theme(strip.text.x = element_text(size = 12, 
        face = "bold"), strip.text.y = element_text(size = 12, 
        face = "bold"), axis.text = element_text(size = 12), 
        axis.title.y = element_text(size = 14, face = "bold"), 
        axis.title.x = element_blank()) + theme(strip.text.y = element_text(angle = 360)) + 
        coord_flip() + theme(panel.spacing = grid::unit(1, 
        "lines")) + theme(legend.position="bottom") + scale_color_manual(values=c("#E69F00", "black"))
plot_0

# Priors dont affect the results much here

