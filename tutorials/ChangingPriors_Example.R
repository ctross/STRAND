########################################
#
#   Poisson Analysis with Custom Priors 
#
########################################

# Clear working space
rm(list = ls()) 

# Load libraries
library(STRAND)

# Load data
data(Bat_Data)

# Number of minutes of blood lickings
nets = list(Lick = round(Bat_Data$Lick/60,0))

# Dyadic variables
dyad = list(Relatedness = Bat_Data$Relatedness, 
            NoOpportunity = Bat_Data$NoOpportunity
              )

# Block variables
group_ids = data.frame(Sex = as.factor(Bat_Data$Sex))

model_dat = make_strand_data(self_report = nets,
                              block_covariates = group_ids, 
                              individual_covariates = NULL, 
                              dyadic_covariates = dyad,
                              outcome_mode = "poisson"
                              )

# Default Priors Model
fit = fit_block_plus_social_relations_model(data=model_dat,
                                            block_regression = ~ Sex ,
                                            focal_regression = ~ 1,
                                            target_regression = ~ 1,
                                            dyad_regression = ~ NoOpportunity + Relatedness,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                                        max_treedepth = NULL, adapt_delta = .98)
                                           )
res = summarize_strand_results(fit)

# Custom, Extra-Vague Priors Model
fit2 = fit_block_plus_social_relations_model(data=model_dat,
                                            block_regression = ~ Sex ,
                                            focal_regression = ~ 1,
                                            target_regression = ~ 1,
                                            dyad_regression = ~ NoOpportunity + Relatedness,
                                            mode="mcmc",
                                            stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                                        max_treedepth = NULL, adapt_delta = .98),
                                            priors=make_priors(
                                             priors_to_change=list(
                                               "B_ingroup"=c(0.01, 3.5), 
                                               "B_outgroup"=c(0.01, 3.5),
                                               "focal_effects"=c(0, 10), 
                                               "target_effects"=c(0, 10), 
                                               "dyad_effects"=c(0, 10),
                                               "sr_sigma"=c(5), 
                                               "dr_sigma"=c(5) 
                                             )
                                           )
                                          )
res2 = summarize_strand_results(fit2)


# Prep results for plotting
tab1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), 
                                 normalized=TRUE,  only_slopes=TRUE, export_as_table=TRUE)
tab2 = strand_caterpillar_plot(res2, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), 
                                 normalized=TRUE,  only_slopes=TRUE, export_as_table=TRUE)

tab1$Priors = "Default, Weakly regularizing"
tab2$Priors = "Modified, Especially vague"

tab0 = rbind(tab1, tab2)

plot_0 = ggplot2::ggplot(tab0, ggplot2::aes(x = Variable, y = Median, 
        ymin = LI, ymax = HI, color=Priors)) + ggplot2::geom_linerange(size = 1,  position = ggplot2::position_dodge(0.33)) + 
        ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(0.33)) + ggplot2::facet_grid(Submodel ~ 
        ., scales = "free", space = "free") + ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + ggplot2::theme(legend.position="bottom") + ggplot2::scale_color_manual(values=c("#E69F00", "black"))


ggplot2::ggsave("Bat_slopes_priors.pdf", plot_0, width=7, height=2.5)

