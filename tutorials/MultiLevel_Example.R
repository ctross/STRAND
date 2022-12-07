
 library(STRAND)
 library(igraph)
 library(ggplot2)
 library(xtable)

 set.seed(1)
 V = 1            # One blocking variable
 G = 3            # Three categories in this variable
 N_id = 100       # Number of people
 N_hh = 33        # Number of households

 B = matrix(-19, nrow=G, ncol=G)
 diag(B) = -10.5 # Block matrix

 B[1,3] = -16.9
 B[3,2] = -13.9

 Clan = sample(c("Darmok", "Jalad", "Tanagra"), N_id, replace=TRUE)
 HH = sample(1:N_hh, size=N_id, replace=TRUE)

 dd = data.frame(ID=c(1:N_id),
                 Income=rnorm(N_id, 0,1), 
                 Emeralds=rnorm(N_id, 0,1),
                 Gold=rbinom(N_id, size=15, prob=0.05), 
                 Kindness=rnorm(N_id,0,1),
                 Freckles=rnorm(N_id,0,1),
                 Clan=Clan,
                 HH=HH
                 )

IncomeDifference = KindnessDifference = matrix(NA, nrow=N_id, ncol=N_id)
  for(i in 1:N_id){
    for(j in 1:N_id){
        IncomeDifference[i,j] = abs(dd$Income[i] - dd$Income[j])
        KindnessDifference[i,j] = dd$Kindness[i] - dd$Kindness[j]
    }
 }

 DyadicPreds = array(NA, c(N_id,N_id,2))
 DyadicPreds[,,1] = IncomeDifference
 DyadicPreds[,,2] = KindnessDifference

 hh = data.frame(ID=c(1:N_hh),
                Land=rnorm(N_hh, 0,1), 
                Trees=rbinom(N_hh, size=10, prob=0.05), 
                Wealth=rnorm(N_hh,0,1),
                Location=rnorm(N_hh,0,1)
                )

Distance = WealthDifference = matrix(NA, nrow=N_hh, ncol=N_hh)

 for(i in 1:N_hh){
    for(j in 1:N_hh){
        Distance[i,j] = abs(hh$Location[i] - hh$Location[j])
        WealthDifference[i,j] = hh$Wealth[i] - hh$Wealth[j]
    }
 }

 HHDyadicPreds = array(NA, c(N_hh,N_hh,2))
 HHDyadicPreds[,,1] = Distance
 HHDyadicPreds[,,2] = WealthDifference
 
 A = simulate_sbm_plus_srm_hh_network(
                          N_id = N_id, 
                          N_hh=N_hh, 
                          HH=HH, 
                          B=list(B=B), 
                          V=V, 
                          groups=data.frame(Clan=factor(Clan)),
                          individual_predictor=data.frame(dd$Emeralds, dd$Gold, dd$Freckles), 
                          individual_effects=matrix(rbind(c(1.7, -2.3), c(-1.1, 0.99), c(1.2, 2.4)),ncol=3, nrow=2, byrow=TRUE),
                          dyadic_predictor=DyadicPreds, 
                          dyadic_effects=c(1.3, -0.8),
                          hh_individual_predictor=data.frame(hh$Land, hh$Trees), 
                          hh_individual_effects=matrix(rbind(c(1.2, -1.3),c(2.5,2.5)),ncol=2, nrow=2, byrow=TRUE),
                          hh_dyadic_predictor=HHDyadicPreds, 
                          hh_dyadic_effects=c(-2.3, 1.8),
                          sr_sigma = c(1.9, 0.8), 
                          sr_rho = -0.5,
                          dr_sigma = 3.2, 
                          dr_rho = 0.7,
                          hh_sr_sigma = c(1.2, 1.8), 
                          hh_sr_rho = 0.5,
                          hh_dr_sigma = 2.9, 
                          hh_dr_rho = 0.7,
                          mode="bernoulli"
                                )

 Net = graph_from_adjacency_matrix(A$network, mode = c("directed"))
 V(Net)$color = c("turquoise4","gray13", "goldenrod3")[A$group_ids$Clan]

 bob = ggnet2(Net, size = 6, color = "color",layout.par = list(repulse.rad = 10, area = 300))
 bob
 ggsave("SampleNet.PDF", bob, width=8.5, height=8.5)

 dat = make_strand_data(self_report = list(Outcome=A$network),
                        hh_ID = HH,
                        block_covariates = data.frame(Clan=factor(Clan)), 
                        individual_covariates = data.frame(Emeralds = dd$Emeralds, Gold = dd$Gold, Freckles = dd$Freckles), 
                        dyadic_covariates = list(IncomeDifference = IncomeDifference, KindnessDifference = KindnessDifference),
                        hh_individual_covariates = data.frame(Land = hh$Land, Trees = hh$Trees), 
                        hh_dyadic_covariates = list(Distance = Distance, WealthDifference = WealthDifference),
                        outcome_mode="bernoulli"
                        )

 fit = fit_block_plus_social_relations_hh_model(data=dat,
                                               block_regression = ~ Clan,
                                               focal_regression = ~ Emeralds + Gold + Freckles,
                                               target_regression = ~ Emeralds + Gold + Freckles,
                                               dyad_regression = ~ IncomeDifference + KindnessDifference,
                                               hh_focal_regression = ~ Land + Trees,
                                               hh_target_regression = ~ Land + Trees,
                                               hh_dyad_regression = ~ Distance + WealthDifference,
                                               mode="mcmc",
                                               stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                                          max_treedepth = NULL, adapt_delta = .98)
)


res = summarize_strand_results(fit)

df = res$summary
df$Type = c("SD", "Slope","Slope", "Slope","SD", "Slope","Slope", "Slope","SD", "Slope", "Slope","Rho", 
            "Rho","Offset", "Offset", "Offset", "Offset","Offset", "Offset", "Offset", "Offset","Offset", 
            "Offset", "SD","Slope", "Slope","SD", "Slope","Slope", "SD","Slope", "Slope","Rho","Rho")

df$Category = c("Deviation", "Sender","Sender", "Sender","Deviation", "Receiver","Receiver", "Receiver","Deviation", "Dyadic", "Dyadic","Correlation", 
            "Correlation","Block", "Block", "Block", "Block","Block", "Block", "Block", "Block","Block", 
            "Block", "Deviation","Sender", "Sender","Deviation", "Receiver","Receiver", "Deviation","Dyadic", "Dyadic","Correlation","Correlation")

df$Category = factor(df$Category)
df$Category = factor(df$Category, levels=c("Block", "Sender", "Receiver", "Dyadic", "Correlation", "Deviation"))

colnames(df) = c("Variable",  "Median",  "L",  "H",    "Mean",    "SD", "Type","Category")

df$Median = as.numeric(as.character(df$Median))
df$L = as.numeric(as.character(df$L))
df$H = as.numeric(as.character(df$H))

# To test model fit, check that generative parameters are recovered
Offset = -7.5   # The data had no intercept parameter, but the model does, so we have to shift these data points
Shrinkage = 0.7 # Shrinkage priors scale everything twoards zero, so scale the data a bit to see comparable points
df$True = c(1.9, 1.7, -1.1, 1.2, 0.8, -2.3, 0.99, 2.4, 3.2, 1.3, -0.8, -0.5, 0.7, NA, B[1,1]-Offset, B[1,2]-Offset, B[1,3]-Offset, 
B[2,1]-Offset, B[2,2]-Offset, B[2,3]-Offset, B[3,1]-Offset, B[3,2]-Offset, B[3,3]-Offset, 1.2, 1.2,  2.5, 1.8, -1.3, 2.5, 2.9, -2.3, 1.8, 0.5,0.7)*Shrinkage

testfit = ggplot(df,aes(x=Variable,y=Median,ymin=L,ymax=H, group=Type, color=Type))+ 
     geom_linerange(size=1.25, position = position_dodge(width = 0.5))+
     geom_point(size=2, position = position_dodge(width = 0.5))+
     geom_point(aes(x=Variable,y=True),size=2, position = position_dodge(width = 0.5), color="darkred")+
     geom_hline(aes(yintercept=0),color="black",linetype="dashed")+coord_flip() + facet_grid(rows = vars(Category),scale="free", space="free") +
     labs(x="Variable", y="Value") + theme(strip.text.x = element_text(size=12,face="bold"), 
     strip.text.y = element_text(size=12,face="bold"),axis.text=element_text(size=12),axis.title=element_text(size=14,
     face="bold"))+theme(strip.text.y = element_text(angle = 360)) +  theme(panel.spacing = unit(1, "lines")) +
      scale_color_manual(
      values = c("#E28B3E","#5DA2A0","#DF6866","#195A41")
     ) + theme(legend.position="none") + theme(legend.text=element_text(size=12)) 
testfit

ggsave("TestFit_HH_SRM.pdf",testfit, width=10, height=8)


print(xtable(df[order(df$Category),c(8,1,2,3,4)]),include.rownames=FALSE)

