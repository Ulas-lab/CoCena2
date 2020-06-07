cutoff_stats_concise = cutoff_stats %>% 
  dplyr::select(R.squared, cutoff, no_edges, no_nodes, no_of_networks) %>% 
  dplyr::distinct()

cutoff_stats_concise = cutoff_stats_concise %>% filter(no_of_networks!=0)
rownames(cutoff_stats_concise) = cutoff_stats_concise$cutoff
cutoff_stats_concise = cutoff_stats_concise[,-2]

#print(cutoff_stats_concise)


crit_minmax = c("max","max", "max", "min" )
names(crit_minmax) = colnames(cutoff_stats_concise)



#plotRadarPerformanceTable(cutoff_stats_concise, crit_minmax,
#                          overlay=FALSE, bw=TRUE, lwd =5)

normalizationTypes <- rep("percentageOfMax", ncol(cutoff_stats_concise))
names(normalizationTypes) = colnames(cutoff_stats_concise)

nPT = normalizePerformanceTable(cutoff_stats_concise[,c("R.squared", "no_edges", "no_nodes", "no_of_networks")], normalizationTypes)
w = c(0.5,0.1,0.5, -1)
names(w) <- colnames(nPT)
ws<-weightedSum(nPT,w)
ranked_ws <- rank(-ws) %>% sort()


calculated_optimal_cutoff <- as.numeric(names(ranked_ws[1]))

print(paste0("The calculated optimal cutoff is: ", calculated_optimal_cutoff))


stats_calculated_optimal_cutoff <- cutoff_stats[cutoff_stats$cutoff == calculated_optimal_cutoff, c("degree", "Probs")]

dd_plot_calculated_optimal = ggplot(stats_calculated_optimal_cutoff,aes(x=log(degree), y= log(Probs))) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw() + 
  ggtitle(paste0("Calculated optimal correlation cut-off [",calculated_optimal_cutoff, "]"))

# print(dd_plot_calculated_optimal)
ggsave(paste0("Degree_distribution_plot_", calculated_optimal_cutoff, "_calculated_optimal_cutoff.pdf"), 
       dd_plot_calculated_optimal, device = cairo_pdf, width = 6, height = 5)

return(dd_plot_calculated_optimal)

# if (manual_cutoff) {
#   
#   optimal_cutoff <- as.numeric(readline(prompt = "Enter the cut-off to be used for the analysis: "))
#   
# } else {
#   
#   
#   optimal_cutoff = calculated_optimal_cutoff
#   
# }




# Other option:

# criteriaMinMax = c("max","max", "max", "min" )
# alternativesRanks = rank(-ws)
# criteriaNumberOfBreakPoints <- c(3,3,2)
# pT = cutoff_stats_concise[,c("R.squared", "no_nodes", "no_of_networks")]
# names(criteriaNumberOfBreakPoints) <- colnames(pT)
# criteriaLBs=apply(pT,2,min)
# names(criteriaLBs) <- colnames(pT)
# 
# criteriaUBs=apply(pT,2,max)
# names(criteriaUBs) <- colnames(pT)
# 
# epsilon <- 0.01
# 
# x<-UTA(pT, criteriaMinMax,
#        criteriaNumberOfBreakPoints, epsilon,
#        alternativesRanks = alternativesRanks,
#        criteriaLBs = criteriaLBs, criteriaUBs = criteriaUBs)
# 
# 
# tPT <- applyPiecewiseLinearValueFunctionsOnPerformanceTable(
#   x$valueFunctions,
#   pT
# )
# mavt <- weightedSum(tPT,rep(1,3))
# 
# 
# plotAlternativesValuesPreorder(mavt, decreasing=TRUE)



##########################################################################################


# Thomas cutoff:

# num_net_vector=cutoff_stats$no_of_networks %>%
#   unique() %>%
#   sort()
# n= num_net_vector %>% length()
# num_subnetworks = num_net_vector[2] #@thomas : why do we do this?
# 
# optimal_cutoff = cutoff_stats %>%
#   dplyr::filter(no_of_networks==num_subnetworks)%>%
#   dplyr::filter(no_edges > as.numeric(quantile(no_edges)["25%"])) %>% # filtering 1st quarter to 3rd quarter
#   dplyr::filter(no_edges < as.numeric(quantile(no_edges)["75%"])) %>%
#   dplyr::filter(no_nodes >0.3*max(no_nodes)) %>%
#   dplyr::arrange(desc(R.squared)) %>%
#   dplyr::filter(dplyr::row_number()==1) %>%
#   subset(select=cutoff) %>%
#   as.numeric()

