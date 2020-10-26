
#######################################################################################################################
# Figure 4A
#######################################################################################################################

volcanoPlotData <- list(
  "low"  = osloMerged$NPXData[osloMerged$score$higherInPatientsQuantiles$groupedBySecondTertile == "low/mid",],
  "high" = osloMerged$NPXData[osloMerged$score$higherInPatientsQuantiles$groupedBySecondTertile == "high",]
) %>% 
  compareTwoGroups(proteinsForScore, pAdjustMethod = "BH")

volcanoPlot <-  ggplot() +
    geom_point(data = volcanoPlotData, mapping = aes(x = medianDiff, y = -log10(wilcoxAdj))) +
    geom_rect(data = NULL,
              mapping = aes(xmin = 0, xmax = 3.5, ymin = -log10(0.05)), ymax = 9.3, color = "red", fill = NA) +
    scale_y_continuous(breaks = (0:5)*2, labels = paste0("10^-", (0:5)*2)) + 
    theme_classic() +
    theme(axis.text = element_text(size = 16, color = "black"),
          aspect.ratio = 1,
          axis.ticks = element_line(color="black"))

volcanoPlotData %>%
  write.table(file = file.path(outputDir, "statistics", "Figure 4A.csv"), sep = ";", dec = ",", row.names = FALSE)

volcanoPlot %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 4A.pdf"), width = 4, height = 4, useDingbats = FALSE)

rm(volcanoPlotData, volcanoPlot)


#######################################################################################################################
# Figure 4B
#######################################################################################################################

scoreGroupsVsControlBoxPlots(osloMerged, "higherInPatientsQuantiles", "groupedBySecondTertile",
                             controlDataset = HD, proteins = scoreOptimizationProteinOrder,
                             colors=c("low/mid" = figureColors("lightBlue"), "high" = figureColors("lightRed"))) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 4B.pdf"), width = 18, height = 6)



#######################################################################################################################
# Figure 4C
#######################################################################################################################

scoreOptimization$plot %>% ggsave(file = file.path(outputDir, "figures", "Figure 4C.pdf"), width = 3, height = 2)



#######################################################################################################################
# Figure 4D
#######################################################################################################################

plotScoreSurvival(osloMerged, "OS", "lowestPossibleFromOslo", "groupedBySecondTertile",
                  colors = colorsForKaplanMeierPlots) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 4D Oslo.pdf"), width = 5.5, height = 8)

plotScoreSurvival(WU, "OS", "lowestPossibleFromOslo", "groupedBySecondTertile",
                  colors = colorsForKaplanMeierPlots) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 4D WU alt.pdf"), width = 5.5, height = 8)




