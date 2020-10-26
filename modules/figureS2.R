
# Use the same protein order as in Figure 1A/B
proteinOrder <- list("healthy" = HD$NPXData, "patients" = osloMerged$NPXData) %>%
  compareTwoGroups(selectedProteins, pAdjustMethod = "BH") %>%
  { .[order(-.$medianDiff), "protein"] }



#######################################################################################################################
# Figure S2A
#######################################################################################################################

sexComparisonHD <- list("male"   = HD$NPXData[HD$metaData$gender == "male", ],
                        "female" = HD$NPXData[HD$metaData$gender == "female", ]) %>%
  compareGroups(proteinOrder, pAdjustMethod = "BH", normalizeHeatmapToMedian = TRUE,
                boxplotColors = figureColors(c("lightRed", "lightBlue")))

sexComparisonHD$statistics %>%
  write.table(file = file.path(outputDir, "statistics", "Figure S2A.csv"), sep = ";", dec = ",", row.names = FALSE)
sexComparisonHD$boxplot %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S2A.pdf"), width = 8, height = 16, useDingbats = FALSE)



#######################################################################################################################
# Figure S2B
#######################################################################################################################

ageComparisonHD <- list("Above 60" = HD$NPXData[HD$metaData$age > 60, ],
                        "Up to 60" = HD$NPXData[HD$metaData$age <= 60, ]) %>%
  compareGroups(proteinOrder, pAdjustMethod = "BH", normalizeHeatmapToMedian = TRUE,
                boxplotColors = figureColors(c("lightRed", "lightBlue")))

ageComparisonHD$statistics %>%
  write.table(file = file.path(outputDir, "statistics", "Figure S2B.csv"), sep = ";", dec = ",", row.names = FALSE)

ageComparisonHD$boxplot %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S2B.pdf"), width = 8, height = 16, useDingbats = FALSE)



rm(proteinOrder, sexComparisonHD, ageComparisonHD)