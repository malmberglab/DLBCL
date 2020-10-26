
# Use the same protein order as in Figure 1A/B
proteinOrder <- list("healthy" = HD$NPXData, "patients" = osloMerged$NPXData) %>%
  compareTwoGroups(selectedProteins, pAdjustMethod="BH") %>%
  { .[order(-.$medianDiff), "protein"] }



#######################################################################################################################
# Figure S3A
#######################################################################################################################

sexComparisonOslo <- list(
  "female" = osloMerged$NPXData[osloMerged$metaData$gender == "female", ],
  "male"   = osloMerged$NPXData[osloMerged$metaData$gender == "male", ]
) %>%
  compareGroups(proteinOrder, pAdjustMethod = "BH",
                normalizeHeatmapToMedian = TRUE,
                boxplotColors = figureColors(c("lightBlue", "lightRed")))

sexComparisonOslo$statistics %>%
  write.table(file = file.path(outputDir, "statistics", "Figure S3A.csv"), sep = ";", dec = ",", row.names = FALSE)

sexComparisonOslo$boxplot %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S3A.pdf"), width = 8, height = 16, useDingbats = FALSE)




#######################################################################################################################
# Figure S3B
#######################################################################################################################

sexComparisonWU <- list(
  "female" = WU$NPXData[WU$metaData$gender == "female", ],
  "male"   = WU$NPXData[WU$metaData$gender == "male", ]
) %>%
  compareGroups(proteinOrder, pAdjustMethod = "BH",
                normalizeHeatmapToMedian = TRUE,
                boxplotColors = figureColors(c("lightBlue", "lightRed")))

sexComparisonWU$statistics %>%
  write.table(file = file.path(outputDir, "statistics", "Figure S3B.csv"), sep = ";", dec = ",", row.names = FALSE)

sexComparisonWU$boxplot %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S3B.pdf"), width = 8, height = 16, useDingbats = FALSE)



rm(proteinOrder, sexComparisonOslo, sexComparisonWU)
