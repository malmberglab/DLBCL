
# Use the same protein order as in Figure 1A/B
proteinOrder <- list("healthy" = HD$NPXData, "patients" = osloMerged$NPXData) %>%
  compareTwoGroups(selectedProteins, pAdjustMethod="BH") %>%
  { .[order(-.$medianDiff), "protein"] }



#######################################################################################################################
# Figure S4A
#######################################################################################################################

subtypeComparisonOslo <- list(
  "Non-GCB DLBCL"    = osloMerged$NPXData[osloMerged$metaData$subtype5 == "non-GCB", ],
  "GCB DLBCL"        = osloMerged$NPXData[osloMerged$metaData$subtype5 == "GCB", ],
  "T/hist.-rich BCL" = osloMerged$NPXData[osloMerged$metaData$subtype5 == "T/hist", ],
  "HGBCL"            = osloMerged$NPXData[osloMerged$metaData$subtype5 == "HGBCL", ],
  "Other"            = osloMerged$NPXData[osloMerged$metaData$subtype5 == "other", ]
) %>%
  compareGroups(proteinOrder, pAdjustMethod = "BH",
                normalizeHeatmapToMedian = TRUE,
                boxplotColors = figureColors(c("lightBlue", "lightRed", "yellow", "orange", "cyan")))

subtypeComparisonOslo$statistics %>%
  write.table(file = file.path(outputDir, "statistics", "Figure S4A.csv"), sep = ";", dec = ",", row.names = FALSE)

subtypeComparisonOslo$boxplot %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S4A.pdf"), width = 8, height = 24, useDingbats = FALSE)




#######################################################################################################################
# Figure S4B
#######################################################################################################################

subtypeComparisonWU <- list(
  "Non-GCB DLBCL"    = WU$NPXData[WU$metaData$subtype5 == "non-GCB", ],
  "GCB DLBCL"        = WU$NPXData[WU$metaData$subtype5 == "GCB", ],
  "T/hist.-rich BCL" = WU$NPXData[WU$metaData$subtype5 == "T/hist", ],
  "Not known"        = WU$NPXData[WU$metaData$subtype5 == "", ]
) %>%
  compareGroups(proteinOrder, pAdjustMethod = "BH",
                normalizeHeatmapToMedian = TRUE,
                boxplotColors = figureColors(c("lightBlue", "lightRed", "yellow", "grey")))

subtypeComparisonWU$statistics %>%
  write.table(file = file.path(outputDir, "statistics", "Figure S4B.csv"), sep = ";", dec = ",", row.names = FALSE)

subtypeComparisonWU$boxplot %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S4B.pdf"), width = 8, height = 24, useDingbats = FALSE)



rm(proteinOrder, subtypeComparisonOslo, subtypeComparisonWU)



