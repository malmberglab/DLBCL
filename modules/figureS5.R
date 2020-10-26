
#######################################################################################################################
# Figure S5A:
#######################################################################################################################

list(
  "20+" = WU$metaData$age_group == "20-29",
  "30+" = WU$metaData$age_group == "30-39",
  "40+" = WU$metaData$age_group == "40-49",
  "50+" = WU$metaData$age_group == "50-59",
  "60+" = WU$metaData$age_group == "60-69",
  "70+" = WU$metaData$age_group == "70-79",
  "80+" = WU$metaData$age_group == "80-89",
  "90+" = WU$metaData$age_group == "90-99"
) %>% { metaDataMultipleGroupComparison(WU, ., selectedProteins, customPlotLimits=c(-5, 5),
                                        statisticsDataset = osloMerged, statisticsGroups = list(
                                          "30+" = osloMerged$metaData$age_group == "30-39",
                                          "40+" = osloMerged$metaData$age_group == "40-49",
                                          "50+" = osloMerged$metaData$age_group == "50-59",
                                          "60+" = osloMerged$metaData$age_group == "60-69",
                                          "70+" = osloMerged$metaData$age_group == "70-79",
                                          "80+" = osloMerged$metaData$age_group == "80-89"
                                        )) } %>% {
  ggsave(.$plot, file = file.path(outputDir, "figures", "Figure S5A.pdf"),
         width = 4.5, height = 20, useDingbats = FALSE)
  write.table(.$statistics, file = file.path(outputDir, "statistics", "Figure S5A.csv"),
              sep = ";", dec = ",", row.names = FALSE)
}



#######################################################################################################################
# Figure S5B:
#######################################################################################################################

list(
  "I"   = WU$metaData$stage == "I",
  "II"  = WU$metaData$stage == "II",
  "III" = WU$metaData$stage == "III",
  "IV"  = WU$metaData$stage == "IV"
) %>% { metaDataMultipleGroupComparison(WU, ., selectedProteins, customPlotLimits = c(-5, 5),
                                        statisticsDataset = osloMerged, statisticsGroups = list(
                                          "I" = osloMerged$metaData$stage == "I",
                                          "II" = osloMerged$metaData$stage == "II",
                                          "III" = osloMerged$metaData$stage == "III",
                                          "IV" = osloMerged$metaData$stage == "IV"
                                        )) } %>% {
  ggsave(.$plot, file = file.path(outputDir, "figures", "Figure S5B.pdf"),
         width = 4.5, height = 20, useDingbats = FALSE)
  write.table(.$statistics, file = file.path(outputDir, "statistics", "Figure S5B.csv"),
              sep = ";", dec = ",", row.names = FALSE)
}



#######################################################################################################################
# Figure S5C:
#######################################################################################################################

list(
  "0" = !is.na(WU$metaData$IPI) & WU$metaData$IPI == 0,
  "1" = !is.na(WU$metaData$IPI) & WU$metaData$IPI == 1,
  "2" = !is.na(WU$metaData$IPI) & WU$metaData$IPI == 2,
  "3" = !is.na(WU$metaData$IPI) & WU$metaData$IPI == 3,
  "4" = !is.na(WU$metaData$IPI) & WU$metaData$IPI == 4,
  "5" = !is.na(WU$metaData$IPI) & WU$metaData$IPI == 5
) %>% { metaDataMultipleGroupComparison(WU, ., selectedProteins, customPlotLimits = c(-5, 5),
                                        statisticsDataset = osloMerged, statisticsGroups = list(
                                          "0" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 0,
                                          "1" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 1,
                                          "2" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 2,
                                          "3" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 3,
                                          "4" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 4,
                                          "5" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 5
                                        )) } %>% {
  ggsave(.$plot, file = file.path(outputDir, "figures", "Figure S5C.pdf"),
         width = 4.5, height = 20, useDingbats = FALSE)
  write.table(.$statistics, file = file.path(outputDir, "statistics", "Figure S5C.csv"),
              sep = ";", dec = ",", row.names = FALSE)
}


