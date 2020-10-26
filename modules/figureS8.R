
#######################################################################################################################
# Figure S8A
#######################################################################################################################

plotScoreSurvival(osloMerged, "OS", "higherInPatientsQuantiles", "groupedByTertiles",
                  query = osloMerged$metaData$subtype5 %in% c("GCB", "non-GCB"),
                  colors = colorsForKaplanMeierPlots) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S8A.pdf"), width = 5.5, height = 8, useDingbats = FALSE)



#######################################################################################################################
# Figure S8B
#######################################################################################################################

plotScoreSurvival(WU, "OS", "higherInPatientsQuantiles", "groupedByTertiles",
                  query = WU$metaData$subtype5 %in% c("GCB", "non-GCB"),
                  colors = colorsForKaplanMeierPlots) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S8B.pdf"), width = 5.5, height = 8, useDingbats = FALSE)



#######################################################################################################################
# Figure S8C
#######################################################################################################################

plotScoreSurvival(osloMerged, "OS", "higherInPatientsQuantiles", "groupedByTertiles",
                  query = osloMerged$metaData$subtype5 == "GCB",
                  colors = colorsForKaplanMeierPlots) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S8C.pdf"), width = 5.5, height = 8, useDingbats = FALSE)



#######################################################################################################################
# Figure S8D
#######################################################################################################################

plotScoreSurvival(WU, "OS", "higherInPatientsQuantiles", "groupedByTertiles",
                  query = WU$metaData$subtype5 == "GCB",
                  colors = colorsForKaplanMeierPlots) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S8D.pdf"), width = 5.5, height = 8, useDingbats = FALSE)

