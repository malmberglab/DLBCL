
#######################################################################################################################
# Table S2
#######################################################################################################################

healthyPatientComparison <- list("healthy" = HD$NPXData, "patients" = osloMerged$NPXData) %>%
  compareTwoGroups(selectedProteins, pAdjustMethod = "BH")
healthyPatientComparison <- healthyPatientComparison[order(-healthyPatientComparison$medianDiff),]

data.frame(
  protein = healthyPatientComparison$protein,
  proteinName = healthyPatientComparison$proteinName,
  higherInPatients = ifelse(healthyPatientComparison$medianDiff > 0 & healthyPatientComparison$wilcoxAdj <= 0.05, ".", ""),
  lowerInPatients = ifelse(healthyPatientComparison$medianDiff < 0 & healthyPatientComparison$wilcoxAdj <= 0.05, ".", "")
) %>% write.table(file = file.path(outputDir, "tables", "Table S2.csv"), sep = ";", dec = ",", row.names = FALSE)

rm(healthyPatientComparison)