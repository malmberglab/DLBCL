
#######################################################################################################################
# Figure 3A
#######################################################################################################################

for (dataset in list(c("Oslo", "osloMerged"), c("St Louis", "WU"))) {
  histogramForScore(get(dataset[2]), "higherInPatientsQuantiles", "groupedByTertiles",
                    colors = c("low"  = figureColors("lightGreen"),
                               "mid"  = figureColors("lightBlue"),
                               "high" = figureColors("lightRed"))) %>%
    ggsave(file = file.path(outputDir, "figures", paste0("Figure 3A ", dataset[1], ".pdf")), width = 4.5, height = 2)
}



#######################################################################################################################
# Figure 3B
#######################################################################################################################

for (dataset in list(c("Oslo", "osloMerged"), c("St Louis", "WU"))) {
    plotScoreSurvival(get(dataset[2]), "OS", "higherInPatientsQuantiles", "groupedByTertiles",
                      colors = colorsForKaplanMeierPlots) %>%
      ggsave(file = file.path(outputDir, "figures", paste0("Figure 3B ", dataset[1], ".pdf")), width = 5.5, height = 8)
}



#######################################################################################################################
# Figure 3C
#######################################################################################################################

for (dataset in list(c("Oslo", "osloMerged"), c("St Louis", "WU"))) {

  clinicalData <- data.frame(
    scoreGroup = get(dataset[2])$score$higherInPatientsQuantiles$groupedByTertiles,
    gender     = get(dataset[2])$metaData$gender,
    age        = ifelse(get(dataset[2])$metaData$age > 60, "above 60", "up to 60"),
    subtype    = get(dataset[2])$metaData$subtype5,
    stage      = get(dataset[2])$metaData$stage,
    Bsymptoms  = get(dataset[2])$metaData$stage_AB,
    IPI        = ifelse(!is.na(get(dataset[2])$metaData$IPI), paste0("IPI ", get(dataset[2])$metaData$IPI), NA),
    RIPIgroup  = sapply(get(dataset[2])$metaData$IPI,
                        function(x) { if (is.na(x)) return(NA); if (x > 2) "poor" else "very good/good" })
  )
  
  clinicalData[clinicalData$Bsymptoms == "", "Bsymptoms"] <- "unknown"
  clinicalData[clinicalData$subtype == "", "subtype"] <- "unknown"
  
  clinicalPlotColorsForOptions <- c(
    "female" = "blue",
    "male" = "red",
    "up to 60" = "blue",
    "above 60" = "red",
    "GCB" = "blue",
    "non-GCB" = "red",
    "T/hist" = "yellow",
    "HGBCL" = "orange",
    "other" = "cyan",
    "I" = "green",
    "II" = "yellow",
    "III" = "orange",
    "IV" = "red",
    "A" = "blue",
    "B" = "red",
    "no X" = "blue",
    "X" = "red",
    "no E" = "blue",
    "E" = "red",
    "IPI 0" = "green",
    "IPI 1" = "greenYellow",
    "IPI 2" = "yellow",
    "IPI 3" = "lightOrange",
    "IPI 4" = "orange",
    "IPI 5" = "red",
    "very good/good" = "blue",
    "poor" = "red",
    "unknown" = "grey"
  )
  
  clinicalDataComparison <- performClinicalDataComparison(clinicalData, clinicalPlotColorsForOptions)
  
  ggsave(clinicalDataComparison$plot,
         file = file.path(outputDir, "figures", paste0("Figure 3C ", dataset[1], ".pdf")), width = 12, height = 2.2)
  
  write.table(clinicalData,
              file = file.path(outputDir, "statistics", paste0("Figure 3C ", dataset[1], " data.csv")),
              sep = ";", dec = ",", row.names = FALSE)
  
  write.table(clinicalDataComparison$statistics,
              file = file.path(outputDir, "statistics", paste0("Figure 3C ", dataset[1], " statistics.csv")),
              sep = ";", dec = ",", row.names = FALSE)
  
  rm(clinicalData, clinicalPlotColorsForOptions, clinicalDataComparison)
}



#######################################################################################################################
# Figure 3D
#######################################################################################################################

for (dataset in list(c("Oslo", "osloMerged"), c("St Louis", "WU"))) {
  plotScoreSurvival(get(dataset[2]), "OS", "higherInPatientsQuantiles", "groupedByTertiles",
                           query = get(dataset[2])$metaData$IPI %in% 3:5,
                           colors = colorsForKaplanMeierPlots) %>%
    ggsave(file = file.path(outputDir, "figures", paste0("Figure 3D ", dataset[1], ".pdf")), width = 5.5, height = 8)
}

