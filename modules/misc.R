
#######################################################################################################################
# Follow-up time
#######################################################################################################################

allOsloStudyIDs <- unique(c(osloMerged$metaData$studyID, "4", "48", "90", "79", "119", "131", "135", "164", "168",
                            "214", "223", "230", "270", "278", "287", "288", "354", "373", "404", "432", "446", "494",
                            "101", "137", "146", "240", "377", "419", "444"))

osloAll <- list(
  "metaData" = metaData[metaData$studyID %in% allOsloStudyIDs, ]
) %>% addSurvival()


calculateFollowUpTime <- function(dataset, name = deparse(substitute(dataset))) {
  df <- data.frame(
    OS = dataset$survival$OS$days,
    status = ifelse(dataset$survival$OS$status == 1, 0, 1),
    group = 1
  )
  cat("\n===  Follow-up time for dataset: ", name, "  ===\n", sep="")
  cat("Median OS:             ", median(df$OS), " days\n", sep="")
  cat("Reverse Kaplan-Meier:  ", summary(survfit(Surv(OS, status) ~ group, data = df))$table["median"],
      " days\n\n", sep = "")
}

sink(file.path(outputDir, "statistics", "Follow-up time.txt"))
calculateFollowUpTime(osloAll)
calculateFollowUpTime(osloMerged)
calculateFollowUpTime(osloOne)
calculateFollowUpTime(WU)
sink()



#######################################################################################################################
# Time from EOT to sample 2
#######################################################################################################################

osloWithPost <- list(
  metaData = osloOne$metaData[osloOne$metaData$studyID %in% intersect(osloOnePost$metaData$studyID,
                                                                      osloOne$metaData$studyID),]
)

sink(file.path(outputDir, "statistics", "Time from EOT.txt"))
print(summary(osloWithPost$metaData$days_EOT_to_post_sample))
sink()



rm(allOsloStudyIDs, osloAll, calculateFollowUpTime)