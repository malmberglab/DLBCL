
#######################################################################################################################
# Table 1
#######################################################################################################################

osloWithPost <- list(
  metaData = osloOne$metaData[osloOne$metaData$studyID %in% intersect(osloOnePost$metaData$studyID,
                                                                      osloOne$metaData$studyID),]
)

cytofStudyIDs <- read.csv(file.path(dataDir, "fcs.csv"), sep = ";", dec = ",") %>%
  { .[.$type == "patient", "studyID"] }
osloCyTOF <- list(
  metaData = metaData[metaData$type == "patient" & metaData$studyID %in% cytofStudyIDs, ]
)

datasets <- c("HD", "osloMerged", "osloWithPost", "osloCyTOF", "WU")

table1 <- do.call(data.frame, sapply(datasets, function(dataset) {
  m <- get(dataset)$metaData
  
  sumWithPercent <- function(x) {
    paste0(sum(x), " (", round(sum(x)/nrow(m)*100, 1), "%)")
  }
  
  c(
    "Total" = nrow(m),
    
    "Gender" = "",
    "Female" = sumWithPercent(m$gender=="female"),
    "Male"   = sumWithPercent(m$gender=="male"),

    "Age (median, range)" = paste0(median(m$age), " (", min(m$age), "-", max(m$age), ")"),

    "Subtype"         = "",
    "GCB"             = sumWithPercent(m$subtype5 == "GCB"),
    "non-GCB"         = sumWithPercent(m$subtype5 == "non-GCB"),
    "High-grade BCL"  = sumWithPercent(m$subtype5 == "HGBCL"),
    "T/hist rich"     = sumWithPercent(m$subtype5 == "T/hist"),
    "Other subtype"   = sumWithPercent(m$subtype5 == "other"),
    "Unknown subtype" = sumWithPercent(m$subtype5 %in% c("GCB", "non-GCB", "HGBCL", "T/hist", "other") == FALSE),

    "Stage" = "",
    "I"   = sumWithPercent(m$stage == "I"),
    "II"  = sumWithPercent(m$stage == "II"),
    "III" = sumWithPercent(m$stage == "III"),
    "IV"  = sumWithPercent(m$stage == "IV"),

    "IPI detailed" = "",
    "IPI 0" = sumWithPercent(m$IPI %in% 0),
    "IPI 1" = sumWithPercent(m$IPI %in% 1),
    "IPI 2" = sumWithPercent(m$IPI %in% 2),
    "IPI 3" = sumWithPercent(m$IPI %in% 3),
    "IPI 4" = sumWithPercent(m$IPI %in% 4),
    "IPI 5" = sumWithPercent(m$IPI %in% 5),

    "First-line treatment" = "",
    "R-CHOP" = sumWithPercent(m$protocol == "R-CHOP"),
    "R-EPOCH/R-CHOEP" = sumWithPercent(m$protocol == "R-CHOEP/R-EPOCH"),
    "Other treatment" = sumWithPercent(m$protocol == "other"),
    "Unknown treatment" = sumWithPercent(m$protocol %in% c("R-CHOP", "R-CHOEP/R-EPOCH", "other") == FALSE)
  )
}, simplify=FALSE, USE.NAMES=TRUE))

table1 %>% write.table(file = file.path(outputDir, "tables", "Table 1.csv"),
                       sep = ";", dec = ",", row.names = TRUE, col.names = NA)



rm(cytofStudyIDs, osloWithPost, osloCyTOF, datasets, table1)