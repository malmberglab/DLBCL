dfSelect <- function(df, var, values) {
  as.data.frame(do.call("rbind", lapply(values, function(x) { df[df[,var]==x,] })))
}



# getProteinValuesForSample()
# Returns values for selected proteins for one sample
#   sample          string          olinkID (same as row names in NPXData)
#   proteins        string vector   refers to proteinName in proteinData
#   ignoreBelowLOD  bool            values below LOD are replaced with NA

getProteinValuesForSample <- function(sample, proteins, ignoreBelowLOD=TRUE) {
  if (length(sample) > 1) stop("getProteinValuesForSample accepts only one sample")
  sapply(proteins, function(protein) {
    value <- NPXData[sample, protein]
    if (ignoreBelowLOD && value < proteinData[protein, paste0("LOD_", sampleData[sample, "LOD"])]) {
      return(NA)
    } else {
      return(value)
    }
  })
}


# fetch()
# Returns list with all data for selected samples
#   studyIDs        int vector      refers to studyID column in metaData
#   proteins        string vector   refers to proteinName in proteinData
#   study           any             refers to study column in sampleData
#   timepoint       any             refers to timepoint column in sampleData
#   removeExcluded  bool            samples marked as excluded in sampleData are removed
#   ignoreBelowLOD  bool            values below LOD will be replaced with NA


fetch <- function(studyIDs, proteins, study=2019, timepoint="pre", removeExcluded=TRUE, ignoreBelowLOD=TRUE) {
  filteredSampleData <- dfSelect(sampleData[sampleData$study==study & sampleData$timepoint==timepoint,],
                                 "studyID", studyIDs)
  if (removeExcluded) {
    excludedSampleData <- filteredSampleData[filteredSampleData$isExcluded,]
    filteredSampleData <- filteredSampleData[!filteredSampleData$isExcluded,]
  }
  
  filteredMetaData <- dfSelect(metaData, "studyID", filteredSampleData$studyID)
  filteredNPXData <- as.data.frame(t(sapply(filteredSampleData$olinkID, "getProteinValuesForSample",
                                            proteins=proteins, ignoreBelowLOD=ignoreBelowLOD)))
  
  list(
    NPXData = filteredNPXData,
    metaData = filteredMetaData,
    proteinData = proteinData[proteins,],
    sampleData = filteredSampleData,
    excludedSampleData = if (removeExcluded) excludedSampleData else NULL
  )
}


# normalizeNPX()

normalizeNPX <- function(fetchedData, normalization) {
  fetchedData$normalizedNPXData <- fetchedData$NPXData
  for (protein in names(normalization)) {
    fetchedData$normalizedNPXData[,protein] <- fetchedData$normalizedNPXData[,protein] + normalization[protein]
  }
  fetchedData$normalization <- normalization
  fetchedData
}


# subsetDataset()
# Subset a dataset based on studyIDs

subsetDataset <- function(dataset, studyIDs) {
  keep <- dataset$metaData$studyID %in% studyIDs
  message("Keeping ", sum(keep), " of ", length(keep), " samples.")
  dataset$NPXData <- dataset$NPXData[keep, ]
  dataset$metaData <- dataset$metaData[keep,]
  dataset$sampleData <- dataset$sampleData[keep,]
  if ("normalizedNPXData" %in% names(dataset)) dataset$normalizedNPXData <- dataset$normalizedNPXData[keep, ]
  if ("survival" %in% names(dataset)) dataset <- addSurvival(dataset)
  
  message("Data not subsetted: ", paste(names(dataset)[names(dataset) %in%
    c("NPXData", "metaData", "sampleData", "normalizedNPXData", "survival") == FALSE], collapse=", "))
  dataset
}




# addSurvival()
# Add survival data to a data set based on meta data

addSurvival <- function(dataset) {
  dataset$survival <- list(
    "OS"  = data.frame(
      days = dataset$metaData$days_treatment_to_last_date,
      status = sapply(dataset$metaData$dead, function(dead) { if (dead == "no") 0 else if (dead=="yes") 1 else NA })
    ),
    "PFS" = data.frame(
      days = dataset$metaData$days_treatment_to_last_control,
      status = ifelse(dataset$metaData$dead == "yes" | dataset$metaData$progression == "yes", 1, 0)
    )
  )

  dataset
}



useNormalizedDataIfAvailable <- function(dataset) {
  if ("normalizedNPXData" %in% names(dataset)) {
    message("Using normalized data")
    return(dataset$normalizedNPXData)
  } else {
    return(dataset$NPXData)
  }
}

