# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gplots)
library(grid)
library(gridExtra)
library(survminer)
library(survival)
library(survivalsvm)
library(stringr)
library(flowCore)
library(Rtsne)
library(MASS)
library(reshape2)
library(xml2)


# Options
options(stringsAsFactors = FALSE)


# Load functions
for (file in list.files(file.path(dir, "functions"), full = TRUE)) {
  source(file)
}


# Load data
NPXData <- read.csv(file.path(dataDir, "data.csv"), sep=";", dec=",", stringsAsFactors=FALSE, row.names=1)
proteinData <- read.csv(file.path(dataDir, "proteins.csv"), sep=";", dec=",", stringsAsFactors=FALSE, row.names=1)
sampleData <- read.csv(file.path(dataDir, "samples.csv"), sep=";", dec=",", stringsAsFactors=FALSE, row.names=1)
olinkIDData <- read.csv(file.path(dataDir, "olinkIDs.csv"), sep=";", dec=",", stringsAsFactors=FALSE, row.names=NULL)
metaData <- read.csv(file.path(dataDir, "meta.csv"), sep=";", dec=",", stringsAsFactors=FALSE, row.names=NULL)
metaData <- metaData[metaData$excluded == "",]

# Attach studyIDs and timepoints to sampleData
sampleData$studyID <- sapply(row.names(sampleData),
                             function(olinkID) { olinkIDData[olinkIDData$olinkID == olinkID, "studyID"] })
sampleData$timepoint <- sapply(row.names(sampleData),
                               function(olinkID) { olinkIDData[olinkIDData$olinkID == olinkID, "timepoint"] })
rm(olinkIDData)
sampleData$olinkID <- row.names(sampleData)

# Exclude duplicates and samples with failed QC
sampleData$isExcluded <- sampleData$QC != "Pass"
sampleData$isExcluded <- sampleData$isExcluded | endsWith(row.names(sampleData), "_P2")

# Select proteins
selectedProteins <- row.names(proteinData)[proteinData$missingFreq_20190358_plasma < 0.5 &
                                           proteinData$missingFreq_20190358_serum  < 0.5 &
                                           proteinData$missingFreq_20170447_serum  < 0.5]



### FETCH, NORMALIZE AND MERGE DATA ###

# Fetch Oslo I data (including healthy donors)
oslo2017 <- fetch(subset(metaData, origin %in% c("HD", "Oslo I"))$studyID,
                  proteins = selectedProteins, study = 2017, timepoint = "pre", ignoreBelowLOD = FALSE)

# Fetch Oslo I+II data from 2019 study (including healthy donors)
oslo2019 <- fetch(subset(metaData, origin %in% c("HD", "Oslo I", "Oslo II"))$studyID,
                  proteins = selectedProteins, study = 2019, timepoint = "pre", ignoreBelowLOD = FALSE)

# Find bridging samples
bridgingSamples <- intersect(oslo2017$metaData$studyID, oslo2019$metaData$studyID)
oslo2017bridging <- fetch(bridgingSamples,
                          proteins = selectedProteins, study = 2017, timepoint = "pre", ignoreBelowLOD = FALSE)
oslo2019bridging <- fetch(bridgingSamples,
                          proteins = selectedProteins, study = 2019, timepoint = "pre", ignoreBelowLOD = FALSE)

# Calculate normalization values
oslo2017normalization <- sapply(selectedProteins, function(protein) {
  median(sapply(bridgingSamples, function(sample) {
    oslo2019bridging$NPXData[oslo2019bridging$metaData$studyID == sample, protein] -
      oslo2017bridging$NPXData[oslo2017bridging$metaData$studyID == sample, protein]
  }))
})

rm(oslo2017, oslo2019, bridgingSamples, oslo2017bridging, oslo2019bridging)


# Fetch Oslo I data
osloOne <- fetch(subset(metaData, origin == "Oslo I")$studyID,
                 proteins = selectedProteins, study = 2017, timepoint = "pre", ignoreBelowLOD = FALSE)

# Fetch Oslo II data
osloTwo <- fetch(subset(metaData, origin=="Oslo II")$studyID,
                 proteins = selectedProteins, study = 2019, timepoint = "pre", ignoreBelowLOD = FALSE)

# Fetch HD data
HD <- fetch(subset(metaData, origin=="HD")$studyID,
            proteins = selectedProteins, study = 2019, timepoint = "pre", ignoreBelowLOD = FALSE)

# Fetch Oslo I post data (432 excluded because post sample is taken before end of treatment)
osloOnePost <- fetch(setdiff(subset(metaData, origin == "Oslo I")$studyID, "432"),
                     proteins = selectedProteins, study = 2017, timepoint = "post", ignoreBelowLOD = FALSE)


# Normalize
osloOne <- osloOne %>% normalizeNPX(oslo2017normalization)
osloOnePost <- osloOnePost %>% normalizeNPX(oslo2017normalization)

rm(oslo2017normalization)


# Add survival
osloOne <- osloOne %>% addSurvival()
osloTwo <- osloTwo %>% addSurvival()


# Merge Oslo data
osloMerged <- list(
  NPXData = rbind(osloOne$normalizedNPXData, osloTwo$NPXData),
  metaData = rbind(osloOne$metaData, osloTwo$metaData),
  proteinData = osloOne$proteinData,
  sampleData = rbind(osloOne$sampleData, osloTwo$sampleData),
  excludedSampleData = rbind(osloOne$excludedSampleData, osloTwo$excludedSampleData)
) %>% addSurvival()


# Fetch WU
WU <- fetch(subset(metaData, origin == "WU")$studyID, selectedProteins, ignoreBelowLOD = FALSE) %>% addSurvival()



# Create scores
proteinsForScore <- list("healthy" = HD$NPXData, "patients" = osloMerged$NPXData) %>%
  compareTwoGroups(selectedProteins, pAdjustMethod = "BH") %>%
  subset(wilcoxAdj <= 0.05 & medianDiff > 0) %>%
  .$protein

osloMerged <- osloMerged %>% addScore(proteinsForScore, quantile = NA, name = "higherInPatientsQuantiles")
WU <- WU %>% addScore(proteinsForScore, quantile = NA, name = "higherInPatientsQuantiles")

# Optimized score
scoreOptimizationProteinOrder <- orderProteinsForOptimization(osloMerged, proteinsForScore)
scoreOptimization <- optimizeProteinCount(osloMerged, scoreOptimizationProteinOrder, "OS")

osloMerged <- osloMerged %>% addScore(scoreOptimization$optimalProteins, quantile = NA, name = "lowestPossibleFromOslo")
WU <- WU %>% addScore(scoreOptimization$optimalProteins, quantile = NA, name = "lowestPossibleFromOslo")


# Define colors
figureColors <- function(colors) {
  colorCodes <- c("blue"        = "#1d0e82",
                  "red"         = "#ed1c24",
                  "green"       = "#00a651",
                  "lightBlue"   = "#8071e3",
                  "lightRed"    = "#ff7378",
                  "lightGreen"  = "#66cc98",
                  "grey"        = "#c6c6c6",
                  "yellow"      = "#fff200",
                  "lightOrange" = "#ffc20e",
                  "greenYellow" = "#a6ce3a",
                  "orange"      = "#f7941d",
                  "cyan"        = "#00aeef",
                  "magenta"     = "#ec008c")
  as.vector(sapply(colors, function(color) {
    colorCodes[color]
  }))
}

colorsForKaplanMeierPlots <- c("score=low/mid"  = figureColors("blue"),
                               "score=low"      = figureColors("green"),
                               "score=mid"      = figureColors("blue"),
                               "score=high"     = figureColors("red"),
                               "score=veryhigh" = figureColors("black"))


