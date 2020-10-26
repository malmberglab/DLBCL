
healthyPatientComparison <- list("healthy" = HD$NPXData, "patients" = osloMerged$NPXData) %>%
  compareTwoGroups(selectedProteins, pAdjustMethod = "BH")
  

#######################################################################################################################
# Figure 1A:
#######################################################################################################################

GOCacheFilePath <- file.path(dir, "cache", "Figure 1A GO keywords.dat")
if (file.exists(GOCacheFilePath)) {
  load(GOCacheFilePath)
} else {
  GOKeywords <- fetchGOKeywordsForProteins(selectedProteins)
  save(GOKeywords, file = GOCacheFilePath)
}

# GO Keywords - biological process (KW-9999)
filteredGOKeywords <- filterGOKeywords(GOKeywords, "9999") %>% formatGOKeywordNames()
selectedKeywords <- as.data.frame(table(unlist(filteredGOKeywords)), stringsAsFactors = FALSE) %>%
  { .[.$Freq > 1,] } %>% { .[order(-.$Freq), "Var1"] }
keywordPlotData <- do.call(rbind, lapply(selectedProteins, function(protein) {
  data.frame(
    protein = protein,
    proteinName = formatProteinName(protein),
    keyword = selectedKeywords,
    included = selectedKeywords %in% filteredGOKeywords[[protein]]
  )
}))

orderedProteins <- healthyPatientComparison[order(healthyPatientComparison$medianDiff), "protein"]

keywordPlotData$protein     <- factor(keywordPlotData$protein,     levels = orderedProteins)
keywordPlotData$proteinName <- factor(keywordPlotData$proteinName, levels = formatProteinName(orderedProteins))
keywordPlotData$keyword     <- factor(keywordPlotData$keyword,     levels = selectedKeywords)

keywordPlot <- ggplot(keywordPlotData, aes(x = keyword, y = proteinName, fill = included)) +
  geom_tile(color = "#ffffff") +
  scale_fill_manual(values = c("TRUE" = figureColors("green"), "FALSE" = "#e6e6e6")) +
  theme_classic() +
  coord_equal() +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.title = element_blank())

keywordPlot %>% ggsave(file = file.path(outputDir, "figures", "Figure 1A.pdf"), width = 6, height = 16, useDingbats = FALSE)


rm(GOCacheFilePath, GOKeywords, filteredGOKeywords, selectedKeywords, keywordPlotData, orderedProteins, keywordPlot)


#######################################################################################################################
# Figure 1B:
#######################################################################################################################

write.table(healthyPatientComparison[order(-healthyPatientComparison$medianDiff),],
            file = file.path(outputDir, "statistics", "Figure 1B.csv"), sep=";", dec=",", row.names = FALSE)

list("healthy" = HD$NPXData, "patients" = osloMerged$NPXData) %>%
  twoGroupsBoxPlotComparison(selectedProteins, normalizeTo = "healthy", colors = figureColors(c("grey", "red"))) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 1B.pdf"), width = 8, height = 16, useDingbats = FALSE)



#######################################################################################################################
# Figure 1C:
#######################################################################################################################

list(
  "30+" = osloMerged$metaData$age_group == "30-39",
  "40+" = osloMerged$metaData$age_group == "40-49",
  "50+" = osloMerged$metaData$age_group == "50-59",
  "60+" = osloMerged$metaData$age_group == "60-69",
  "70+" = osloMerged$metaData$age_group == "70-79",
  "80+" = osloMerged$metaData$age_group == "80-89"
) %>% { metaDataMultipleGroupComparison(osloMerged, ., selectedProteins, customPlotLimits = c(-4, 4)) } %>% {
  ggsave(.$plot, file = file.path(outputDir, "figures", "Figure 1C.pdf"), width = 4.5, height = 20, useDingbats = FALSE)
  write.table(.$statistics, file = file.path(outputDir, "statistics", "Figure 1C.csv"), sep=";", dec=",", row.names = FALSE)
}



#######################################################################################################################
# Figure 1D:
#######################################################################################################################

list(
  "I"   = osloMerged$metaData$stage == "I",
  "II"  = osloMerged$metaData$stage == "II",
  "III" = osloMerged$metaData$stage == "III",
  "IV"  = osloMerged$metaData$stage == "IV"
) %>% { metaDataMultipleGroupComparison(osloMerged, ., selectedProteins, customPlotLimits = c(-4, 4)) } %>% {
  ggsave(.$plot, file = file.path(outputDir, "figures", "Figure 1D.pdf"), width = 4.5, height = 20, useDingbats = FALSE)
  write.table(.$statistics, file = file.path(outputDir, "statistics", "Figure 1D.csv"), sep=";", dec=",", row.names = FALSE)
}



#######################################################################################################################
# Figure 1E:
#######################################################################################################################

list(
  "0" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 0,
  "1" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 1,
  "2" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 2,
  "3" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 3,
  "4" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 4,
  "5" = !is.na(osloMerged$metaData$IPI) & osloMerged$metaData$IPI == 5
) %>% { metaDataMultipleGroupComparison(osloMerged, ., selectedProteins, customPlotLimits = c(-4, 4)) } %>% {
  ggsave(.$plot, file = file.path(outputDir, "figures", "Figure 1E.pdf"), width = 4.5, height = 20, useDingbats = FALSE)
  write.table(.$statistics, file = file.path(outputDir, "statistics", "Figure 1E.csv"), sep=";", dec=",", row.names = FALSE)
}



#######################################################################################################################
# Figure 1F
#######################################################################################################################

BSymptomsComparison <- do.call(rbind, lapply(selectedProteins, function(protein) {
  osloA <- osloMerged$NPXData[osloMerged$metaData$stage_AB == "A", protein]
  osloB <- osloMerged$NPXData[osloMerged$metaData$stage_AB == "B", protein]
  WUA <- WU$NPXData[WU$metaData$stage_AB == "A", protein]
  WUB <- WU$NPXData[WU$metaData$stage_AB == "B", protein]
  
  data.frame(
    protein = protein,
    proteinName = formatProteinName(protein),
    osloDiff = mean(osloB) - mean(osloA),
    WUDiff = mean(WUB) - mean(WUA),
    osloP = wilcox.test(osloA, osloB)$p.value,
    WUP = wilcox.test(WUA, WUB)$p.value,
    osloCount = length(osloA) + length(osloB),
    osloCountA = length(osloA),
    osloCountB = length(osloB),
    WUCount = length(WUA) + length(WUB),
    WUCountA = length(WUA),
    WUCountB = length(WUB)
  )
}))

BSymptomsComparison$osloPAdj <- p.adjust(BSymptomsComparison$osloP, method = "BH")
BSymptomsComparison$WUPAdj   <- p.adjust(BSymptomsComparison$WUP, method = "BH")

BSymptomsComparison$significant <- sapply(BSymptomsComparison$protein, function(protein) {
  osloPAdj <- BSymptomsComparison[BSymptomsComparison$protein == protein, "osloPAdj"]
  WUPAdj <- BSymptomsComparison[BSymptomsComparison$protein == protein, "WUPAdj"]
  
  if (osloPAdj <= 0.05 && WUPAdj <= 0.05) {
    return("both")
  } else if (osloPAdj <= 0.05) {
    return("only Oslo")
  } else if (WUPAdj <= 0.05) {
    return("only St. Louis")
  } else {
    return("none")
  }
})

plot <- ggplot(BSymptomsComparison, aes(x = osloDiff, y = WUDiff, color = osloDiff >= 1 & WUDiff >= 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  scale_color_manual(values = c("TRUE" = figureColors("red"), "FALSE" = figureColors("grey")), guide = FALSE) +
  geom_text_repel(data = subset(BSymptomsComparison, osloDiff >= 1 & WUDiff >= 1),
                  aes(label = proteinName), size = 5, color = "black", box.padding = unit(0.8, "lines"),
                  point.padding = unit(0.5, "lines")) +
  scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5), limits = c(-0.5, 2.5),
                     name = "B symptoms - no B symptoms (NPX)\nOslo") +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5), limits = c(-0.5, 2.5),
                     name = "St. Louis\nB-symptoms - no B symptoms (NPX)") +
  theme_pubr() +
  theme(aspect.ratio = 1, axis.ticks = element_line(color = "black"))

plot %>% ggsave(file = file.path(outputDir, "figures", "Figure 1F.pdf"), width = 4, height = 5, useDingbats = FALSE)

rm(BSymptomsComparison, plot)


#######################################################################################################################
# Figure 1G
#######################################################################################################################

commonPrePostPatients <- intersect(osloOnePost$metaData$studyID, osloOne$metaData$studyID)

prePostComparisonPaired <- do.call(rbind, lapply(selectedProteins, function(protein) {
  preValues <- sapply(commonPrePostPatients, function(studyID) {
    osloOne$normalizedNPXData[osloOne$metaData$studyID == studyID, protein]
  })
  postValues <- sapply(commonPrePostPatients, function(studyID) {
    osloOnePost$normalizedNPXData[osloOnePost$metaData$studyID == studyID, protein]
  })
  
  wilcox <- wilcox.test(preValues, postValues, paired = TRUE)
  
  data.frame(
    protein = protein,
    proteinName = formatProteinName(protein),
    preCount = length(preValues),
    postCount = length(postValues),
    preMean = mean(preValues),
    postMean = mean(postValues),
    preMedian = median(preValues),
    postMedian = median(postValues),
    meanDiff = mean(postValues) - mean(preValues),
    medianDiff = median(postValues) - median(preValues),
    preSD = sd(preValues),
    postSD = sd(postValues),
    wilcoxPaired = wilcox$p.value,
    wilcoxPairedSign = pToAsterisks(wilcox$p.value)
  )
}))
prePostComparisonPaired$wilcoxPairedAdj <- p.adjust(prePostComparisonPaired$wilcoxPaired, method="BH")
prePostComparisonPaired$wilcoxPairedAdjSign <- pToAsterisks(prePostComparisonPaired$wilcoxPairedAdj)

# Make a data frame from all protein values from all common patients. Normalize to healthy donors
prePostDataNorm <- do.call(rbind, lapply(selectedProteins, function(protein) {
  healthyMedian <- median(HD$NPXData[, protein])
  preValues <- sapply(commonPrePostPatients, function(studyID) {
    osloOne$normalizedNPXData[osloOne$metaData$studyID == studyID, protein]
  }) - healthyMedian
  postValues <- sapply(commonPrePostPatients, function(studyID) {
    osloOnePost$normalizedNPXData[osloOnePost$metaData$studyID == studyID, protein]
  }) - healthyMedian
  
  data.frame(
    protein = protein,
    proteinName = formatProteinName(protein),
    group = c(rep("pre", length(preValues)), rep("post", length(postValues))),
    value = c(preValues, postValues)
  )
}))

healthyMedians <- apply(HD$NPXData[,selectedProteins], 2, median)

prePostDotData <- data.frame(
  protein = rep(prePostComparisonPaired$protein, 2),
  proteinName = formatProteinName(rep(prePostComparisonPaired$protein, 2)),
  timepoint = rep(c("pre", "post"), each=nrow(prePostComparisonPaired)),
  median = c(prePostComparisonPaired$preMedian - healthyMedians, prePostComparisonPaired$postMedian - healthyMedians), 
  sign = c(prePostComparisonPaired$wilcoxPairedAdjSign, rep("", nrow(prePostComparisonPaired)))
)

plot <- ggplot(prePostDotData, aes(
  x = factor(proteinName,
             levels = prePostComparisonPaired[order(-prePostComparisonPaired$preMedian+healthyMedians), "proteinName"]),
  color = factor(timepoint, levels=c("pre", "post")),
  y = median, 
  label = sign
)) +
  geom_hline(yintercept = 0, color = "gray80") +
  geom_point() +
  geom_text(y = 4.7, color = "black", angle = 90, size = 3) + 
  scale_color_manual(values = c("pre"  = figureColors("blue"), "post" = figureColors("red")), name = "Group") +
  scale_y_continuous(breaks = (-20:20), limits = c(NA, 5)) +
  theme_classic() +
  labs(y = "NPX (normalized to healthy median)", x = NULL) +
  theme(axis.text = element_text(color="black"),
        legend.position = "top",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

plot %>% ggsave(file = file.path(outputDir, "figures", "Figure 1G.pdf"), width = 12, height = 4, useDingbats = FALSE)

prePostComparisonPaired[order(-prePostComparisonPaired$preMedian+healthyMedians),] %>%
  write.table(file = file.path(outputDir, "statistics", "Figure 1G.csv"), sep = ";", dec = ",", row.names = FALSE)


rm(commonPrePostPatients, prePostComparisonPaired, prePostDataNorm, healthyMedians, prePostDotData, plot)


