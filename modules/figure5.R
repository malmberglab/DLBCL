
cytofSamples <- read.csv(file.path(dataDir, "fcs.csv"), sep = ";", dec = ",")
cytofSamples$sampleID <- 1:nrow(cytofSamples)
cytofSamples$sampleGroup <- ifelse(cytofSamples$type == "healthy", 1, 2)



#######################################################################################################################
# Load data for figure 5A-B
#######################################################################################################################

seed <- 2020
eventsPerSample <- 1000
params <- c("CD3-170", "CD4-145", "CD8-141", "CD14-EQ-175", "CD16-148", "CD19-171", "CD28-160", "CD38-144",
            "CD45RA-155", "CD56-EQ-176", "CD57-EQ-142", "CCR7-167", "HLA-DR-173", "NKG2A-169", "IgD-146", "TCRgd-152")


# Read FCS files
fcsCacheFilePath <- file.path(dir, "cache", "Figure 5 t-SNE.dat")
if (file.exists(fcsCacheFilePath)) {
  load(fcsCacheFilePath)
} else {
  cytofData <- do.call(rbind, lapply(1:nrow(cytofSamples), function(i) {
    cat(cytofSamples[i, "fcsFile"], "\n")
    fcs <- read.FCS(file.path(dataDir, "fcs", cytofSamples[i, "fcsFile"]))
    fcsData <- as.data.frame(exprs(fcs))
    names(fcsData) <- as.character(parameters(fcs)$desc)
    fcsData <- fcsData[, params]
    fcsData <- asinh(fcsData / 5)
    
    set.seed(seed * cytofSamples[i, "fcsID"])
    fcsData <- fcsData[sample(1:nrow(fcsData), eventsPerSample), ]
    fcsData$sampleID <- cytofSamples[i, "fcsID"]
    fcsData$sampleGroup <- cytofSamples[i, "sampleGroup"]
    
    return(fcsData)
  }))

  set.seed(seed)
  cytofTsne <- Rtsne(cytofData[, params], perplexity = 30, theta = 0.5, verbose = TRUE)$Y

  cytofData$tsne1 <- cytofTsne[, 1]
  cytofData$tsne2 <- cytofTsne[, 2]
  
  rm(cytofTsne)

  save(cytofData, file = fcsCacheFilePath)
}

rm(fcsCacheFilePath)



#######################################################################################################################
# Figure 5A
#######################################################################################################################

set.seed(seed)
eventOrder <- sample(1:nrow(cytofData), nrow(cytofData))

intensityPlots <- lapply(params, function(param) {
  ggplot() +
    geom_point(data = cytofData[eventOrder, ],
               mapping = aes(x = tsne1, y = tsne2, color = cytofData[eventOrder, param]),
               alpha = I(0.2), size = 0.3, shape = 16) +
    geom_density2d(data = cytofData[eventOrder,],
                   mapping = aes(x = tsne1, y = tsne2), color = "black", size = 0.1) +
    scale_color_gradientn(colors = c("blue", "turquoise", "yellow", "orange", "red"),
                          guide = FALSE, limits = c(0, NA)) +
    xlim(min(cytofData$tsne1) * 1.15, max(cytofData$tsne1) * 1.15) +
    ylim(min(cytofData$tsne2) * 1.15, max(cytofData$tsne2) * 1.15) +
    ggtitle(param) +
    theme_void() +
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          aspect.ratio = 1,
          plot.margin = unit(rep(-0.15, 4), "cm"))
})

arrangeGrob(grobs = intensityPlots, ncol = 6) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 5A.png"), width = 8, height = 5.2)



#######################################################################################################################
# Figure 5B
#######################################################################################################################

createTSNEDifferencePlot(cytofData[, c("tsne1", "tsne2")], cytofData$sampleGroup) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 5B.png"), width = 3, height = 3)



rm(seed, eventsPerSample, params, cytofData, eventOrder, intensityPlots, calculateTSNEDifference,
   createTSNEDifferencePlot)



#######################################################################################################################
# Load data for figure 5C-E
#######################################################################################################################

cytofReadouts <- read.csv(file.path(dataDir, "fcsData.csv"), sep = ";", dec = ",", row.names = 1,
                          stringsAsFactors = FALSE)

cytofReadouts$type <- sapply(row.names(cytofReadouts), function(fcsFile) {
  cytofSamples[cytofSamples$fcsFile == fcsFile, "type"]
}) %>% factor(c("healthy", "patient"))

cytofReadouts$studyID <- sapply(row.names(cytofReadouts), function(fcsFile) {
  cytofSamples[cytofSamples$fcsFile == fcsFile, "studyID"]
})

cytofReadouts$score <- sapply(cytofReadouts$studyID, function(studyID) {
  if (studyID %in% osloMerged$metaData$studyID) {
    osloMerged$score$lowestPossibleFromOslo$numeric[osloMerged$metaData$studyID == studyID]
  } else {
    NA
  }
})

cytofReadouts$scoreGroup <- sapply(cytofReadouts$studyID, function(studyID) {
  if (studyID %in% osloMerged$metaData$studyID) {
    as.character(osloMerged$score$lowestPossibleFromOslo$groupedBySecondTertile[osloMerged$metaData$studyID == studyID])
  } else {
    NA
  }
}) %>% factor(c("low/mid", "high"))



#######################################################################################################################
# Figure 5C
#######################################################################################################################


plot <- ggplot(cytofReadouts, aes(x = type, y = MDSCs, fill = type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("healthy" = figureColors("grey"), "patient" = figureColors("lightBlue")),
                    name = NULL, guide = guide_legend(nrow = 2)) +
  scale_y_continuous(breaks = (0:10)/10, labels = paste0((0:10)*10, "%"),
                     name = "% CD14+ HLA-DRlow/-\nCD3- CD19- in CD45+") +
  ggtitle(paste0("p = ", round(wilcox.test(MDSCs ~ type, cytofReadouts)$p.value, 8))) +
  theme_pubr() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"))

plot %>% ggsave(file = file.path(outputDir, "figures", "Figure 5C.pdf"), width = 3, height = 4, useDingbats = FALSE)

rm(plot)



#######################################################################################################################
# Figure 5D
#######################################################################################################################

TCellMarkers <- c("Percent_CD57",
                  "Mean_CD38",
                  "Percent_HLADR",
                  "Mean_CD2",
                  "Mean_CD27",
                  "Percent_PD1",
                  "Percent_Ki67",
                  "Percent_CD71",
                  "Percent_CD25",
                  "Mean_CD98",
                  "Mean_CD127")

plotReadouts <- list(
  c("CD56bright_NK", "CD56dim_NK", "CD4_T", "CD8_T", "gd_T", "B", "Monocytes", "MDSCs"),
  paste0("Percent_", c("Tn", "Tcm", "Ttm", "Tem", "Temra"), "_in_CD4_T_cells"),
  paste0("Percent_", c("Tn", "Tcm", "Ttm", "Tem", "Temra"), "_in_CD8_T_cells"),
  paste0(TCellMarkers, "_in_CD4_T_cells"),
  c(paste0(TCellMarkers, "_in_CD8_T_cells"), "Percent_PD1_in_Ki67_CD8_T_cells")
)

correlationPlots <- lapply(plotReadouts, function(readouts) {
  correlateReadoutsToScore(cytofReadouts, readouts, method = "spearman") %>%
    plotReadoutScoreCorrelations(orderByEstimate = TRUE, limits = 0.62 * c(-1,1))
})

arrangeGrob(grobs = correlationPlots, nrow = 1) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 5D.pdf"), width = 18, height = 5.2, useDingbats = FALSE)




rm(cytofReadouts, correlations)
