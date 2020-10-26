
#######################################################################################################################
# Figure 2A:
#######################################################################################################################

pdf(file.path(outputDir, "figures", "Figure 2A.pdf"), height = 14, width = 14)
data <- do.call(data.frame, sapply(selectedProteins, function(protein) {
  values <- c(HD$NPXData[, protein], osloMerged$NPXData[, protein])
  return(values - median(values))
}, simplify = FALSE, USE.NAMES = TRUE))
largeHeatmap <- gplots::heatmap.2(t(data),
                                  trace         = "none",
                                  scale         = "none",
                                  col           = colorRampPalette(c("blue", "white", "red"))(100),
                                  labRow        = formatProteinName(selectedProteins),
                                  ColSideColors = c(rep("#CCCCCC", nrow(HD$NPXData)),
                                                    rep("black", nrow(osloMerged$NPXData))),
                                  distfun       = function(x) dist(x, method = "euclidean"),
                                  hclustfun     = function(x) hclust(x, method = "complete"),
                                  margins       = c(10, 10))
dev.off()


#######################################################################################################################
# Figure 2B:
#######################################################################################################################

clusters <- largeHeatmap$colDendrogram
clusterStudyIDs <- lapply(seq_along(clusters), function(cluster) {
  sapply(unlist(clusters[[cluster]]), function(i) {
    c(HD$metaData$studyID, osloMerged$metaData$studyID)[i]
  })
})

osloMerged$score$fromHeatmap <- list(cluster = sapply(osloMerged$metaData$studyID, function(studyID) {
  for (cluster in seq_along(clusterStudyIDs)) {
    if (studyID %in% clusterStudyIDs[[cluster]])
      return(paste("Cluster", cluster))
  }
}))

plotScoreSurvival(osloMerged, "OS", "fromHeatmap", "cluster",
                  colors=c("score=Cluster 1" = figureColors("blue"), "score=Cluster 2" = figureColors("red"))) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure 2B.pdf"), width = 5.5, height = 8)


rm(data, largeHeatmap, clusters, clusterStudyIDs)