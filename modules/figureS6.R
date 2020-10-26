
#######################################################################################################################
# Figure S6:
#######################################################################################################################

data <- do.call(data.frame, sapply(selectedProteins, function(protein) {
  values <- c(HD$NPXData[, protein], osloMerged$NPXData[, protein])
  return(values - median(values))
}, simplify = FALSE, USE.NAMES = TRUE))


hclustCut <- function(x, k) {
  list(cluster = cutree(hclust(dist(x, method = "euclidean"), method = "complete"), k = k))
}

gapPlot <- factoextra::fviz_gap_stat(cluster::clusGap(data, FUNcluster = hclustCut, K.max = 10),
                          maxSE = list(method = "firstmax", SE.factor = 1), linecolor = "black") +
  ggtitle(NULL) +
  theme(aspect.ratio = 1)

silhouettePlot <- factoextra::fviz_nbclust(data, FUN = hclustCut, method = "silhouette", k.max = 10,
	                                       linecolor = "black") +
  ggtitle(NULL) +
  theme(aspect.ratio = 1)

arrangeGrob(grobs = list(gapPlot, silhouettePlot), ncol = 2) %>%
  ggsave(file = file.path(outputDir, "figures", "Figure S6.pdf"), width = 8.5, height = 3.5)


rm(data, hclustCut, gapPlot, silhouettePlot)
