calculateTSNEDifference <- function(tSNESample1, tSNESample2, n = 200) {
  range <- lapply(1:2, function(d) {
    c(
      min(min(tSNESample1[, d]), min(tSNESample2[, d])),
      max(max(tSNESample1[, d]), max(tSNESample2[, d]))
    )*1.15
  })
  group1 <- kde2d(tSNESample1[, 1], tSNESample1[, 2], lims = c(range[[1]], range[[2]]), n = n)
  group2 <- kde2d(tSNESample2[, 1], tSNESample2[, 2], lims = c(range[[1]], range[[2]]), n = n)
  
  if (!identical(group1$x, group2$x) || !identical(group1$y, group2$y)) stop("Coordinates not identical")
  
  diff12 <- group1
  diff12$z <- group1$z - group2$z
  
  rownames(diff12$z) <- diff12$x
  colnames(diff12$z) <- diff12$y
  
  diff12_melt <- melt(diff12$z, id.var = row.names(diff12))
  names(diff12_melt) <- c("x", "y", "diff")
  
  return(diff12_melt)
}

createTSNEDifferencePlot <- function(tSNE, group, colors = c("red", "blue")) {
  if (length(group) != nrow(tSNE)) stop("tSNE and groups should have the same length.")
  groups <- group %>% unique() %>% setdiff(NA)
  if (length(groups) != 2) stop("Exactly 2 groups should be provided.")
  
  cat(colors[1], " = more in group ", groups[2], "\n", sep="")
  cat(colors[2], " = more in group ", groups[1], "\n", sep="")
  
  diff <- calculateTSNEDifference(tSNESample1 = tSNE[!is.na(group) & group==groups[1], ],
                                  tSNESample2 = tSNE[!is.na(group) & group==groups[2], ])
  range <- max(abs(range(diff$diff))) * c(-1, 1)
  cat("Range: ", range[1], " to ", range[2], "\n", sep="")
  
  ggplot() +
    geom_tile(data = diff, aes(x = x, y = y, fill = diff)) +
    geom_density2d(data = NULL, mapping = aes(x = tSNE[,1], y = tSNE[,2]), color = "grey20", size = 0.2) +
    scale_fill_gradient2(low = colors[1], mid = "white", high = colors[2], midpoint = 0, limits = range) +
    xlim(range(tSNE[,1]) * 1.2) + 
    ylim(range(tSNE[,2]) * 1.2) +
    guides(fill = FALSE) +
    theme_void() +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
}


correlateReadoutsToScore <- function(data, readouts, method = "spearman") {
  do.call(rbind, lapply(readouts, function(readout) {
    corTest <- cor.test(data[!is.na(data$score), "score"],
                        data[!is.na(data$score), readout],
                        method = method)
    data.frame(
      readout = readout,
      estimate = corTest$est,
      method = method,
      estimateType = names(corTest$est),
      p = corTest$p.value
    )
  }))
}


plotReadoutScoreCorrelations <- function(correlations, orderByEstimate = TRUE, padTo = 16, padLength = 35,
                                         limits = c(NA, NA)) {
  if (orderByEstimate) {
    readoutFactors <- correlations$readout[order(-correlations$estimate)]
  } else {
    readoutFactors <- correlations$readout
  }
  
  if (padTo > nrow(correlations)) {
    readoutFactors <- c(readoutFactors, paste0("pad-", 1:(padTo - nrow(correlations))))
    readoutFactors[length(readoutFactors)] <- rep("x", padLength) %>% paste(collapse = "")
  }
  
  correlations$readout <- factor(correlations$readout, readoutFactors)
  
  
  asteriskColors <- c(
    "****" = figureColors("blue"),
    "***" = figureColors("magenta"),
    "**" = figureColors("red"),
    "*" = figureColors("orange"),
    "ns." = figureColors("grey")
  )
  
  ggplot(correlations, aes(x = readout, y = estimate, fill = factor(pToAsterisks(p), rev(names(asteriskColors))))) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_hline(yintercept = 0) +
    
    scale_fill_manual(values = asteriskColors, drop = FALSE, name = "p value") +
    scale_y_continuous(limits = limits, breaks = ((-5):5)/5) +
    scale_x_discrete(drop = FALSE) +
    theme_pubr() +
    labs(x = NULL, y = paste(correlations$method[1], correlations$estimateType[1])) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}


createProteinCytofCorrelationPlot <- function(dataset, readoutData, readouts, proteins, method = "spearman") {
  readoutData <- readoutData[readoutData$studyID %in% dataset$metaData$studyID, ]
  
  correlations <- do.call(rbind, lapply(readouts, function(readout) {
    do.call(rbind, lapply(proteins, function(protein) {
      readoutValues <- readoutData[, readout]
      proteinValues <- sapply(readoutData$studyID, function(studyID) {
        dataset$NPXData[dataset$metaData$studyID == studyID, protein]
      })
      corTest <- cor.test(readoutValues, proteinValues, method = method)
      data.frame(
        readout = readout,
        protein = protein,
        proteinName = formatProteinName(protein),
        estimate = corTest$estimate,
        p = corTest$p.value
      )
    }))
  }))

  proteinOrder <- correlations %>%
    pivot_wider(id_cols = "protein", names_from = "readout", values_from = "estimate") %>%
    as.data.frame() %>%
    { row.names(.) <- .$protein; .[, names(.) != "protein"] } %>%
    as.matrix() %>%
    { rownames(.)[hclust(dist(.))$order] }
  
  plot <- ggplot(correlations, aes(x = factor(proteinName, formatProteinName(proteinOrder)),
                                   y = factor(readout, rev(readouts)),
                                   fill = estimate)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-0.8, 0.8)) +
    coord_equal() +
    theme_classic() +
    labs(x = NULL, y = NULL) +
    theme(
      axis.ticks = element_line(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5)
    )
  
  return(list(
    statistics = correlations,
    plot = plot
  ))
}


