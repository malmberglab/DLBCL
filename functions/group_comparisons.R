
compareGroups <- function(groups, proteins, pAdjustMethod = NULL, normalizeHeatmapToMedian = FALSE,
                          normalizeHeatmapToMedianOfDataset = NULL, boxplotColors = NULL) {

  if (normalizeHeatmapToMedian && !is.null(normalizeHeatmapToMedianOfDataset))
    stop("normalizeHeatmapToMedian cannot be TRUE if normalizeHeatmapToMedianOfDataset is not NULL")

  return <- list()

  result <- do.call(rbind, lapply(proteins, function(protein) {
    values <- sapply(names(groups), function(groupName) {
      groups[[groupName]][, protein]
    }, simplify = FALSE, USE.NAMES = TRUE)
    
    df <- data.frame(protein = protein, proteinName = formatProteinName(protein), groupCount = length(groups))
    for (group in names(groups)) {
      df[, paste0(group, "Count")] <- length(values[[group]])
      df[, paste0(group, "Mean")] <- mean(values[[group]])
      df[, paste0(group, "Median")] <- median(values[[group]])
      df[, paste0(group, "SD")] <- sd(values[[group]])
    }

    if (length(groups) == 2) {

      df$meanDiff <- df[, paste0(names(groups)[2], "Mean")] - df[, paste0(names(groups)[1], "Mean")]
      df$medianDiff <- df[, paste0(names(groups)[2], "Median")] - df[, paste0(names(groups)[1], "Median")]

      wilcox <- wilcox.test(values[[1]], values[[2]])
      tTest <- t.test(values[[1]], values[[2]])
  
      df$wilcox <- wilcox.test(values[[1]], values[[2]])$p.value
      df$wilcoxSign <- pToAsterisks(df$wilcox)
      df$tTest <- t.test(values[[1]], values[[2]])$p.value
      df$tTestSign <- pToAsterisks(df$tTest)

    } else if (length(groups) > 2) {

      data <- do.call(rbind, lapply(names(values), function(group) {
        if (sum(values[[group]]) < 1) return(data.frame())
        data.frame(
          x = values[[group]],
          g = group
        )
      }))
  
      kruskal <- kruskal.test(data$x, data$g)
      
      df$kruskal <- kruskal$p.value
      df$kruskalSign <- pToAsterisks(kruskal$p.value)
    }

    df
  }))

  if (!is.null(pAdjustMethod)) {
    if ("wilcox" %in% names(result)) {
      result$wilcoxAdj <- p.adjust(result$wilcox, method = pAdjustMethod)
      result$wilcoxAdjSign <- pToAsterisks(result$wilcoxAdj)
      result$tTestAdj <- p.adjust(result$tTest, method = pAdjustMethod)
      result$tTestAdjSign <- pToAsterisks(result$tTestAdj)
    } else if ("kruskal" %in% names(result)) {
      result$kruskalAdj <- p.adjust(result$kruskal, method = pAdjustMethod)
      result$kruskalAdjSign <- pToAsterisks(result$kruskalAdj)
    }
  }

  return$statistics <- result


  
  longData <- do.call(rbind, lapply(proteins, function(protein) {
    df <- do.call(rbind, lapply(names(groups), function(groupName) {
      data.frame(group = groupName, value = groups[[groupName]][, protein])
    }))
    if (!is.null(normalizeHeatmapToMedianOfDataset)) {
      df$value <- df$value - median(normalizeHeatmapToMedianOfDataset[, protein])
    }
    else if (normalizeHeatmapToMedian) {
      df$value <- df$value - median(df$value)
    }
    df$protein <- protein
    df$proteinName <- formatProteinName(protein)
    return(df)
  }))


  # Heat map
  heatmap <- ggplot(longData, aes(x = factor(group, names(groups)),
                                  y = factor(proteinName, rev(formatProteinName(proteins))),
                                  fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.ticks = element_line(color = "black"),
      axis.text.x = element_text(color="black"),
      axis.text.y = element_text(color="black"))

  return$heatmap <- heatmap



  if (is.null(boxplotColors))
    boxplotColors <- figureColors(c("lightRed", "lightBlue", "lightGreen"))

  # Box plot
  boxplot <- ggplot(longData, aes(y = value,
                                  x = factor(proteinName, levels = rev(formatProteinName(proteins))),
                                  fill = factor(group, levels = rev(names(groups)))))
  
  if (!normalizeHeatmapToMedian && is.null(normalizeHeatmapToMedianOfDataset)) {
    boxplot <- boxplot + ylab("NPX")
  } else {
    boxplot <- boxplot + geom_hline(yintercept=0, color="gray80")
    if (normalizeHeatmapToMedian) {
      boxplot <- boxplot + ylab(paste0("NPX (normalized to median of all samples)"))
    } else {
      boxplot <- boxplot + ylab(paste0("NPX (normalized to median of separate dataset)"))
    }
  }

  boxplot <- boxplot +
    geom_boxplot(outlier.size=0.3) +
    scale_y_continuous(breaks=(-20:20)) +
    scale_fill_manual(values=boxplotColors[1:length(groups)] %>% setNames(names(groups)), name="Group") +
    theme_classic() +
    coord_flip() +
    xlab(NULL) +
    theme(axis.text = element_text(color="black"),
          legend.position = "top",
          axis.ticks = element_line(color = "black"))
  

  return$boxplot <- boxplot





  return(return)
}




compareTwoGroups <- function(groups, proteins, pAdjustMethod=NULL) {
  if (length(groups) != 2) stop("Exactly two groups must be provided")

  result <- do.call(rbind, lapply(proteins, function(protein) {
    values <- sapply(names(groups), function(groupName) {
      groups[[groupName]][, protein]
    }, simplify=FALSE, USE.NAMES=TRUE)
    
    wilcox <- wilcox.test(values[[1]], values[[2]])
    tTest <- t.test(values[[1]], values[[2]])
    
    df <- data.frame(protein = protein, proteinName = formatProteinName(protein))
    for (group in names(groups)) {
      df[, paste0(group, "Count")] <- length(values[[group]])
      df[, paste0(group, "Mean")] <- mean(values[[group]])
      df[, paste0(group, "Median")] <- median(values[[group]])
      df[, paste0(group, "SD")] <- sd(values[[group]])
    }

    df$meanDiff <- df[, paste0(names(groups)[2], "Mean")] - df[, paste0(names(groups)[1], "Mean")]
    df$medianDiff <- df[, paste0(names(groups)[2], "Median")] - df[, paste0(names(groups)[1], "Median")]

    df$wilcox <- wilcox.test(values[[1]], values[[2]])$p.value
    df$wilcoxSign <- pToAsterisks(df$wilcox)
    df$tTest <- t.test(values[[1]], values[[2]])$p.value
    df$tTestSign <- pToAsterisks(df$tTest)
    df
  }))

  if (!is.null(pAdjustMethod)) {
    result$wilcoxAdj <- p.adjust(result$wilcox, method=pAdjustMethod)
    result$wilcoxAdjSign <- pToAsterisks(result$wilcoxAdj)
    result$tTestAdj <- p.adjust(result$tTest, method=pAdjustMethod)
    result$tTestAdjSign <- pToAsterisks(result$tTestAdj)
  }

  result
}


twoGroupsBoxPlotComparison <- function(groups, proteins, normalizeTo = NULL, order = NULL,
                                       colors = c("gray80", "red")) {
  if (length(groups) != 2) stop("Exactly two groups must be provided")

  df <- do.call(rbind, lapply(proteins, function(protein) {
    values <- sapply(names(groups), function(groupName) {
      groups[[groupName]][, protein] - if (!is.null(normalizeTo)) median(groups[[normalizeTo]][, protein]) else 0
    }, simplify=FALSE, USE.NAMES=TRUE)
    
    data.frame(
      protein = protein,
      proteinName = formatProteinName(protein),
      group = c(rep(names(groups)[1], length(values[[1]])), rep(names(groups)[2], length(values[[2]]))),
      value = c(values[[1]], values[[2]])
    )
  }))
  
  proteinOrder <- groups %>%
    compareTwoGroups(selectedProteins, pAdjustMethod = "BH") %>%
    { .[order(.$medianDiff), "proteinName"] }

  plot <- ggplot(df, aes(x = factor(proteinName, levels = proteinOrder),
                         fill = factor(group, levels = rev(names(groups))),
                         y = value)) +
    geom_boxplot(outlier.size=0.3) +
    scale_y_continuous(breaks=(-20:20)) +
    scale_fill_manual(values=colors %>% setNames(names(groups)), name="Group") +
    theme_classic() +
    coord_flip() +
    xlab(NULL) +
    theme(axis.text=element_text(color="black"), legend.position="top")
  if (is.null(normalizeTo)) {
    plot <- plot + ylab("NPX")
  } else {
    plot <- plot + geom_hline(yintercept=0, color="gray80") +
      ylab(paste0("NPX (normalized to ", normalizeTo, " median)"))
  }

  plot
}


metaDataMultipleGroupComparison <- function(dataset, groups, proteins=selectedProteins, customPlotLimits=c(NA, NA),
  groupLevelCount = 8, statisticsDataset = NULL, statisticsGroups = NULL) {

  if ("normalizedNPXData" %in% names(dataset)) {
    message("Using normalized data")
    dataset$NPXData <- dataset$normalizedNPXData
  }
  
  if (is.null(statisticsDataset)) {
    statisticsDataset <- dataset
    statisticsGroups <- groups
  }
  if ("normalizedNPXData" %in% names(statisticsDataset)) {
    message("Using normalized data for statistics")
    statisticsDataset$NPXData <- statisticsDataset$normalizedNPXData
  }

  df <- do.call(rbind, lapply(proteins, function(protein) {
    medians <- sapply(names(groups), function(group) {
      median(dataset$NPXData[groups[[group]], protein], na.rm=TRUE)
    })
    
    data.frame(
      protein = protein,
      proteinName = formatProteinName(protein),
      group = names(groups),
      median = medians,
      #medianNorm = medians - medians[1]
      medianNorm = medians - median(dataset$NPXData[, protein])
      #medianNorm = medians - median(HD$normalizedNPXData[, protein])
    )
  }))
  
  statistics <- do.call(rbind, lapply(proteins, function(protein) {
    data <- do.call(rbind, lapply(names(statisticsGroups), function(group) {
      if (sum(statisticsGroups[[group]]) < 1) return(data.frame())
      data.frame(
        x = statisticsDataset$NPXData[statisticsGroups[[group]], protein],
        g = group
      )
    }))

    kruskal <- kruskal.test(data$x, data$g)
    
    data.frame(
      protein = protein,
      proteinName = formatProteinName(protein),
      pKruskal = kruskal$p.value,
      pKruskalSign = pToAsterisks(kruskal$p.value)
    )
  }))
  
  includedProteins <- statistics[statistics$pKruskal <= 0.05, "protein"]
  
  dfMatrix <- as.matrix(do.call(rbind, lapply(includedProteins, function(protein) {
    sapply(names(groups), function(group) {
      median(statisticsDataset$NPXData[statisticsGroups[[group]], protein], na.rm=TRUE)
    })
  })))
  row.names(dfMatrix) <- includedProteins
  dfMatrix <- dfMatrix - apply(statisticsDataset$NPXData[,includedProteins], 2, median)
  
  proteinOrder <- includedProteins[hclust(dist(dfMatrix))$order]
  proteinOrderWithPadding <- c(formatProteinName(proteinOrder), "padpadpadpad")
  
  groupLevels <- names(groups)
  i <- 1
  while (length(groupLevels) < groupLevelCount) {
    groupLevels <- c(groupLevels, paste(rep(" ", i), collapse=""))
    i <- i+1
  }


  
  plot <- ggplot(df[df$protein %in% includedProteins,], aes(x = factor(group, levels = groupLevels),
                                                            y = factor(proteinName, levels = proteinOrderWithPadding),
                                                            fill = medianNorm)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, name="diff",
                         limits=customPlotLimits, na.value="black") +
    theme_classic() +
    ylab(NULL) + xlab(NULL) +
    coord_equal() +
    theme(axis.text=element_text(size=16, color="black")) +
    scale_x_discrete(drop=FALSE) +
    scale_y_discrete(drop=FALSE)
  
  list(
    data = df,
    statistics = statistics[sapply(rev(proteinOrder), function(protein) { which(statistics$protein==protein) }),],
    plot = plot,
    matrix = dfMatrix
  )
}


