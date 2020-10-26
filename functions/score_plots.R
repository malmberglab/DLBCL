
histogramForScore <- function(dataset, scoreName, group, colors = NULL) {
  score <- dataset$score[[scoreName]]
  density <- density(score$numeric)
  groupLevels <- levels(score[[group]])
  cutoffs <- sapply(1:(length(groupLevels)-1), function(i) {
    mean(
      max(score$numeric[score[[group]] == groupLevels[i]]),
      min(score$numeric[score[[group]] == groupLevels[i+1]])
    )
  })

  df <- data.frame(x = density$x,
                   y = density$y,
                   part = factor(groupLevels[findInterval(density$x, cutoffs) + 1], levels = groupLevels))

  ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    geom_ribbon(aes(ymin = 0, ymax = y, fill = part)) +
    theme_classic() +
    scale_x_continuous(breaks = (0:10)*10, name = "Score", limits = c(0, length(dataset$score[[scoreName]]$proteins))) +
    scale_y_continuous(name = "Patients") +
    scale_fill_manual(values = colors, guide = FALSE) +
    theme(axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 16, color = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

scoreGroupsVsControlBoxPlots <- function(dataset, scoreName, scoreGroup, proteins = NULL, controlDataset, colors) {
  if (is.null(proteins)) proteins <- dataset$score[[scoreName]]$proteins
  
  plotData <- do.call(rbind, lapply(proteins, function(protein) {
    df <- rbind(
      data.frame(group = "healthy", value = controlDataset$NPXData[, protein]),
      data.frame(group = dataset$score[[scoreName]][[scoreGroup]], value = dataset$NPXData[, protein])
    )
    df$value <- df$value - median(controlDataset$NPXData[, protein])
    df$protein <- protein
    df$proteinName <- formatProteinName(protein)
    df
  }))

  colors <- c("healthy" = "gray80", colors)
  
  ggplot(plotData, aes(x = factor(proteinName, levels = formatProteinName(proteins)),
                       fill = factor(group, levels = names(colors)),
                       y = value)) +
    geom_hline(yintercept = 0, color = "gray80") +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual(values = colors, name = "Group") +
    scale_y_continuous(breaks = (-20:20)*2) +
    theme_classic() +
    labs(y = "NPX (normalized to healthy)", x = NULL) +
    theme(axis.text.y = element_text(size = 18, color = "black"),
          axis.text.x = element_text(size = 16, color = "black", angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top")
}
