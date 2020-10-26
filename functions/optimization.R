orderProteinsForOptimization <- function(dataset, proteins) {
  list("low"=dataset$NPXData[dataset$score$higherInPatientsQuantiles$groupedBySecondTertile=="low/mid",],
    "high"=dataset$NPXData[dataset$score$higherInPatientsQuantiles$groupedBySecondTertile=="high",]) %>% 
    compareTwoGroups(proteins, pAdjustMethod="BH") %>%
    subset(wilcoxAdj <= 0.05 & medianDiff > 0) %>% { .[order(-.$medianDiff),]} %>%
    .$protein
}

optimizeProteinCount <- function(dataset, orderedProteins, survivalType, type = "groupScoresBySecondTertile") {
  optimization <- do.call(rbind, lapply(1:length(orderedProteins), function(p) {
    proteins <- orderedProteins[1:p]
    scores <- createScores(dataset, proteins, quantile=NA)
    if (type == "groupScoresBySecondTertile")
      scores <- scores %>% groupScoresBySecondTertile()
      
    coxreg <- summary(coxph(Surv(days, status) ~ score,
                            cbind(dataset$survival[[survivalType]], data.frame(score = scores))))
    data.frame(
      proteins = paste(proteins, collapse=", "),
      proteinCount = length(proteins),
      concordance = coxreg$concordance["C"],
      logrank = coxreg$sctest["pvalue"]
    )
  }))

  optimalCount <- optimization %>% subset(concordance == max(concordance)) %>% .$proteinCount %>% min()
  
  optimizationPlot <- ggplot(optimization, aes(x=proteinCount, y=concordance)) +
    geom_line() +
    ylim(c(0.5, 1)) +
    geom_vline(xintercept=optimalCount, color="red") +
    scale_x_continuous(breaks=(0:100)*5, minor_breaks=0:100) +
    theme_linedraw()
  
  list(
    optimization = optimization,
    optimalCount = optimalCount,
    optimalProteins = orderedProteins[1:optimalCount],
    plot = optimizationPlot
  )
}
