performClinicalDataComparison <- function(clinicalData, colors) {
  groups <- levels(clinicalData$scoreGroup)

  clinicalDataFractions <- do.call(rbind.data.frame, lapply(names(clinicalData)[names(clinicalData) != "scoreGroup"],
                                                            function(var) {
    df <- as.data.frame(table(clinicalData[, c(var, "scoreGroup")]), stringsAsFactors=FALSE)
    names(df) <- c("option", "group", "count")
    
    df$fraction <- sapply(1:nrow(df), function(row) {
      df[row, "count"] / sum(df[df$group==df[row, "group"], "count"])
    })
    
    df$variable <- var
    df
  }))
  
  plots <- lapply(unique(clinicalDataFractions$variable), function(variable) {
    ggplot(clinicalDataFractions[clinicalDataFractions$variable==variable,],
           aes(x=factor(group, levels=rev(groups)), y=fraction, fill=factor(option, levels=rev(names(colors))))) +
      geom_col(width=0.8, position="stack") +
      coord_flip() +
      theme_void() +
      scale_fill_manual(values=setNames(figureColors(colors), names(colors)), guide=FALSE) +
      theme(legend.position = "bottom",
            legend.direction="horizontal",
            strip.text=element_blank(),
            plot.margin = unit(c(0,0.1,2,0), "cm")) +
      ggtitle(variable)
  })
  
  statistics <- do.call(rbind, lapply(unique(clinicalDataFractions$variable), function(variable) {
    test <- do.call(rbind, lapply(groups, function(group) {
      sapply(unique(clinicalDataFractions[clinicalDataFractions$variable == variable, "option"]), function(option) {
        clinicalDataFractions[clinicalDataFractions$variable == variable &
                              clinicalDataFractions$group == group &
                              clinicalDataFractions$option == option, "count"]
      })
    })) %>% fisher.test()

    testWithoutUnknown <- do.call(rbind, lapply(groups, function(group) {
      sapply(setdiff(unique(clinicalDataFractions[clinicalDataFractions$variable == variable, "option"]), "unknown"),
             function(option) {
        clinicalDataFractions[clinicalDataFractions$variable == variable &
                              clinicalDataFractions$group == group &
                              clinicalDataFractions$option == option, "count"]
      })
    })) %>% fisher.test()
    
    df <- data.frame(
      variable = variable,
      p = test$p.value,
      sign = pToAsterisks(test$p.value),
      method = test$method,
      pWithoutUnknown = testWithoutUnknown$p.value,
      signWithoutUnknown = pToAsterisks(testWithoutUnknown$p.value),
      methodWithoutUnknown = testWithoutUnknown$method
    )

    df <- cbind(df, sapply(groups, function(group) {
      options <- sapply(unique(clinicalDataFractions[clinicalDataFractions$variable == variable, "option"]),
                        function(option) {
        count <- clinicalDataFractions[clinicalDataFractions$variable == variable &
                                       clinicalDataFractions$group == group &
                                       clinicalDataFractions$option == option, "count"]
        return(paste0(option, " (", count, ")"))
      })
      return(paste(options, collapse = ", "))
    }, simplify = FALSE, USE.NAMES = TRUE))

    return(df)
  }))

  list(
    plot = arrangeGrob(grobs=plots, nrow=1),
    statistics = statistics
  )

}
