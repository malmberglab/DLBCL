
plotScoreSurvival <- function(dataset, survivalType = "OS", scoreName, scoreGroup, query = NULL, colors = NULL,
                              onlyKM = FALSE) {
  
  df <- data.frame(months = dataset$survival[[survivalType]]$days / 365 * 12,
                   status = dataset$survival[[survivalType]]$status,
                   score  = dataset$score[[scoreName]][[scoreGroup]])
  
  if (!is.null(query)) {
    df <- df[query, ]
  }

  timeBreaks  <- if (max(df$months) > 60) 24 else 12

  plots <- ggsurvplot(survfit(Surv(months, status) ~ score, data=df), data = df, pval = TRUE, conf.int = FALSE,
                      risk.table = TRUE, break.time.by = timeBreaks)
  
  bothPlots <- sapply(c("plot", "table"), function(plotName) {
    plots[[plotName]] +
      xlab("Months") +
      scale_x_continuous(breaks = (0:5) * timeBreaks) +
      theme_classic() +
      theme(axis.text = element_text(size = 16, color = "black"),
            axis.title = element_text(size = 16, color = "black"),
            legend.position = "top",
            panel.grid.major = element_line(size = 0.5, color = "gray90"))
  }, simplify = FALSE, USE.NAMES = TRUE)

  bothPlots$plot <- bothPlots$plot + scale_y_continuous(breaks = (0:5) / 5, name = survivalType)
  bothPlots$table <- bothPlots$table + scale_y_discrete(labels = format((0:19) / 2, nsmall = 1))

  if (!is.null(colors)) {
    bothPlots$plot <- bothPlots$plot + scale_color_manual(values = colors, name = NULL)
  }

  if (onlyKM) return(bothPlots$plot)

  arrangeGrob(grobs = bothPlots, ncol = 1)
}


createVariablesForCox <- function(dataset, scoreName, scoreGroup) {
  list(
    score = dataset$score[[scoreName]]$numeric,
    IPI = dataset$metaData$IPI / 5,
    Bsymptoms = dataset$metaData$stage_AB,
    stageIIIorIV = sapply(dataset$metaData$stage, function(x) {
      if (is.na(x)) NA else if (x %in% c("I", "II")) "I/II" else if (x %in% c("III", "IV")) "III/IV" else NA
    }),
    gender = dataset$metaData$gender,
    subtype5 = sapply(dataset$metaData$subtype5, function(x) { if (x == "") NA else x })
  )
}


# Used by multivariateScoreAnalysis():
groupIPIhighLow <- function(scores, factor = TRUE) {
  res <- sapply(as.numeric(scores), function(score) {
    if (is.na(score)) {
      NA
    } else if (score > 2) {
      "high"
    } else {
      "low"
    }
  })
  
  if (factor) {
    factor(res, levels = c("low", "high"))
  } else {
    res
  }
}

convertAnnArborStageToInteger <- function(x) {
  as.numeric(sapply(x, function(y) {
    if (is.na(y)) { return(NA) }
    else if (y == "I") { return(1) }
    else if (y == "II") { return(2) }
    else if (y == "III") { return(3) }
    else if (y == "IV") { return(4) }
    else { return(NA) }
  }))
}


prepareDataForCoxRegression <- function(dataset, survivalType, variables) {
  variablesWithSurvival <- c(variables, list(
    days = dataset$survival[[survivalType]]$days,
    status = dataset$survival[[survivalType]]$status
  ))
  
  do.call(data.frame, sapply(names(variablesWithSurvival), function(variable) {
    sapply(variablesWithSurvival[[variable]], function(x) {
      if (is.na(x) || x == "") NA else x
    })
  }, simplify = FALSE, USE.NAMES = TRUE))
}


univariateCox <- function(dataset, survivalType = "OS", variables) {
  df <- prepareDataForCoxRegression(dataset, survivalType, variables)

  univariateSummary <- do.call(rbind, lapply(names(variables), function(variable) {
    coxreg <- summary(coxph(as.formula(paste("Surv(days, status) ~", variable)), data = df))

    if (variable == row.names(coxreg$coefficients)[1]) {
      # This is a continuous variable
      levels <- NA
      reference <- NA
    } else {
      levels <- gsub(paste0("^", variable), "", row.names(coxreg$coefficients))
      reference <- unique(df[, variable])[unique(df[, variable]) %in% c(NA, levels) == FALSE]
    }
    
    
    ret <- data.frame(
      variable    = variable,
      level       = levels,
      reference   = reference,
      count       = if (all(is.na(levels))) { sum(!is.na(df[,variable])) } else {
                      sapply(levels, function(x) { sum(df[,variable] == x, na.rm=TRUE) }) },
      p           = coxreg$coefficients[, 5],
      sign        = pToAsterisks(coxreg$coefficients[, 5]),
      logrankP    = coxreg$sctest["pvalue"],
      logrankSign = pToAsterisks(coxreg$sctest["pvalue"]),
      coef        = coxreg$coefficients[, 1],
      HR          = coxreg$coefficients[, 2],
      lower95     = coxreg$conf.int[, 3],
      upper95     = coxreg$conf.int[, 4],
      isReference = FALSE
    )

    if (!is.na(reference)) {
      ret <- rbind(data.frame(
        variable    = variable,
        level       = reference,
        reference   = NA,
        count       = sum(df[, variable] == reference, na.rm = TRUE),
        p           = NA,
        sign        = NA,
        logrankP    = NA,
        logrankSign = NA,
        coef        = 0,
        HR          = 1,
        lower95     = 1,
        upper95     = 1,
        isReference = TRUE
      ), ret)
    }

    ret
  }))
  row.names(univariateSummary) <- NULL

  return(list(
    data = df,
    summary = univariateSummary[!univariateSummary$isReference, names(univariateSummary) != "isReference"],
    summaryWithReference = univariateSummary,
    survivalType = survivalType
  ))
}



plotUnivariateCox <- function(univariateCoxResult, customPlotLimits=NULL) {
  plotCox(univariateCoxResult$summaryWithReference, paste0(univariateCoxResult$survivalType, ": univariate Cox"),
          customPlotLimits)
}

plotUnivariateCoxWithMultivariate <- function(univariateCoxResult, pValueCutoff = 0.1, customPLotLimits = NULL) {
  multivariateVariables <- unique(univariateCoxResult$summary[univariateCoxResult$summary$logrankP <= pValueCutoff,
                                                              "variable"])
  multivariateCoxResult <- multivariateCox(univariateCoxResult$data, univariateCoxResult$survivalType,
                                           multivariateVariables)

  arrangeGrob(grobs = list(plotUnivariateCox(univariateCoxResult),
                           plotMultivariateCox(multivariateCoxResult)), ncol = 2)
}

multivariateCox <- function(dataset, survivalType = "OS", variables) {
  multivariateCoxWithData(prepareDataForCoxRegression(dataset, survivalType, variables), survivalType)
}

multivariateCoxWithData <- function(data, survivalType = "OS") {
  df <- data
  variables <- names(data)[names(data) %in% c("days", "status") == FALSE]

  coxph <- coxph(as.formula(paste("Surv(days, status) ~", paste(variables, collapse=" + "))), data = df)
  coefficients <- cbind(summary(coxph)$coefficients, summary(coxph)$conf.int)
  statistics <- list(
    n                   = summary(coxph)$n,
    nevent              = summary(coxph)$nevent,
    concordance         = as.numeric(summary(coxph)$concordance["C"]),
    likelihoodRatioTest = summary(coxph)$logtest,
    waldTest            = summary(coxph)$waldtest,
    logrankTest         = summary(coxph)$sctest
  )

  multivariateSummary <- do.call(rbind, lapply(variables, function(variable) {
    if (variable %in% row.names(coefficients)) {
      # This is a continuous variable
      data.frame(
        variable    = variable,
        level       = NA,
        reference   = NA,
        count       = sum(!is.na(df[,variable])),
        p           = coefficients[variable, "Pr(>|z|)"],
        sign        = pToAsterisks(coefficients[variable, "Pr(>|z|)"]),
        coef        = coefficients[variable, "coef"],
        HR          = coefficients[variable, "exp(coef)"],
        lower95     = coefficients[variable, "lower .95"],
        upper95     = coefficients[variable, "upper .95"],
        isReference = FALSE
      )
    } else {
      levels <- as.character(row.names(coefficients))
      levels <- levels[startsWith(levels, variable)]
      levels <- gsub(paste0("^", variable), "", levels)
      reference <- as.character(unique(df[,variable])[unique(df[,variable]) %in% c(NA, levels) == FALSE])

      do.call(rbind, lapply(c(reference, as.character(levels)), function(level) {
        row <- paste0(variable, level)
        if (row %in% row.names(coefficients)) {
          data.frame(
            variable    = variable,
            level       = level,
            reference   = reference,
            count       = sum(df[,variable] == level, na.rm = TRUE),
            p           = coefficients[row, "Pr(>|z|)"],
            sign        = pToAsterisks(coefficients[row, "Pr(>|z|)"]),
            coef        = coefficients[row, "coef"],
            HR          = coefficients[row, "exp(coef)"],
            lower95     = coefficients[row, "lower .95"],
            upper95     = coefficients[row, "upper .95"],
            isReference = FALSE
          )
        } else {
          if (level != reference) stop("Non-reference level not found in multivariate Cox: ", level,
                                       " (ref: ", reference, ", variable: ", variable, ")")
          data.frame(
            variable    = variable,
            level       = level,
            reference   = reference,
            count       = sum(df[,variable] == level, na.rm = TRUE),
            p           = NA,
            sign        = NA,
            coef        = 0,
            HR          = 1,
            lower95     = NA,
            upper95     = NA,
            isReference = TRUE
          )
        }
      }))
    }
  }))
  row.names(multivariateSummary) <- NULL

  return(list(
    data = df,
    coxph = coxph,
    statistics = statistics,
    summary = multivariateSummary[!multivariateSummary$isReference, names(multivariateSummary) != "isReference"],
    summaryWithReference = multivariateSummary,
    survivalType = survivalType
  ))
}

multivariateCoxFromUnivariate <- function(univariateCoxResult, pValueCutoff=0.1) {
  multivariateVariables <- c("days", "status", unique(univariateCoxResult$summary[univariateCoxResult$summary$p <=
                                                                                    pValueCutoff, "variable"]))
  multivariateCoxResult <- multivariateCoxWithData(univariateCoxResult$data[,multivariateVariables],
                                                   univariateCoxResult$survivalType)
  multivariateSummary <- univariateCoxResult$summaryWithReference[, c("variable", "level", "reference", "isReference")]
  columnsToFill <- c("count", "p", "sign", "coef", "HR", "lower95", "upper95")
  multivariateSummary[, columnsToFill] <- NA
  for (row in 1:nrow(multivariateSummary)) {
    if (is.na(multivariateSummary[row, "level"])) {
      oldRow <- which(multivariateCoxResult$summaryWithReference$variable == multivariateSummary[row, "variable"] &
                      is.na(multivariateCoxResult$summaryWithReference$level))
    } else {
      oldRow <- which(multivariateCoxResult$summaryWithReference$variable == multivariateSummary[row, "variable"] &
                      multivariateCoxResult$summaryWithReference$level == multivariateSummary[row, "level"])
    }
    
    if (length(oldRow) != 0) {
      multivariateSummary[row, columnsToFill] <- multivariateCoxResult$summaryWithReference[oldRow, columnsToFill]
    }
  }

  multivariateCoxResult$summaryWithReference <- multivariateSummary
  multivariateCoxResult$summary <- multivariateSummary[!multivariateSummary$isReference,
                                                       names(multivariateSummary) != "isReference"]

  multivariateCoxResult
}

plotMultivariateCox <- function(multivariateCoxResult, customPlotLimits=NULL) {
  plotCox(multivariateCoxResult$summaryWithReference, paste0(multivariateCoxResult$survivalType, ": multivariate Cox"),
          customPlotLimits)
}

plotCox <- function(summary, title="", customPlotLimits=NULL) {
  logBase <- if (!is.null(customPlotLimits) && any(abs(customPlotLimits) >= 100)) 10 else 2

  # Add some texts. Need to do this before we start messing with lower95/upper95 for the arrows
  summary$HRText <- ifelse(summary$isReference, "(ref)", paste0(round(summary$HR, 2), " (", round(summary$lower95, 2),
                                                                "-", round(summary$upper95, 2), ")"))
  summary$pValueText <- ifelse(summary$isReference, "", pValueFormat(summary$p))
  summary$countText <- paste0("n = ", summary$count)

  summary[is.na(summary$HR), c("HRText", "pValueText", "countText")] <- ""


  summary$upperCut <- NA
  summary$lowerCut <- NA

  if (!is.null(customPlotLimits)) {
    xlims <- customPlotLimits

    lowerCut <- !is.na(summary$HR) & !is.na(summary$lower95) & summary$HR > xlims[1] & summary$lower95 < xlims[1]
    summary[lowerCut, "lowerCut"] <- xlims[1]
    summary[lowerCut, "lower95"] <- summary[lowerCut, "HR"]

    upperCut <- !is.na(summary$HR) & !is.na(summary$upper95) & summary$HR < xlims[2] & summary$upper95 > xlims[2]
    summary[upperCut, "upperCut"] <- xlims[2]
    summary[upperCut, "upper95"] <- summary[upperCut, "HR"]
  } else {
    xlims <- c(min(summary$lower95, na.rm=TRUE), max(summary$upper95, na.rm=TRUE))
  }

  originalXlims <- xlims
  xlims <- c((logBase^(log(xlims[1] * 0.8, base = logBase))),
             (logBase^(log(xlims[2] * 1.2, base = logBase) + log(10^6, base = logBase))))

  breaks <- logBase^((floor(log(originalXlims[1], base=logBase)/(if (logBase < 10) 2 else 1)):
                        ceiling(log(originalXlims[2], base=logBase)/2))*(if (logBase < 10) 2 else 1))
  rows <- paste0(summary$variable, ": ", summary$level)

  ggplot(summary, aes(y = factor(paste0(variable, ": ", level), levels = rev(rows)), x = HR)) +
    geom_rect(xmin = log(xlims[1], base = logBase), xmax = log(xlims[2], base = logBase), ymin=seq_along(rows) - 0.5,
              ymax = seq_along(rows) + 0.5, fill = rep(c("#ffffff", "#f0f0f0"), length.out = length(rows))) +
    geom_vline(xintercept = 1, color = "black", linetype = "dotted") +
    geom_point(shape = 15, size = 3) +

    geom_segment(aes(x = HR, xend = lowerCut, yend = factor(paste0(variable, ": ", level))),
                 arrow = arrow(length = unit(0.15, "cm"), ends = "last"), na.rm = TRUE) +
    geom_segment(aes(x = HR, xend = upperCut, yend = factor(paste0(variable, ": ", level))),
                 arrow = arrow(length = unit(0.15, "cm"), ends = "last"), na.rm = TRUE) +

    geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.35, na.rm = TRUE) +
    scale_x_continuous(breaks = breaks, labels = as.character(breaks),
                       trans = paste0("log", logBase), expand = c(0,0), limits = xlims) +
    geom_text(aes(label = countText, x = (logBase^(log(originalXlims[2] * 1.2, base = logBase) + 0))),
              hjust = 0, size = 4) +
    geom_text(aes(label = HRText, x = (logBase^(log(originalXlims[2] * 1.2, base = logBase) +
                                                  log(32, base = logBase)))), hjust = 0, size = 4) +
    geom_text(aes(label = pValueText, x = (logBase^(log(originalXlims[2] * 1.2, base = logBase) + 
                                                      log(16000, base = logBase)))), hjust = 0, size = 4) +
    ylab(NULL) + xlab("HR (95% conf.int.)") +
    theme_classic() +
    theme(
      axis.text    = element_text(colour = "black", size = 14),
      axis.text.y  = element_text(hjust=0),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank()
    ) +
    ggtitle(title)
}



