# addScore()
# add a score to a data set
# if quantile is set to NA, the score will be the sum of quantiles

addScore <- function(dataset, proteins=NULL, scores=NULL, quantile, name) {
  if (!is.null(proteins)) {
    scores <- createScores(dataset, proteins, quantile)
  }

  if ("score" %in% names(dataset) == FALSE) dataset$score <- list()

  dataset$score[[name]] <- list(
    numeric = scores,
    groupedByMedian = groupScoresByMedian(scores),
    groupedByTertiles = groupScoresByTertiles(scores),
    groupedBySecondTertile = groupScoresBySecondTertile(scores),
    quantile = quantile,
    proteins = proteins
  )

  return(dataset)
}


groupScoresByMedian <- function(scores) {
  factor(sapply(scores, function(x) {8
    if (x < quantile(scores, 1/2)) {
      "low"
    } else {
      "high"
    }
  }), levels=c("low", "high"))
}


groupScoresByTertiles <- function(scores) {
  factor(sapply(scores, function(x) {
    if (x < quantile(scores, 1/3)) {
      "low"
    } else if (x < quantile(scores, 2/3)) {
      "mid"
    } else {
      "high"
    }
  }), levels=c("low", "mid", "high"))
}


groupScoresBySecondTertile <- function(scores) {
  factor(sapply(scores, function(x) {
    if (x < quantile(scores, 2/3)) {
      "low/mid"
    } else {
      "high"
    }
  }), levels=c("low/mid", "high"))
}


createScores <- function(dataset, proteins, quantile) {
  dataset$NPXData <- useNormalizedDataIfAvailable(dataset)

  if (length(proteins) == 1) {
    if (is.na(quantile)) {
      ecdf(dataset$NPXData[,proteins])(dataset$NPXData[,proteins])
    } else {
      ifelse(dataset$NPXData[,proteins] > quantile(dataset$NPXData[,proteins], quantile), 1, 0)
    }
  } else {
    rowSums(apply(dataset$NPXData[,proteins], 2, function(x) {
      if (is.na(quantile)) {
        ecdf(x)(x)
      } else {
        x > quantile(x, quantile)
      }
    }))
  }
}

