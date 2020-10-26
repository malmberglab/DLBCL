
#######################################################################################################################
# Table 2
#######################################################################################################################

formatRegressionTable <- function(obj) {
  data.frame(
    variable = ifelse(is.na(obj$reference), obj$variable, ""),
    level = ifelse(is.na(obj$level), "(cont.)", obj$level),
    HR = ifelse(obj$isReference, "(ref)",
                paste0(round(obj$HR, 2), " (", round(obj$lower95, 2), ", ", round(obj$upper95, 2), ")")),
    p = sapply(obj$p, function(p){ if (is.na(p)) "" else if (p < 0.0001) "< 0.0001" else round(p, 4) })
  )
}

variablesForRegression <- c("score", "IPI", "Bsymptoms", "gender", "stageIIIorIV", "subtype5")

osloVariableDataForUnivariateCox <- createVariablesForCox(osloMerged, "higherInPatientsQuantiles",
                                                          "groupedBySecondTertile")[variablesForRegression]

univariateCoxOslo <- univariateCox(osloMerged, "OS", osloVariableDataForUnivariateCox)
univariateCoxOslo$summaryWithReference %>% formatRegressionTable() %>%
  write.table(file = file.path(outputDir, "tables", "Table 2 Oslo univar.csv"),
              sep = ";", dec = ",", row.names = FALSE)

univariateCoxOslo %>%
  multivariateCoxFromUnivariate(pValueCutoff = 1) %>%
  { .$summaryWithReference } %>%
  formatRegressionTable() %>%
  write.table(file = file.path(outputDir, "tables", "Table 2 Oslo multivar.csv"),
              sep = ";", dec = ",", row.names = FALSE)


WUVariableDataForUnivariateCox <- createVariablesForCox(WU, "higherInPatientsQuantiles",
                                                        "groupedBySecondTertile")[variablesForRegression]

univariateCoxWU <- univariateCox(WU, "OS", WUVariableDataForUnivariateCox)
univariateCoxWU$summaryWithReference %>%
  formatRegressionTable() %>%
  write.table(file = file.path(outputDir, "tables", "Table 2 St Louis univar.csv"),
              sep = ";", dec = ",", row.names = FALSE)

univariateCoxWU %>%
  multivariateCoxFromUnivariate(pValueCutoff = 1) %>%
  { .$summaryWithReference } %>%
  formatRegressionTable() %>%
  write.table(file = file.path(outputDir, "tables", "Table 2 St Louis multivar.csv"),
              sep = ";", dec = ",", row.names = FALSE)




rm(formatRegressionTable, variablesForRegression, osloVariableDataForUnivariateCox, univariateCoxOslo,
   WUVariableDataForUnivariateCox, univariateCoxWU)