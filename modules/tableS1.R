
#######################################################################################################################
# Table S1: Protein list
#######################################################################################################################

proteinList <- proteinData
proteinList$uniprotIDs <- ifelse(proteinList$uniprotID2 == "",
	                             proteinList$uniprotID1,
	                             paste0(proteinList$uniprotID1, "/", proteinList$uniprotID2))
proteinList$excluded <- ifelse(proteinList$proteinName %in% HD$proteinData$proteinName, "no", "yes")

formattedTable <- data.frame(Protein = proteinList$proteinName,
                             Gene    = proteinList$geneName,
                             Uniprot = proteinList$uniprotIDs,
                             Excl    = proteinList$excluded)
formattedTable <- formattedTable[order(formattedTable$Protein),]

formattedTable %>% write.table(file = file.path(outputDir, "tables", "Table S1.csv"),
	                           sep = ";", dec = ",", row.names = FALSE)


rm(proteinList, formattedTable)