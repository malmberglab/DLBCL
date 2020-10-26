
fetchGOKeywordsForProteins <- function(proteins) {
  sapply(proteins, function(protein) {
    terms <- c()
    for (i in 1:2) {
      uniprotID <- proteinData[protein, paste0("uniprotID", i)]
      if (uniprotID != "") {
        uniprotXML <- read_xml(paste0("http://www.uniprot.org/uniprot/", uniprotID, ".xml"))
        keywordsEntries <- uniprotXML %>% xml_find_all(".//d1:entry/d1:keyword")
        terms <- c(terms,
                   keywordsEntries %>% { paste0(xml_attr(., "id"), "//", xml_text(.)) })

      }
    }
    return(unique(terms))
  }, simplify=FALSE, USE.NAMES=TRUE)
}

fetchUniprotKeywordsForParent <- function(parent) {
  read_xml(paste0("https://www.uniprot.org/keywords/", parent, ".rdf")) %>%
    xml_find_all(".//skos:narrower") %>%
    xml_attr("resource") %>%
    sapply(function(x) { strsplit(x, "/")[[1]] %>% tail(n=1) }) %>%
    as.character() %>%
    str_pad(4, side = "left", pad = "0") %>%
    { paste0("KW-", .) }
}

filterGOKeywords <- function(keywordsForProteins, parent) {
  parentKeywords <- fetchUniprotKeywordsForParent(parent)
  sapply(names(keywordsForProteins), function(protein) {
    keywords <- keywordsForProteins[[protein]]
    keywords[substr(keywords, 1, 9) %in% paste0(parentKeywords, "//")]
  }, simplify = FALSE, USE.NAMES = TRUE)
}

formatGOKeywordNames <- function(keywordsForProteins) {
  sapply(names(keywordsForProteins), function(protein) {
    sapply(keywordsForProteins[[protein]], function(keyword) {
      strsplit(keyword, "//")[[1]][2]
    })  %>% as.character()
  }, simplify = FALSE, USE.NAMES = TRUE)
}

