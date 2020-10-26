pToAsterisks <- function(p) {
  return(sapply(p, function(x) {
    if (is.na(x)) { "" }
    else if (x <= 0.0001) { "****" }
    else if (x <= 0.001) { "***" }
    else if (x <= 0.01) { "**" }
    else if (x <= 0.05) { "*" }
    else { "ns." }
  }))
}

pValueFormat <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) {
      return(NA)
    } else {
      if (x < 0.0001) "p < 0.0001" else paste("p =", format(round(x, 4), scientific=FALSE))
    }
  })    
}

formatProteinName <- function(proteins, to="proteinName", data=proteinData) {
  sapply(proteins, function(protein) {
    ret <- data[protein, to]
    if (length(ret) < 1) {
      warning("Did not find ", to, " for protein: ", protein)
      return(protein)
    } else if (length(ret) > 1) {
      warning("Multiple ", to, " for protein: ", protein)
      return(ret[1])
    } else {
      return(ret)
    }
  })
}
