#
# A Systemic Protein Deviation Score Linked to PD-1+ CD8+ T-cell Expansion That
# Predicts Overall Survival in Diffuse Large B-cell Lymphoma
#
# Ask, Tschan-Plessl et al.
#


# Directory of this file
dir <- ""

# Directory where data is located
dataDir <- file.path(dir, "data")

# Directory where figures and statistics should be saved
outputDir <- file.path(dir, "output")


###


if (dir.exists(outputDir)) {
  stop("Output directory already exists!")
} else {
  dir.create(outputDir)
  dir.create(file.path(outputDir, "figures"))
  dir.create(file.path(outputDir, "tables"))
  dir.create(file.path(outputDir, "statistics"))
  
  if (!dir.exists(file.path(dir, "cache"))) dir.create(file.path(dir, "cache"))
  
  source(file.path(dir, "load.R"))
  
  # Generate figures and statistics
  for (file in list.files(file.path(dir, "modules"), full = TRUE)) {
    cat("\n== ", file, " ==\n\n", sep = "")
    source(file)
  }
}

  