## This is an R script for the conversion of mzTab to a better readable tsv format.

# clear entire workspace
rm(list = ls())

input.file <- commandArgs(TRUE)[1]
output.file <- commandArgs(TRUE)[2]

# find start of the section
startSection <- function(file, section.identifier) {
  data <- file(file, "r")
  row = 0
  while (TRUE) {
    row = row + 1
    line = readLines(data, n=1)
    if (substr(line, 1, 3)==section.identifier) {
      break
    }
  }
  close(data)
  return (row)
}

# function describing how to collapse rows
# In the case of string columns (e.g. accessions), the row entries are first made unique and then written to a comma-separated string.
# In all other cases, the entry of the first row is returned.
collapseRows <- function(x) {
  if (is.character(x)) {
    x <- paste(unique(x[!is.na(x)]), collapse=",")
    if (x=="") {
      return(NA)
    }
    else {
      return (x)
    }
  }
  else {
    return (x[1])
  }
}

# count the occurrences of character c in string s
countOccurrences <- function(c,s) {
  s2 <- gsub(c,"",s)
  return (nchar(s) - nchar(s2))
}

# check that the protein accession is of the format *|*|*
# Note that NA returns TRUE.
checkAccessionFormat <- function(accession) {
  if (is.na(accession)) {
    return (TRUE)
  }
  n <- length(accession)
  count <- countOccurrences("[|]",accession)
  m <- length(which(count==2))
  return (n==m)
}

# Extracts the second entry from a string of the form *|*|*.
getAccession <- function(string) {
  if (is.na(string)) {
    return (NA)
  }
  return (unlist(strsplit(string, "[|]"))[2])
}

# Extracts the third entry from a string of the form *|*|*.
getGene <- function(string) {
  if (is.na(string)) {
    return (NA)
  }
  return (unlist(strsplit(string, "[|]"))[3])
}

# read the PSM section of an mzTab file
readMzTabPSM <- function(file) {
  
  # find start of the PSM section
  first.row <- startSection(file, "PSH")
  
  # read entire mzTab
  data <- read.table(file, sep="\t", skip=first.row-1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)
  
  # extract PSM data
  psm.data <- data[which(data[,1]=="PSM"),]
  psm.data$PSH <- NULL
  
  # In case the accession column is of the format *|*|*, we split this column into an accession and a gene column.
  if (all(sapply(psm.data$accession, checkAccessionFormat))) {
    psm.data$gene <- sapply(psm.data$accession, getGene)
    psm.data$accession <- sapply(psm.data$accession, getAccession)
  }
  
  # In the mzTab format, PSMs with multiple protein accessions are written to multiple rows.
  # Here we collapse these rows and separate multiple accessions/genes by comma.
  psm.data <- aggregate(psm.data, by=list(temp = psm.data$PSM_ID), FUN=collapseRows)
  psm.data$temp <- NULL
  
  return (psm.data)
}

psm.data <- readMzTabPSM(input.file)
write.table(psm.data, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
