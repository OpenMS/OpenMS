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

# count the occurrences of character c in string s
countOccurrences <- function(char,s) {
	s2 <- gsub(char,"",s)
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

# read the PEP section of an mzTab file
readMzTabPEP <- function(file) {
  
  # find start of the PEP section
  first.row <- startSection(file, "PEH")
  
  # read entire mzTab
  data <- read.table(file, sep="\t", skip=first.row-1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)
  
  # extract PEP data
  peptide.data <- data[which(data[,1]=="PEP"),]
  peptide.data$PEH <- NULL
  
  # In case the accession column is of the format *|*|*, we split this column into an accession and a gene column.
  if (all(sapply(peptide.data$accession, checkAccessionFormat))) {
    peptide.data$gene <- sapply(peptide.data$accession, getGene)
    peptide.data$accession <- sapply(peptide.data$accession, getAccession)
  }
  
  return (peptide.data)
}

peptide.data <- readMzTabPEP(input.file)
write.table(peptide.data, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
