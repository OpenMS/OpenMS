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

# get index in protein groups list containing protein x
getIndex <- function(x, members) {
	g <- gsub(x, "", members, fixed=TRUE)
	d <- nchar(members)-nchar(g)
	return (which(d>0)[1])
}

# returns first entry of a comma-separated list
firstEntry <- function(x) {
	list <- strsplit(as.character(x),",",fixed=TRUE)
	return (unlist(lapply(list, '[[', 1)))
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

# read the PRT section of an mzTab file
readMzTabPRT <-function(file) {
  
  # find start of the PRT section
  first.row <- startSection(file, "PRH")
  
  # read entire mzTab
  data <- read.table(file, sep="\t", skip=first.row-1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)
  
  # extract PRT data
  protein.data <- data[which(data[,1]=="PRT"),]
  protein.data$PRH <- NULL
  
  # In case the accession column is of the format *|*|*, we split this column into an accession and a gene column.
  if (all(sapply(protein.data$accession, checkAccessionFormat))) {
    protein.data$gene <- sapply(protein.data$accession, getGene)
    protein.data$accession <- sapply(protein.data$accession, getAccession)
  }
  
  # split into different types
  proteins <- protein.data[which(protein.data$opt_global_protein_group_type=="single_protein"),]
  protein.groups <- protein.data[which(protein.data$opt_global_protein_group_type=="protein_group"),]
  indistinguishable.groups <- protein.data[which(protein.data$opt_global_protein_group_type=="indistinguishable_group"),]
  
  if ((dim(protein.groups)[1] > 0) && (dim(indistinguishable.groups)[1] > 0)) {
    # match indistinguishable groups to protein groups
    group.index <- as.vector(sapply(firstEntry(indistinguishable.groups.members), getIndex, members=protein.groups.members))
    table <- data.frame(cbind(group.index, indistinguishable.groups.members))

    # merge information from the protein list
    colnames(table) <- c("protein group","accessions")
    table$accession <- firstEntry(table$accessions)
    table <- merge(table, proteins, by="accession")
    table$accession <- NULL

    # order table by protein.group
    table$"protein group" <- as.integer(table$"protein group")
    table <- table[order(table$"protein group"),]
  }
  else {
    table <- proteins
    colnames(table) <- gsub("accession","accessions", colnames(table))
  }
  
  return (table)
}

protein.data <- readMzTabPRT(input.file)
write.table(protein.data, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
