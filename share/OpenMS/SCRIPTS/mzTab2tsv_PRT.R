## This is an R script for the conversion of mzTab to a better readable tsv format

# options
options(digits=10)

input.file <- commandArgs(TRUE)[1]
output.file <- commandArgs(TRUE)[2]

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

# count the occurences of character c in string s
countOccurrences <- function(char,s) {
	s2 <- gsub(char,"",s)
	return (nchar(s) - nchar(s2))
}

# check that all protein accessions are of the format *|*|*
checkAccessionFormat <- function(accessions) {
	n <- length(accessions)
	count <- countOccurrences("[|]",accessions)
	m <- length(which(count>=2))
	return (n==m)
}

getAccessions <- function(string) {
	accessions <- strsplit(string,",")
	accessions <- lapply(accessions,strsplit,"[|]")
	accessions <- data.frame(matrix(unlist(accessions),ncol=3,byrow=TRUE))
	return (paste(as.character(accessions[,2]), collapse=","))
}

getGenes <- function(string) {
	accessions <- strsplit(string,",")
	accessions <- lapply(accessions,strsplit,"[|]")
	accessions <- data.frame(matrix(unlist(accessions),ncol=3,byrow=TRUE))
	return (paste(as.character(accessions[,3]), collapse=","))
}

# read the PRT section of an mzTab file
readMzTabPRT <- function(file) {
  
  # read entire mzTab
  no.col <- max(count.fields(file, sep = "\t", quote=""))
  data <- read.table(file, sep="\t", fill=TRUE, quote="", col.names=1:no.col)
  
  # extract protein data
  protein.data <- data[which(data[,1]=="PRT"),]
  colnames(protein.data) <- unlist(data[which(data[,1]=="PRH")[1],])
  protein.data$PRH <- NULL
  columns.to.keep <- which(colnames(protein.data)!="")
  protein.data <- protein.data[,columns.to.keep]
  
  proteins <- protein.data[which(protein.data$opt_global_protein_group_type=="single_protein"),]
  protein.groups <- protein.data[which(protein.data$opt_global_protein_group_type=="protein_group"),]
  indistinguishable.groups <- protein.data[which(protein.data$opt_global_protein_group_type=="indistinguishable_group"),]

  protein.groups.members <- as.character(protein.groups$ambiguity_members)
  
  indistinguishable.groups.members <- as.character(indistinguishable.groups$ambiguity_members)
  
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
  } else {
    table <- proteins
    colnames(table) <- gsub("accession","accessions", colnames(table))
  }
  
  # simplify accessions (in case they are of the format *|*|* )
  table$accessions <- as.character(table$accessions)
  if (checkAccessionFormat(table$accessions)) {
    x <- unlist(lapply(table$accessions,getAccessions))
    y <- unlist(lapply(table$accessions,getGenes))
    table$accessions <- x
    table$gene <- y
  }

  return (protein.groups)
}

table <- readMzTabPRT(input.file)
write.table(table, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
