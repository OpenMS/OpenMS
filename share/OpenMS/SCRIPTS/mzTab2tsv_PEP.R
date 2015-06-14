## This is an R script for the conversion of mzTab to a better readable tsv format

# options
options(digits=10)

input.file <- commandArgs(TRUE)[1]
output.file <- commandArgs(TRUE)[2]

# read entire mzTab
no.col <- max(count.fields(input.file, sep = "\t", quote=""))
data <- read.table(input.file,sep="\t",fill=TRUE, quote="", col.names=1:no.col)

# extract peptide data
peptide.data <- data[which(data[,1]=="PEP"),]
colnames(peptide.data) <- unlist(data[which(data[,1]=="PEH")[1],])
peptide.data$PEH <- NULL

countOccurrences <- function(char,s) {
	s2 <- gsub(char,"",s)
	return (nchar(s) - nchar(s2))
}

checkAccessionFormat <- function(accessions) {
	n <- length(accessions)
	count <- countOccurrences("[|]",accessions)
	m <- length(which(count==2))
	return (n==m)
}

# simplify accession (in case it is of the format *|*|* )
peptide.data$accession <- as.character(peptide.data$accession)
if (checkAccessionFormat(peptide.data$accession)) {
	list <- strsplit(peptide.data$accession,"[|]")
	peptide.data$accession <- unlist(lapply(list, '[[', 2))
	peptide.data$gene <- unlist(lapply(list, '[[', 3))
}

write.table(peptide.data, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)