## This is an R script for the conversion of mzTab to a better readable tsv format

# options
options(digits=10)

input.file <- commandArgs(TRUE)[1]
output.file <- commandArgs(TRUE)[2]

# count the occurences of character c in string s
countOccurrences <- function(c,s) {
  s2 <- gsub(c,"",s)
  return (nchar(s) - nchar(s2))
}

# check that all protein accessions are of the format *|*|* 
checkAccessionFormat <- function(accessions) {
  n <- length(accessions)
  count <- countOccurrences("[|]",accessions)
  m <- length(which(count==2))
  return (n==m)
}

# collapse rows
# (In mzTab, PSMs with multiple protein accessions are reported in multiple rows. This function collapses them to a single row.)
collapseRows <- function(psm.data) {
  
  # generate index vector idx
  tmp.psm.id <- 0
  idx <- c()
  accessions.tmp <- c()
  accessions.strings <- c()
  for (i in 1:length(psm.data$PSM_ID)) {
    
    if (psm.data$PSM_ID[i] == tmp.psm.id) {
      if (length(accessions.tmp) > 0) {
        accessions.strings <- c(accessions.strings, toString(accessions.tmp, sep=','))
        accessions.tmp <- c()
      }
      
      idx <- c(idx,i)
      tmp.psm.id <- tmp.psm.id + 1
    }
    
    accessions.tmp <- c(accessions.tmp, psm.data$accession[i])
  }
  accessions.strings <- c(accessions.strings, toString(accessions.tmp, sep=','))
  
  psm.data <- psm.data[idx,]
  psm.data$accession <- accessions.strings
  
  return (psm.data)
}

# read the PSM section of an mzTab file
readMzTabPSM <- function(file) {
  
  # read entire mzTab
  no.col <- max(count.fields(file, sep = "\t", quote=""))
  data <- read.table(file, sep="\t", fill=TRUE, quote="", col.names=1:no.col)
  
  # extract PSM data
  psm.data <- data[which(data[,1]=="PSM"),]
  colnames(psm.data) <- unlist(data[which(data[,1]=="PSH")[1],])
  psm.data$PSH <- NULL
  
  # simplify accession (in case it is of the format *|*|* )
  psm.data$accession <- as.character(psm.data$accession)
  if (checkAccessionFormat(psm.data$accession)) {
    list <- strsplit(psm.data$accession, "[|]")
    psm.data$accession <- unlist(lapply(list, '[[', 2))
    psm.data$gene <- unlist(lapply(list, '[[', 3))
  }
  
  psm.data <- collapseRows(psm.data)
  
  return (psm.data)
}

psm.data <- readMzTabPSM(input.file)
write.table(psm.data, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
