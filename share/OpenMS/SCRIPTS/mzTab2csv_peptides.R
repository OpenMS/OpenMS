## This is an R script reformats mzTab files (*.tsv) to excel files (*.csv)

#source("http://bioconductor.org/biocLite.R")
#biocLite("MSnbase")
library(MSnbase)

input.file<-commandArgs(TRUE)[1]
output.file<-commandArgs(TRUE)[2]

# MSnSet
x <- readMzTabData(input.file, what = "PEP")

# data frame
peptide.data <- fData(x)
peptide.data[] <- lapply(peptide.data, as.character)

# column names
columns <- colnames(peptide.data)

unique.stripped.sequences <- unique(peptide.data$sequence)
unique.modified.sequences <- unique(peptide.data$opt_global_modified_sequence)

# find best (lowest PEP score) modified sequences
n <- length(unique.modified.sequences)
indices <- c()
for (i in 1:n) {
	j <- which(peptide.data$"opt_global_modified_sequence" == unique.modified.sequences[i])
	PEPs <- peptide.data[j,]$"opt_psm_Posterior Error Probability_score"
	k <- j[which(PEPs == min(PEPs))]
	k <- k[1]
	
	indices <- c(indices, k)
}

# simplify accession
list <- strsplit(peptide.data$accession,"[|]")
peptide.data$accession <- unlist(lapply(list, '[[', 2))

# add fold change
peptide.data$fold_change <- log(as.numeric(peptide.data$"peptide_abundance_study_variable[2]")/as.numeric(peptide.data$"peptide_abundance_study_variable[1]"), base =2)

# prepare output
output.columns <- c(1,24,2,3,28,25,12,10,13,15,18,29)
output <- peptide.data[indices,output.columns]
colnames(output) <- c("sequence_stripped","sequence_modified","accession","proteotypic","target_decoy","posterior_error_probability","charge","retention_time","mass_to_charge","abundance_light","abundance_heavy","fold_change")
write.table(output, output.file, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
