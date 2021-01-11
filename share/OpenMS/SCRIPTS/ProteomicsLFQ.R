# example script to perform an MSstats analysis
Sys.setenv(LANG = "en")
library(MSstats)
library(dplyr)
library(tibble)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

MSstats_input <- args[1] 
mzTab_input <- args[2] 
mzTab_output <- args[3]

Sys.setenv(LANG = "en")

data <- read.csv(MSstats_input, sep=",", header=T, stringsAsFactors=T)

quant <- OpenMStoMSstatsFormat(data, removeProtein_with1Feature=T)
processed.quant <- dataProcess(quant, censoredInt = 'NA')

# generating these plots takes quite a while. Disable if needed.
dataProcessPlots(data=processed.quant, type="QCPlot",which.Protein="allonly")
dataProcessPlots(data=processed.quant, type="ProfilePlot",which.Protein="all")

# iPRG2015 study matrix
comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

############ also calculate missingness on condition level

# input: ProcessedData matrix of MSstats
# output: 
#   calculate fraction of na in condition (per protein)
# Groups:   PROTEIN [762]
#   PROTEIN                 `1`   `2`
#   <fct>                 <dbl> <dbl>
# 1 sp|A1ANS1|HTPG_PELPD   0    0.5  
# 2 sp|A2I7N3|SPA37_BOVIN  0    0.5  
# 3 sp|A2VDF0|FUCM_HUMAN   0    0.5  
# 4 sp|A6ND91|ASPD_HUMAN   0.5  0.5  
# 5 sp|A7E3W2|LG3BP_BOVIN  0.5  0.5  
# 6 sp|B8FGT4|ATPB_DESAA   0    0.5


getMissingInCondition <- function(processedData)
{
p <- processedData

# count number of samples per condition
n_samples = p %>% group_by(GROUP) %>% summarize(n_samples = length(unique((as.numeric(SUBJECT))))) 

p <- p %>% 
   filter(!is.na(INTENSITY)) %>% # remove rows with INTENSITY=NA
   select(PROTEIN, GROUP, SUBJECT) %>%
   distinct() %>% 
   group_by(PROTEIN, GROUP) %>% 
   summarize(non_na = n())  # count non-NA values for this protein and condition
   
p <- left_join(p, n_samples) %>% 
       mutate(missingInCondition = 1 - non_na/n_samples) # calculate fraction of missing values in condition

# create one column for every condition containing the missingness
p <- spread(data = p[,c("PROTEIN", "GROUP", "missingInCondition")], key = GROUP, value = missingInCondition)
return(p)
}
mic <- getMissingInCondition(processed.quant$ProcessedData)
#filtered.quant <- processed.quant
#filtered.quant$RunlevelData <- merge(x=processed.quant$RunlevelData, y=mic, by.y="PROTEIN", by.x="Protein")
#filtered.quant$RunlevelData[is.na(filtered.quant$RunlevelData)] <- 1 # set completely missing to 1.0 (had no matching entry in join and were set to NA)
#filtered.quant$ProcessedData <- merge(x=processed.quant$ProcessedData, y=mic, by="PROTEIN")
#filtered.quant$ProcessedData[is.na(filtered.quant$ProcessedData)] <- 1 # set completely missing to 1.0 (had no matching entry in join and were set to NA)

groupcomp <- groupComparison(contrast.matrix=comparison, data=processed.quant)
# for plotting, remove proteins with infinite fold change / p-value NA (e.g., those only present in one condition)
groupcomp$Volcano = groupcomp$ComparisonResult[!is.na(groupcomp$ComparisonResult$pvalue),]
groupComparisonPlots(data=groupcomp$Volcano, type="VolcanoPlot", width=12, height=12,dot.size = 2,ylimUp = 7)

# annotate how often the protein was quantified in each condition (NA values introduced by merge of completely missing are set to 1.0)
groupcomp$ComparisonResult <- merge(x=groupcomp$ComparisonResult, y=mic, by.x="Protein", by.y="PROTEIN")
commoncols <- intersect(colnames(mic), colnames(groupcomp$ComparisonResult))
groupcomp$ComparisonResult[, commoncols]<-groupcomp$ComparisonResult %>% select(commoncols) %>% mutate_all(funs(replace(., is.na(.), 1))) 
				 			 
#write comparison to CSV (one CSV per contrast)							 
writeComparisonToCSV <- function(DF) 
{
write.table(DF, file=paste0("comparison_",unique(DF$Label),".csv"), quote=FALSE, sep='\t', row.names = FALSE)
return(DF)
}
groupcomp$ComparisonResult %>% group_by(Label) %>% do(writeComparisonToCSV(as.data.frame(.)))  


################# MzTab
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

# find start of the mzTab section tables
MTD.first_row <- startSection(mzTab_input, "MTD")
PRT.first_row <- startSection(mzTab_input, "PRH")
PEP.first_row <- startSection(mzTab_input, "PEH")
PSM.first_row <- startSection(mzTab_input, "PSH")

# read entire mzTab and extract protein data
MTD <- read.table(mzTab_input, sep="\t", 
  skip=MTD.first_row-1,
  nrows=PRT.first_row - MTD.first_row - 1 -1, # one extra empty line
  fill=TRUE, 
  header=TRUE, 
  quote="", 
  na.strings=c("null","NA"), 
  stringsAsFactors=FALSE, 
  check.names=FALSE)  


PRT <- read.table(mzTab_input, sep="\t", 
  skip=PRT.first_row-1,
  nrows=PEP.first_row - PRT.first_row - 1 -1, # one extra empty line
  fill=TRUE, 
  header=TRUE, 
  quote="", 
  na.strings=c("null","NA"), 
  stringsAsFactors=FALSE, 
  check.names=FALSE)  

PEP <- read.table(mzTab_input, sep="\t", 
  skip=PEP.first_row-1,
  nrows=PSM.first_row - PEP.first_row - 1 - 1, # one extra empty line
  fill=TRUE, 
  header=TRUE, 
  quote="", 
  na.strings=c("null","NA"), 
  stringsAsFactors=FALSE, 
  check.names=FALSE)  

PSM <- read.table(mzTab_input, sep="\t", 
  skip=PSM.first_row-1,
  fill=TRUE, 
  header=TRUE, 
  quote="", 
  na.strings=c("null","NA"), 
  stringsAsFactors=FALSE, 
  check.names=FALSE)  

#### Insert quantification data from MSstats into PRT section
# first we create a run level protein table form MSstats output
# then we merge the values into the mzTab PRT table


# Input: MSstats RunLevelData
# Output: Run level quantification
# Create a run level protein table
#   PROTEIN                 `1`   `2`   `3`   `4`   `5`   `6`   `7`   `8`   `9`  `10`  `11`  `12`  `13`  `14`  `15`  `16`  `17`  `18`  `19`  `20`
#   <fct>                 <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 sp|A1ANS1|HTPG_PELPD   24.2  24.9  22.8  25.3  24.7  22.9  24.6  25.1  24.0  22.1  25.0  24.3  23.6  NA    NA    NA    NA    NA    NA    NA  
# 2 sp|A2I7N1|SPA35_BOVIN  22.9  23.6  22.4  23.8  23.4  NA    23.6  23.9  22.5  NA    23.7  23.5  22.5  22.5  23.0  23.0  22.6  22.2  22.1  22.8
getRunLevelQuant <- function(runLevelData)
{
runlevel.long <- tibble(RUN=as.numeric(runLevelData$RUN), PROTEIN=runLevelData$Protein, INTENSITY=runLevelData$LogIntensities)
runlevel.wide <- spread(data = runlevel.long, key = RUN, value = INTENSITY)
return(runlevel.wide)
}
quant.runLevel=getRunLevelQuant(processed.quant$RunlevelData)
colnames(quant.runLevel)[1] = "accession"

quant.runLevel$accession<-as.character(quant.runLevel$accession)

for (col_nr in seq(from=2, to=length(colnames(quant.runLevel))))
{
  colnames(quant.runLevel)[col_nr]=(paste0("protein_abundance_assay[", colnames(quant.runLevel)[col_nr] , "]"))
}

# TODO: check if assays in MzTab match to runs. Also true for fractionated data?

# clear old quant values from ProteinQuantifier
PRT[,grepl( "protein_abundance_assay" , names(PRT))] = NA
PRT[,grepl( "protein_abundance_study_variable" , names(PRT))] = NA

# merge in quant.runLevel values into PRT
PRT_assay_cols <- grepl("protein_abundance_assay", names(PRT))
PRT_stdv_cols <- grepl("protein_abundance_study_variable", names(PRT))
RL_assay_cols <- grepl("protein_abundance_assay", names(quant.runLevel))

for (acc in quant.runLevel$accession)
{
  q<-which(quant.runLevel$accession==acc)

  # acc from MSstats might be a group e.g., "A;B" so 
  # we check the single leader protein in mzTab PRT$accession against both A and B
  w<-which(PRT$accession %in% strsplit(acc, ";", fixed=TRUE)[[1]])

  if (length(w) == 0) 
  { 
    # TODO: check why not all summarized protein accessions are in the mzTab. Minimum number of peptides/features different?
    print(paste("Warning: ", acc, " not in mzTab but reported by MSstats"))
  }
  else
  {
    PRT[w, PRT_assay_cols] <- quant.runLevel[q, RL_assay_cols]
    PRT[w, PRT_stdv_cols] <- quant.runLevel[q, RL_assay_cols] # we currently store same data in stdv and assay column
  }
}

write.table(MTD, mzTab_output, sep = "\t", quote=FALSE, row.names = FALSE, na = "null")
write("\n",file=mzTab_output,append=TRUE)
write.table(PRT, mzTab_output, sep = "\t", quote=FALSE, row.names = FALSE, append=TRUE, na = "null")
write("\n",file=mzTab_output,append=TRUE)
write.table(PEP, mzTab_output, sep = "\t", quote=FALSE, row.names = FALSE, append=TRUE, na = "null")
write("\n",file=mzTab_output,append=TRUE)
write.table(PSM, mzTab_output, sep = "\t", quote=FALSE, row.names = FALSE, append=TRUE, na = "null")
