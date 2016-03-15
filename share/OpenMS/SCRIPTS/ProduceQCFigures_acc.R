library("ggplot2")
library(scales)
options(warn=-1) #suppress warnings

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
#file<-"/tmp/TOPPAS_out/023-QCExtractor-out_csv/old1.csv"
knime.in<-read.csv(file=file,head=TRUE,sep="\t")
names(knime.in)<- c("RT", "MZ", "Score", "PeptideSequence", "Charge", "TheoreticalWeight", "DeltaPpm")

png(post)
##########################
###Mass accuracy
##########################
ggplot(knime.in, aes(x=DeltaPpm)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.1, fill="green") + # Overlay with transparent density plot
  geom_vline(aes(xintercept=median(DeltaPpm, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1) + 
  xlim(c(-10,10)) + 
  ylab("Density")
######################################
garbage<-dev.off()
