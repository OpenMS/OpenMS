## This is an R script to produce the figures that are attached to the qcML format
library("ggplot2")
library(scales)
options(warn=-1) #suppress warnings

#options
options(digits=10)

file_p<-commandArgs(TRUE)[1]
file_id<-commandArgs(TRUE)[2]
post<-commandArgs(TRUE)[3]
png(post)
#file_p<-"/tmp/TOPPAS_out/024-QCExtractor-out_csv/old1.csv"
#file_id<-"/tmp/TOPPAS_out/023-QCExtractor-out_csv/old1.csv"
#precs<-read.table(file_p, header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
#ids<-read.table(file_id, header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
precs<-read.csv(file=file_p,head=TRUE,sep="\t")
ids<-read.csv(file=file_id,head=TRUE,sep="\t")
names(precs)<- c("RT", "MZ")
names(ids)<- c("RT", "MZ", "Score", "PeptideSequence", "Charge", "TheoreticalWeight", "DeltaPpm")

##########################
###IDs on rt/mz map vs precursors
##########################

spec<-cbind(precs[])
id<-cbind(ids[,1:2])
spec$color<-"is_recorded"
id$color<-"is_identified"
spec$rt<-as.POSIXct(as.character(0),format="%S")+spec$RT
id$rt<-as.POSIXct(as.character(0),format="%S")+id$RT

ggplot(spec, aes(rt, MZ, color=color)) + 
  geom_point() + 
  geom_point(data=id, aes(rt, MZ, color=color)) +
  xlab("RT (HH:MM)") 
######################################
garbage<-dev.off()

