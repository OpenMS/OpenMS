## This is an R script to produce the figures that are attached to the qcML format
library("ggplot2")
library(scales)
options(warn=-1) #suppress warnings

#options
options(digits=10)

file_features<-commandArgs(TRUE)[1]
file_masses<-commandArgs(TRUE)[2]
post<-commandArgs(TRUE)[3]

#file_features<-"/tmp/TOPPAS_out/022-QCExtractor-out_csv/old1.csv"
knime.in<-read.csv(file=file_features,head=TRUE,sep="\t")
theoretical_masses<-read.csv(file=file_masses,head=TRUE,sep="\t")
theoretical_masses <- theoretical_masses[complete.cases(theoretical_masses),]
#names(knime.in)<- c("MZ", "RT", "Intensity", "Charge", "Quality", "FWHM", "IDs")
png(post)

##########################
###dupidu
##########################

#knime.in$rt <- as.POSIXct(as.character(0),format="%S")+knime.in$RT
knime.in$featM <- knime.in$MZ * knime.in$Charge - knime.in$Charge*1.007276467
knime.in$featN <- floor(knime.in$featM)
knime.in$featF <- knime.in$featM - floor(knime.in$featM)

png<-smoothScatter(theoretical_masses,nrpoints=0,ylab="fractional mass",xlab="nominal mass", xlim=c(0,5000),pch=4,col="blue")
png<-points(cbind(knime.in$featN,knime.in$featF),col="#88000011", pch=20,cex=1)
png<-legend("topleft",c("theoretical","experimental"),col=c("darkblue","darkred"),pch=19,bty='n')

######################################
garbage<-dev.off()

