#install.packages("scales")
#install.packages("ggplot2")
## This is an R script to produce the figures that are attached to the qcML format
library("ggplot2")
library(scales)
options(warn=-1) #suppress warnings

file<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]

#file<-"/tmp/2015-10-28_171407_Cetirizin_2990_1/TOPPAS_tmp/QCWorkflow_fr/019_QCExtractor/out_csv/old1.dta"
knime.in<-read.csv(file=file,head=TRUE,sep="\t")
names(knime.in)<-c("rt","it")
knime.in$rt <- as.POSIXct(as.character(0),format="%S")+knime.in$rt
png(post)
##########################
###Injection time
##########################
ggplot(data=knime.in, aes(x=rt, y=it)) + 
  geom_point(shape=2) + 
  geom_line(y=0, colour="blue") + 
  stat_smooth(colour="red", method=loess) +
  #scale_x_datetime( breaks="10 mins",  minor_breaks="1 mins", labels=date_format("%H:%M")) + 
  xlab("RT (HH:MM)") +
  ylab("Injection time (ms)")
######################################
garbage<-dev.off()
