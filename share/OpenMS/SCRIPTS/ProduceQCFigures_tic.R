## This is an R script to produce the figures that are attached to the qcML format
library("ggplot2")
library(scales)
options(warn=-1) #suppress warnings

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
knime.in<-read.csv(file=file,head=TRUE,sep="\t")
names(knime.in)<- c("RT", "TIC")
png(post)
######################################
###TIC
######################################
knime.in$rt <- as.POSIXct(as.character(0),format="%S")+knime.in$RT
ggplot(data=knime.in, aes(x=rt, y=TIC)) + 
  geom_line() + 
  #scale_x_datetime( breaks="5 mins",  minor_breaks="1 mins", labels=date_format("%H:%M")) + 
  xlab("RT (HH:MM)") 
######################################
garbage<-dev.off()
