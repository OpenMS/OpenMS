#install.packages("scales")
#install.packages("ggplot2")
## This is an R script to produce the figures that are attached to the qcML format
library("ggplot2")
library(scales)
options(warn=-1) #suppress warnings
##########################
###Injection time
##########################
knime.in<-read.csv(file=file,head=TRUE,sep="\t")
knime.in$rt <- as.POSIXct(as.character(0),format="%S")+knime.in$"MS:1000894_[sec]"
knime.in$it <- knime.in$"MS:1000927"
ggplot(data=knime.in, aes(x=rt, y=it)) + 
 geom_point(shape=2) + 
 geom_smooth() +
 scale_x_datetime( breaks="10 mins",  minor_breaks="1 mins", labels=date_format("%H:%M")) + 
 xlab("RT (HH:MM)") +
 ylab("Injection time (ms)")
##########################
