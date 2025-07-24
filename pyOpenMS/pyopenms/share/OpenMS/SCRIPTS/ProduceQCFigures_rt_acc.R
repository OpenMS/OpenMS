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
if(nrow(knime.in) < 2){
  df <- data.frame()
  ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10)
}else{
  knime.in$rt <- as.POSIXct(as.character(0),format="%S")+knime.in$RT
  ggplot(data=knime.in, aes(x=rt , y=DeltaPpm)) + 
    geom_point(alpha=0.5) + 
    ylim(c(-10,10)) + 
    geom_line(y=0, colour="blue") + 
    stat_smooth(colour="red", method=loess, span=1/5) +
    xlab("RT (HH:MM)") 
}
######################################
garbage<-dev.off()
