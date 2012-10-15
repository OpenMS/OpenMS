ls ## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
path<-paste(commandArgs(TRUE)[2],"/",sep="")
##########################
###Mass accuracy
##########################
a<-read.delim(file)
png(paste(path,"accuracy.png",sep=""))
tmp<-(1-a$MZ/a$TheoreticalWeight)*1e6
hist(tmp,xlim=c(-10,10),breaks=10000,xlab="ppm",main=paste(strsplit(file,split="_id")[[1]][1],"_mass_accuracy",sep=""))
abline(v=median(tmp),col="red", lwd=2)
mtext(paste("median(accuracy)=",round(median(tmp),3)," ppm",sep=""))
######################################
dev.off()