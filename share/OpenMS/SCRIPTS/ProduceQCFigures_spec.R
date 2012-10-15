ls ## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

file_p<-commandArgs(TRUE)[1]
file_id<-commandArgs(TRUE)[2]
path<-paste(commandArgs(TRUE)[3],"/",sep="")
##########################
###Spectra recordings
##########################
a<-read.delim(file_p)
b<-read.delim(file_id)
png(paste(path,"spec.png",sep=""))
plot(a$RT..sec./60,a$Precursor,pch=16,,xlab="RT (min)",ylab="m/z",main=paste(strsplit(file_p,split=".")[[1]][1],"_rec_id_specs",sep=""),cex=0.3)
points(b$RT/60,b$MZ,col="red",pch=4,cex=0.3)
legend("topleft",c("recorded spectra","identified spectra"),pch=19,col=c(1,2))
######################################
dev.off()
#
#
#
#
#