ls ## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

path<-commandArgs(TRUE)[1]
######
###TIC
######
file = dir(pattern="_TIC.tsv")
a<-read.delim(dir(pattern="_TIC.tsv"))
######################################
png(paste(path,"tic.png",sep=""))
res = barplot(t(a$TIC),main=paste(strsplit(file,split="_TIC")[[1]][1],"_TIC",sep=""), xlab="RT (min)",ylab="Intensity")
time_seq = (seq(min(a$RT..sec.),max(a$RT..sec.),1)) 
time_seq<-round(time_seq/60)
t<-which(time_seq %% 1==0 & duplicated(time_seq)==F)
tmp = seq(min(res),max(res),max(res)/(length(time_seq)))
axis(1,at=tmp[t],labels=time_seq[t])
######################################
dev.off()
#
#
#
#
#
##########################
###Spectra recordings
##########################
file = dir(pattern="_PRECURSOR.tsv")
a<-read.delim(dir(pattern="_PRECURSOR.tsv"))
b<-read.delim(dir(pattern="_id.tsv"))
png(paste(path,"spec.png",sep=""))
plot(a$RT..sec./60,a$Precursor,pch=16,,xlab="RT (min)",ylab="m/z",main=paste(strsplit(file,split="_PRECURSOR")[[1]][1],"_rec_id_specs",sep=""),cex=0.3)
points(b$RT/60,b$MZ,col="red",pch=4,cex=0.3)
legend("topleft",c("recorded spectra","identified spectra"),pch=19,col=c(1,2))
######################################
dev.off()
#
#
#
#
#
##########################
###Mass accuracy
##########################
file = dir(pattern="_id.tsv")
a<-read.delim(dir(pattern="_id.tsv"))
png(paste(path,"accuracy.png",sep=""))
tmp<-(1-a$MZ/a$TheoreticalWeight)*1e6
hist(tmp,xlim=c(-10,10),breaks=10000,xlab="ppm",main=paste(strsplit(file,split="_id")[[1]][1],"_mass_accuracy",sep=""))
abline(v=median(tmp),col="red", lwd=2)
mtext(paste("median(accuracy)=",round(median(tmp),3)," ppm",sep=""))
######################################
dev.off()