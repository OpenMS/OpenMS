## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

##########################
###IDs on rt/mz map vs precursors
##########################
file_p<-read.table(commandArgs(TRUE)[1], header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
file_id<-read.table(commandArgs(TRUE)[2], header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
post<-commandArgs(TRUE)[3]

file_p$RT_.sec. = file_p$RT_.sec./60
file_id$RT = file_id$RT/60

spec<-cbind(file_p,col=rep(1,length(file_p$RT_.sec.)))
id<-cbind(file_id[,1:2],rep(2,length(file_id$RT)))
all<-rbind(as.matrix(spec),as.matrix(id))

png(post)
plot(all,col=all[,3], xlab="RT[sec]", ylab="mz", pch=4, cex=0.3)
legend("topleft",c("recorded spectra","identified spectra"),pch=19,col=c(1,2))
######################################
dev.off()
