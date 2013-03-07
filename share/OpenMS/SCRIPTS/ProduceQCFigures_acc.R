ls ## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
##########################
###Mass accuracy
##########################
a<-read.delim(file)
png(post)
if ('delta_ppm' %in% colnames(a)) {
	hist(a$delta_ppm,xlim=c(-10,10),breaks=seq(min(a$delta_ppm)-0.01, max(a$delta_ppm)+0.01, 0.01),xlab="ppm",main=paste("delta ppm"))
	abline(v=median(a$delta_ppm),col="red", lwd=2)
	mtext(paste("median(accuracy)=",round(median(a$delta_ppm),3)," ppm",sep=""))
} else {
	tmp<-(1-a$MZ/a$TheoreticalWeight)*1e6
	hist(tmp,xlim=c(-10,10),breaks=seq(min(tmp)-0.01, max(tmp)+0.01, 0.01),xlab="ppm",main=paste("delta ppm"))
	abline(v=median(tmp),col="red", lwd=2)
	mtext(paste("median(accuracy)=",round(median(tmp),3)," ppm",sep=""))
}
######################################
dev.off()
