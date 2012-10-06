ls ## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
path<-commandArgs(TRUE)[2]
######
###TIC
######
a<-read.delim(file)
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