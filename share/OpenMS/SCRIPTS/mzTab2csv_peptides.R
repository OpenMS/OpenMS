## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
######
###TIC
######
a<-read.csv(file=file,head=TRUE,sep="\t")
######################################
#png(post)
#res = barplot(t(a$TIC), xlab="RT (min)",ylab="Intensity")
#time_seq = seq(min(a$RT),max(a$RT),1200)
#time_seq = round(time_seq)/60
#t<-which(time_seq %% 1==0 & duplicated(time_seq)==F)
#tmp = seq(min(res),max(res),max(res)/(length(time_seq)))
#axis(1,at=tmp[t],labels=time_seq[t])
######################################
#dev.off()
#
#
#
#
#
write.table(a, post, sep=",")