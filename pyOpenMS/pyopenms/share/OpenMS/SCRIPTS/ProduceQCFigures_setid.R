## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

file<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
######
###setid
######

a<-read.table(file=file, header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
######################################
png(post)
bxpo=list()
bxpo$names=a[,1]
a <- as.matrix(a[,-1])
a <- t(a[,c("min","Q1","Q2","Q3","max")])
bxpo$stats = a
bxp(bxpo)
######################################
dev.off()
#
#
#
#
#
