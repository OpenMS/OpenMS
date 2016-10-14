library("ggplot2")
library("reshape2")
library("plyr")

file.table.in = commandArgs(TRUE)[1] ## file.table.in = "residuals.csv"
file.plot.out = commandArgs(TRUE)[2] ## file.plot.out = "residuals.png"

cat(paste0("Reading file '", file.table.in, "' to plot residual masses ..."))

d = read.csv(file.table.in, comment.char = "#", strip.white = TRUE, fill = FALSE)
head(d)

dpm = melt(d, id.vars = c("RT", "mz.ref", "intensity"))
head(dpm)
#dpm$mz.ref = factor(dpm$mz.ref)
## annotate mz.ref with average intensity
dpm2 = ddply(dpm, "mz.ref", function(x) {
  x$masstrace = paste0("m/z ",
                       # we want zero/space padded masses, such that ggplot will sort them by mass automatically
                       format(round(x$mz.ref, 5), nsmall = 5, width = 9, zero.print = TRUE),
                       " ~ int ",
                       format(median(x$intensity), scientific = TRUE, digits = 2))
  return(x)
})
head(dpm2)

#getOption("device") ## RStudioGD
#options(device = "pdf")
#dev.new(filename = file.plot.out, file = file.plot.out)
png(filename = file.plot.out, width=1920)

pl = ggplot(dpm2[grep("ppm", dpm2$variable),]) + 
       geom_hline(yintercept = 0, colour="grey") +
       geom_hline(yintercept = c(-1,1), colour="grey", linetype = "dotdash") +
       facet_wrap(~ masstrace) + 
       geom_point(aes(x = RT, y=value, color=variable)) +
       guides(color=guide_legend(title=NULL)) +
       ggtitle("Calibrant's ppm error over time") +
       xlab("RT [sec]")
print(pl)

dev.off()
