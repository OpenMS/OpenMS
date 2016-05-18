library("ggplot2")
library("reshape2")

file.table.in = commandArgs(TRUE)[1] ## file.table.in = "residuals.csv"
file.plot.out = commandArgs(TRUE)[2] ## file.plot.out = "residuals.png"

d = read.csv(file.table.in, comment.char = "#", strip.white = TRUE)
head(d)

dpm = melt(d, id.vars = c("RT", "mz.ref"))
head(dpm)
dpm$mz.ref = factor(dpm$mz.ref)

#getOption("device") ## RStudioGD
#options(device = "pdf")
#dev.new(filename = file.plot.out, file = file.plot.out)
png(filename = file.plot.out)

pl = ggplot(dpm[grep("ppm", dpm$variable),])
if (length(unique(d$mz.ref)) < 8) {
  ## only a few mass traces? --> use different point shapes
  geas = aes(x = RT, y=value, col=variable, shape = mz.ref)
}  else {
  ## too many point types (probably from PepID's) -- just use one shape
  geas = aes(x = RT, y=value, col=variable)
} 
pl = pl + geom_point(geas) +
  ggtitle("Calibrant's ppm error over time") +
  xlab("RT [sec]") 

print(pl)
dev.off()
