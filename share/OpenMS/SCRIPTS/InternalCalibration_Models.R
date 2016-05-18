library("ggplot2")
library("reshape2")

file.table.in = commandArgs(TRUE)[1] ## file.table.in = "models.csv"
file.plot.out = commandArgs(TRUE)[2] ## file.plot.out = "output.png"

d = read.csv(file.table.in, check.names = FALSE, comment.char = "#", strip.white = TRUE)
model_count = sum(d$source == "local")

dm = melt(d, id.vars = c("RT", "source"))
head(dm)

#options(device = "pdf")
#dev.new(filename = file.plot.out, file = file.plot.out)
png(filename = file.plot.out)

pl = ggplot(dm) + 
      geom_point(aes(x = RT, y=value, col=source)) +
      ggtitle(paste("Model coefficients over time\n", model_count, "model(s)", collapse="")) +
      xlab("RT [sec]") +
      facet_grid(  variable ~ ., scales="free_y")

print(pl)

dev.off()


