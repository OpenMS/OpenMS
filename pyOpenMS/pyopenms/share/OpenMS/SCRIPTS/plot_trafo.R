#!/usr/bin/env Rscript

library(XML)

## utility function:
"%within%" <- function(x, range) {
  (x >= range[1]) & (x <= range[2])
}

## read pairs of data points from trafoXML file:
read.pairs <- function(filename) {
  pairs <- matrix(nrow=0, ncol=2)
  pair.handler <- function(name, attrs) {
    pairs <<- rbind(pairs, as.numeric(c(attrs["from"], attrs["to"])))
  } 
  xmlEventParse(filename, list("Pair"=pair.handler))
  pairs
}

## create short, but unique names from trafoXML files:
unique.names <- function(paths) {
  stopifnot(!any(duplicated(paths)))
  paths <- sub("\\.trafoXML", "", paths, ignore.case=TRUE)
  labels <- basename(paths)
  if (!any(duplicated(labels)))
    return(labels)

  parts <- strsplit(paths, .Platform$file.sep, fixed=TRUE)
  parts <- lapply(parts, rev)
  i <- 2
  repeat {
    labels <- file.path(sapply(parts, function(p) p[i]), labels)
    if (!any(duplicated(labels)))
      return(labels)
    i <- i + 1
  }
}

## plot data points:
plot.pairs <- function(filenames, percent=90, pch=1, legend.loc="topleft",
                       legend.ncol=2) {
  filenames <- unique(filenames)
  pairs <- lapply(filenames, read.pairs)
  lens <- sapply(pairs, nrow)
  pairs <- do.call(rbind, pairs)
  diffs <- pairs[, 2] - pairs[, 1]
  diffs.range <- range(diffs)
  if (percent < 100) {
    frac <- (100 - percent) / 2 / 100
    q <- quantile(diffs, c(frac, 1 - frac))
    ## double the quantile range:
    yrange <- q + (diff(q) / 2) * c(-1, 1)
    ## ...unless the data range is smaller:
    yrange[1] <- max(yrange[1], diffs.range[1])
    yrange[2] <- min(yrange[2], diffs.range[2])

    xrange <- range(pairs[diffs %within% yrange, 1])
  }
  else {
    yrange <- xrange <- NULL
  }

  colors <- rainbow(length(filenames))
  plot(pairs[, 1], diffs, xlim=xrange, ylim=yrange, col=rep(colors, lens),
       pch=pch, main="Retention time transformation", xlab="original RT [s]",
       ylab=expression(paste(Delta, "RT [s]", sep="")))
  abline(h=0, col="grey")
  if (legend.loc != "none")
    legend(legend.loc, legend=unique.names(filenames), pch=20, col=colors,
           ncol=legend.ncol, cex=0.8)
}

## command line parameters:
opt <- data.frame(
    c("percent", "pch", "legend.loc", "legend.ncol"),
    desc=c("Percentage of data points to define (half the) visible range",
           "Plotting character", "Location of legend",
           "Number of columns for legend"),
    value=c("90", ".", "topleft", "2"), row.names=1,
    stringsAsFactors=FALSE)

params <- commandArgs(trailingOnly=TRUE)

if (length(params) < 2) {
  cat("Usage: Rscript Plot_trafoXML.R",
      paste0("[", rownames(opt), "=?]", collapse=" "),
      "in1.trafoXML [in2.trafoXML ...] out.pdf\n\n")
  cat("Generate a plot of RT transformation data.\n\n")
  cat("Input: trafoXML file(s)\n")
  cat("Output: PDF file with plot\n")
  cat("Optional parameters:\n")
  width <- max(nchar(rownames(opt)))
  cat(paste0("  ", format(rownames(opt), width=width), " ", opt$desc,
             " (default: ", opt$value, ")", collapse="\n"), "\n")
  quit("no")
}

## no R package for handling command line parameters installed by default :-(
params.split <- strsplit(params, "=", fixed=TRUE)
for (i in 1:length(params)) {
  parts <- params.split[[i]]
  if (length(parts) == 1)
    break # no "=", therefore no optional parameter
  if (!(parts[[1]] %in% rownames(opt))) {
    cat("Unknown parameter:", parts[[1]], "- ignored.\n")
    next
  }
  parts[[2]] <- sub("^['\"](.*)['\"]$", "\\1", parts[[2]]) # remove quotes
  opt[parts[[1]], "value"] <- parts[[2]]
}

filenames <- params[i:(length(params) - 1)]
outfile <- params[length(params)]
for (i in 1:nrow(opt)) {
  assign(rownames(opt)[i], opt[i, "value"])
}
percent <- as.numeric(percent)
legend.ncol <- as.numeric(legend.ncol)
if (pch %in% as.character(1:25))
  pch <- as.numeric(pch)

pdf(outfile)
plot.pairs(filenames, percent, pch, legend.loc, legend.ncol)
invisible(dev.off())

cat("Done.\n")
