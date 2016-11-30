## This is an exemplary R-Script which can be used in conjunction with TOPP:GenericWrapper (of type: RScript_General)
## In this mode, the GenericWrapper provides four 'in' and six 'out' slots (four single files, two lists), which the user can couple to in/out files as desired
## Slots can be empty, and depending on who is invoking this script, you should not rely on
## argument strings being present (even empty) or not.
## e.g. a user may write
## ...... -out3 "" ...
## or leave it out completely.

## grabbing command args
## you might want to use a dedicated R package to do this
## The script will be invoked like this when used with GenericWrapper, where <inX> and <outX> might be missing completely:
## <thisScript> -in1 <in1> -in2 <in2> -in3 <in3> -in4 <in4> -out1 <out1> -out2 <out2> -out3 <out3> -out4 <out4> -outlist1 <outlist1> -outlist2 <outlist2>
argv = commandArgs(TRUE)
#argv = c("-in1", "bla", "-in3", "-out1", "o1", "-out3", "") ## internal debug, worst combination of arguments.. and we should be able to deal with it
cat("Arguments passed are:")
cat(argv)

cat("\n\nLooking at parameters now ... \n\n")

## sanity check for input. This script (arbitrarily demands that the first input file (in1) is provided plus an optional output (out1) 
##  while assuming that the outer GenericWrapper node provides up to four inputs plus six outputs)
## everything that starts with a "-" is assumed to be a parameter name (not a value)
idx_in1 = which(argv == "-in1") + 1
if (length(idx_in1)!=1 | is.na(argv[idx_in1]) | nchar(argv[idx_in1])==0 | substr(argv[idx_in1],1,1)=="-")
{
  stop("This script requires one input file for slot 'in1' for arbitrary reasons. The value must not start with '-'\n",
       "Usage:", "<thisScript> -in1 <in1> -in2 <list> [[-in3 <ignored> -in4 <ignored>] -out1 <optional> -out2 <ignored> -out3 <ignored> -out4 <ignored>] -outlist1 <optional> [-outlist2 <ignored>]", " \n");
}


in1 = argv[2]
cat(paste0("Argument -in1: '", in1, "'\n"))

idx_in2 = which(argv == "-in2") + 1
if (length(idx_in2)!=1 | is.na(argv[idx_in2]) | nchar(argv[idx_in2])==0 | substr(argv[idx_in2],1,1)=="-")
{
  stop("This script requires a second input in list format (in2) for arbitrary reasons. The values must not start with '-'\n",
       "Usage:", "<thisScript> -in1 <in1> -in2 <list> [[-in3 <ignored> -in4 <ignored>] -out1 <optional> -out2 <ignored> -out3 <ignored> -out4 <ignored>] -outlist1 <optional> [-outlist2 <ignored>]", " \n");
}
idx_in2_end = idx_in2 + 1
while (!(length(idx_in2_end)!=1 | is.na(argv[idx_in2_end]) | nchar(argv[idx_in2_end])==0 | substr(argv[idx_in2_end],1,1)=="-"))
{ ## consume as many files as present until a new parameter shows up
  idx_in2_end = idx_in2_end + 1
}
idx_in2_end = idx_in2_end - 1

in2 = argv[idx_in2:idx_in2_end]
cat(paste0("Argument -in2 (list): '", paste0(in2, collapse=" + "), "'\n"))


## do something with input ...
## ...



## deal with output (here we only look at -out1 and -outlist1 ...)
idx_out1 = which(argv == "-out1") + 1
if (length(idx_out1)==1 && !is.na(argv[idx_out1]) && nchar(argv[idx_out1])>0 && substr(argv[idx_out1],1,1)!="-")
{
  out1 = argv[idx_out1]
  cat(paste0("Argument -out1 provided as: '", out1, "'\n"))
  ## if the file is requested, we need to deliver
  cat(file=out1, "The R script wrote some output here...")
} else {
  cat("No output out1 requested!\n")
}

## deal with output (here we only look at -out1 ...
idx_outlist1 = which(argv == "-outlist1") + 1
if (length(idx_outlist1)==1 && !is.na(argv[idx_outlist1]) && nchar(argv[idx_outlist1])>0 && substr(argv[idx_outlist1],1,1)!="-")
{
  idx_outlist1_end = idx_outlist1 + 1
  while (!(length(idx_outlist1_end)!=1 | is.na(argv[idx_outlist1_end]) | nchar(argv[idx_outlist1_end])==0 | substr(argv[idx_outlist1_end],1,1)=="-"))
  { ## consume as many files as present until a new parameter shows up
    idx_outlist1_end = idx_outlist1_end + 1
  }
  idx_outlist1_end = idx_outlist1_end - 1
  outlist1 = argv[idx_outlist1:idx_outlist1_end]
  cat(paste0("Argument -outlist1 provided as: '", paste0(outlist1, collapse=" + "), "'\n"))
  ## if the file is requested, we need to deliver
  for (outlist_entry in outlist1)
  {
    cat(paste0("Writing some content to : '", outlist_entry, "' ...\n"))
    cat(file=outlist_entry, "The R script wrote some output here...")
  }
} else {
  cat("No output outlist1 requested!\n")
}
