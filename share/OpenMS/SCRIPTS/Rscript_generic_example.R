## This is an exemplary R-Script which can be used in conjunction with TOPP:GenericWrapper (of type: RScript_General)
## In this mode, the GenericWrapper provides four 'in' and four 'out' slots, which the user can couple to in/out files as desired
## Slots can be empty, and depending on who is invoking this script, you should not rely on
## argument strings being present (even empty) or not.
## e.g. a user may write
## ...... -out3 "" ...
## or leave it out completely.

## grabbing command args
## you might want to use a dedicated R package to do this
## The script will be invoked like this when used with GenericWrapper, where <inX> and <outX> might be missing completely:
## <thisScript> -in1 <in1> -in2 <in2> -in3 <in3> -in4 <in4> -out1 <out1> -out2 <out2> -out3 <out3> -out4 <out4>
argv = commandArgs(TRUE)
#argv = c("-in1", "bla", "-in3", "-out1", "o1", "-out3", "") ## internal debug, worst combination of arguments.. and we should be able to deal with it
print("Arguments passed are:")
print(argv)

## sanity check for input. This script (arbitrarily demands that the first input file (in1) is provided plus an optional output (out1) 
##  while assuming that the outer GenericWrapper node provides up to four inputs plus four outputs)
## everything that starts with a "-" is assumed to be a
idx_in1 = which(argv == "-in1") + 1
if (length(idx_in1)!=1 | is.na(argv[idx_in1]) | nchar(argv[idx_in1])==0 | substr(argv[idx_in1],1,1)=="-")
{
  stop("This script requires one input file (in1), which must not start with '-'\n",
       "All other input files (2-4) will be ignored.\n",
       "Usage:", "<thisScript> -in1 <in1> [[-in2 <ignored> -in3 <ignored> -in4 <ignored>] -out1 <optional> -out2 <ignored>]", " \n");
}

in1 = argv[2]
print(paste0("Argument in1: '", in1, "'"))
## do something with input ...
## ...


## deal with output (here we only look at -out1 ...
idx_out1 = which(argv == "-out1") + 1
if (length(idx_out1)==1 && !is.na(argv[idx_out1]) && nchar(argv[idx_out1])>0 && substr(argv[idx_out1],1,1)!="-")
{
  out1 = argv[idx_out1]
  print(paste0("Argument out1 provided as: '", out1, "'"))
  ## if the file is requested, we need to deliver
  cat(file=out1, "The R script wrote some output here...")
} else {
  print("No output requested!")
}
