## This is an exemplary R-Script which can be used in conjunction with TOPP:GenericWrapper (of type: RScript_General)
## In this mode, the GenericWrapper provides four 'in' and two 'out' slots, which the user can couple to in/out files as desired
## (slots can be empty)

## grabbing command args
## you might want to use a dedicated R package to do this
## The script will be invoked like this when used with GenericWrapper:
## <thisScript> -in1 <in1> -in2 <in2> -in3 <in3> -in4 <in4> -out1 <out1> -out2 <out2>
argv = commandArgs(TRUE)
#argv = c("-in1", "bla", "-out1", "o1") ## internal debug
print("Arguments passed are:")
print(argv)

## sanity check for input. This script (arbitrarily demands that the first input file (in1) is provided plus an optional output (out1) 
##  while assuming that the outer GenericWrapper node provides up to four inputs followed by two outputs)
if (is.na(argv[2]) | nchar(argv[2])==0)
{
  stop("This script requires one input file (in1).\n",
       "All other input files (2-4) will be ignored.\n",
       "Usage:", "<thisScript> -in1 <in1> [[-in2 <ignored> -in3 <ignored> -in4 <ignored>] -out1 <optional> -out2 <ignored>]", " \n");
}

in1 = argv[2]
print(paste0("Argument in1: '", in1, "'"))
## do something with input ...
## ...


## deal with output (here we only look at -out1 ...
idx_out1 = which(argv == "-out1") + 1
if (length(idx_out1)==1 && !is.na(argv[idx_out1]) && nchar(argv[idx_out1])>0)
{
  out1 = argv[idx_out1]
  print(paste0("Argument out1 provided as: '", out1, "'"))
  ## if the file is requested, we need to deliver
  cat(file=out1, "The R script wrote some output here...")
} else {
  print("No output requested!")
}
