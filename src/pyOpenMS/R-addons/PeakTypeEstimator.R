
PeakTypeEstimator$$set("public","estimateType",function(spec){
  if(!(is.R6(spec) && class(spec)[1]=="MSSpectrum")) { stop("arg spec wrong type") }
  return(private$$py_obj$$estimateType(r_to_py(spec)))
})