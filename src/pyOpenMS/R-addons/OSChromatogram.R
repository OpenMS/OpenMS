
OSChromatogram$$set("public","getTimeArray",function() {
  return(private$$py_obj$$getTimeArray())
})

OSChromatogram$$set("public","getIntensityArray",function() {
  return(private$$py_obj$$getIntensityArray())
})

OSChromatogram$$set("public","setTimeArray",function(data) {
  if(!(is.vector(data) && all(sapply(data,is_scalar_double)))) { stop("arg data wrong type") }
  private$$py_obj$$setTimeArray(as.list(data))
  invisible()
})

OSChromatogram$$set("public","setIntensityArray",function(data) {
  if(!(is.vector(data) && all(sapply(data,is_scalar_double)))) { stop("arg data wrong type") }
  private$$py_obj$$setIntensityArray(as.list(data))
  invisible()
})

