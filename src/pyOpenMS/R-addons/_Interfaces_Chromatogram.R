
# Returns a vector of double values.
Chromatogram$$set("public","getTimeArray",function() {
  return(private$$py_obj$$getTimeArray())
} )

# Returns a vector of double values.
Chromatogram$$set("public","getIntensityArray",function() {
  return(private$$py_obj$$getIntensityArray())
} )

Chromatogram$$set("public","setTimeArray",function(data) {
  if ( !(is.vector(data) && all(sapply(data, function(x) {is_scalar_integer(x) || is_scalar_double(x)})) ) ) { stop("arg transitions wrong type") }
  private$$py_obj$$setTimeArray(as.list(data))
  invisible()
} )

Chromatogram$$set("public","setIntensityArray",function(data) {
  if ( !(is.vector(data) && all(sapply(data, function(x) {is_scalar_integer(x) || is_scalar_double(x)})) ) ) { stop("arg transitions wrong type") }
  private$$py_obj$$setIntensityArray(as.list(data))
  invisible()
} )
