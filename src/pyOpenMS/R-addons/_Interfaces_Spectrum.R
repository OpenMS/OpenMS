
# Returns a vector of double values.
Spectrum$$set("public","getMZArray",function() {
  return(private$$py_obj$$getMZArray())
} )

# Returns a vector of double values.
Spectrum$$set("public","getIntensityArray",function() {
  return(private$$py_obj$$getIntensityArray())
} )

Spectrum$$set("public","setMZArray",function(data) {
  if ( !(is.vector(data) && all(sapply(data, function(x) {is_scalar_integer(x) || is_scalar_double(x)})) ) ) { stop("arg transitions wrong type") }
  private$$py_obj$$setMZArray(as.list(data))
  invisible()
} )

Spectrum$$set("public","setIntensityArray",function(data) {
  if ( !(is.vector(data) && all(sapply(data, function(x) {is_scalar_integer(x) || is_scalar_double(x)})) ) ) { stop("arg transitions wrong type") }
  private$$py_obj$$setIntensityArray(as.list(data))
  invisible()
} )



