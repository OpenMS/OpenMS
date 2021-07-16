
# Cython signature: numpy_vector, numpy_vector get_peaks()
# Will return a list of two vectors (m/z, intensity) corresponding
# to the peaks in the MSSpectrum. Provides fast access to peaks.
MSSpectrum$$set("public","get_peaks", function() {
  res <- private$$py_obj$$get_peaks()
  res[[1]] <- as.vector(res[[1]])
  res[[2]] <- as.vector(res[[2]])
  return(res)
} )

# Cython signature: set_peaks((numpy_vector, numpy_vector))
# Takes a list of two vectors/arrays(1d) (m/z, intensity) and populates the MSSpectrum.
MSSpectrum$$set("public","set_peaks", function(peaks) {
  if(!(is_list(peaks) && length(peaks) == 2)) { stop("Input for set_peaks needs to be a list of size 2 (mz and intensity vector))") }
  if( !(is_vector(peaks[[1]]) && is_vector(peaks[[2]]) && length(peaks[[1]]) == length(peaks[[2]])) ) { stop("Length of mz and intensity vector must be equal!!") }
  if( !( is.null(dim(peaks[[1]])) || length(dim(peaks[[1]]))==1) ) { stop("Input mz array must be one dimensional") }
  if( !( is.null(dim(peaks[[2]])) || length(dim(peaks[[2]]))==1) ) { stop("Input intensity array must be one dimensional") }
  if( !((is_integer(peaks[[1]]) || is_double(peaks[[1]])) && (is_integer(peaks[[2]]) || is_double(peaks[[2]]))) ) { stop("mz and intensity vector can only be of integer/double type")}
  private$$py_obj$$set_peaks(list(as.array(peaks[[1]]),as.array(peaks[[2]])))
  invisible()
}
)

MSSpectrum$$set("public","intensityInRange", function(mzmin,mzmax) {
  if (!is_scalar_double(mzmin)) { stop("arg mzmin should be a scalar double value") }
  if (!is_scalar_double(mzmax)) { stop("arg mzmax should be a scalar double value") }
  return(private$$py_obj$$intensityInRange(mzmin,mzmax))
}
)
