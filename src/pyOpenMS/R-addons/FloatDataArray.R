
# Gets the raw data for the float data array
# Example usage:
# fd = FloatDataArray$$new()
# data = fd$$get_data()
FloatDataArray$$set("public","get_data",function(){
  ans <- private$$py_obj$$get_data()
  return(as.vector(ans))
}
)


# Sets the raw data for the float data array
#
# Example usage:
#
# fd = FloatDataArray$$new()
# data = as.double(1:3)
# fd$$set_data(data)
FloatDataArray$$set("public","set_data",function(data){
  if (!( is_vector(data) && is_double(data) && (is.null(ncol(data)) || is.na(ncol(data))) )) { stop(paste0("Wrong argument ",data)) }
  private$$py_obj$$set_data(npy$$asarray(data)$$astype(npy$$float32))
}
)

#' @export
`[.FloatDataArray` <- function(x,ix){
  stopifnot(R6::is.R6(x))
  if(!(isTRUE(all.equal(ix,as.integer(ix))))) { stop("index must be integer") }
  tryCatch({
    return(x$$.__enclos_env__$$private$$py_obj[as.integer(ix-1)])
           }, error = function(e) { stop(paste0("invalid index ",ix)) }
  )
}

#' @export
`[<-.FloatDataArray` <- function(x,ix,value){
  stopifnot(R6::is.R6(x))
  if(!(isTRUE(all.equal(ix,as.integer(ix))))) { stop("index must be integer") }
  if(!(is_scalar_double(value) || is_scalar_integer(value))) { stop("value must be numeric") }
    tryCatch({
    x$$.__enclos_env__$$private$$py_obj$$`__setitem__`(as.integer(ix-1), value)
    invisible(x)
           }, error = function(e) { stop(paste0("invalid index ",ix)) }
  )
}

